import os
import argparse
import subprocess
from typing import Dict
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
import yaml
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Run validation scripts arguments.")
    parser.add_argument(
        '-f',
        '--force_recompute',
        action='store_true',
        help="Force re-calculation of validation results instead of using cached results."
    )

    args = parser.parse_args()
    force_recompute = args.force_recompute

    yaml_path = Path("validation/validation_scenarios.yaml")

    with yaml_path.open("r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    scenarios = config["validation_scenarios"]

    out_dir = 'validation/results'
    for scenario in scenarios:
        run_antimony_validation_scenarios(
            out_dir,
            scenario,
            force_recompute
        )
        run_r_validation_scenarios(
            scenario,
            force_recompute
        )
        plot_scenario_results(out_dir, scenario)

def run_antimony_validation_scenarios(
    out_dir: str,
    scenario: Dict,
    force_recompute: bool
):
    run_validaton_scenario_antimony(
        scenario = scenario,
        out_file = os.path.join(out_dir, f'results_ant_{scenario['id']}.csv'),
        force_recompute = force_recompute
    )

def run_r_validation_scenarios(
    scenario: Dict,
    force_recompute: bool
):
    # Check if out file exists or whether we want to force recalculation
    if os.path.exists(scenario['r_output_file']) and not force_recompute:
        print("Skipping R validation scenarios: results already available")
        return

    # Run R validation scenarios
    print("Running R validation scenarios")
    subprocess.run(['Rscript', scenario['r_validation_script']])

def run_validaton_scenario_antimony(
  scenario: Dict,
  out_file: str,
  force_recompute: bool
):
    # Skip if output already available and no forced recalculation
    if os.path.exists(out_file) and not force_recompute:
        print(f"Skipping Antimony validation scenario {scenario['id']}: results already available")
        return

    print(f"Running Antimony validation scenario {scenario['id']}")

    # Load the PBK model
    f_ant = scenario['model_path']
    rr_model = te.loada(f_ant)

    # Set chemical parameters
    load_parametrisation(rr_model, scenario['param_file'])

    # Set age and BW parameters
    rr_model.AgeInit = 0
    rr_model.AgeRef = 0
    rr_model.BWRef = rr_model.BWBirth

    # Set initial amounts according to scenario
    if 'initial_amounts' in scenario.keys():
        for item in scenario['initial_amounts']:
            rr_model.setInitAmount(item['target'], item['amount'])

    # Define dosing regimen
    num_days = int(scenario['duration'])

    # Set initial amounts according to scenario
    if 'dosing_events' in scenario.keys():
        for item in scenario['dosing_events']:
            target = item['target']
            eid = item['id']
            daily_intake = item['amount'] # ug/day
            until = item['until']
            rr_model.addEvent(eid, False, f"time % 1 == 0 && time < {until}", False)
            rr_model.addEventAssignment(eid, target, f"{target} + BW * {daily_intake}", False)

    # Regenerate model after setting initial amounts and adding all events
    rr_model.regenerateModel(True, True)

    # Define the output selections
    selections = ['time'] + [output['output_id_ant'] for output in scenario['comparisons']]

    # Simulate the PBPK model
    evaluation_resolution = int(scenario['evaluation_resolution'])
    results = rr_model.simulate(0, num_days, evaluation_resolution * (num_days + 1), selections)

    # Create output folder if not exists
    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    # Write results file
    df = pd.DataFrame(results, columns=selections)
    df.to_csv(out_file, index=False)

def load_parametrisation(model, filename):
    df = pd.read_csv(filename)
    df['Value'] = df['Value'].astype(float)
    for index, row in df.iterrows():
        model[str(row['Parameter'])] = row['Value']

def plot_scenario_results(out_dir: str, scenario: Dict):
    # Load data
    out_file_ant = os.path.join(out_dir, f'results_ant_{scenario['id']}.csv')
    df_out_ant = pd.read_csv(out_file_ant)
    df_out_r = pd.read_csv(scenario['r_output_file'])

    for comparison in scenario['comparisons']:
        # Extract time and output variable from output Antimony model
        time_ant = df_out_ant['time'].values
        y_ant = df_out_ant[comparison['output_id_ant']].values

        output_variable = comparison['label']

        # Extract time and output variable from output R model
        time_r = df_out_r['time'].values
        y_r = df_out_r[comparison['output_id_r']].values

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 7))

        # Plot lines
        ax.plot(time_ant, y_ant, 'r-', linewidth=1, label='Antimony model')
        ax.plot(time_r, y_r, 'b--', linewidth=1, label='R model')

        # Customize plot
        ax.set_xlabel('Time', fontsize=12, fontweight='bold')
        ax.set_ylabel(output_variable, fontsize=12, fontweight='bold')
        ax.set_title(f'Output timeseries comparison Antimony and R model: {output_variable}', fontsize=14)
        ax.legend(loc='best', fontsize=11, framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)

        # Format y-axis
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        out_file = os.path.join(out_dir, f'results_{scenario['id']}_{comparison['id']}.png')
        plt.tight_layout()
        plt.savefig(out_file)

if __name__ == "__main__":
    main()
