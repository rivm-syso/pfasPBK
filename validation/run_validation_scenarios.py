import os
import argparse
import subprocess
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt

PATH_VALIDATION_R_SCRIPT = 'validation/R_model/run_westerhout_2024.R'

# Note: this is the path of the file that is written by the R script
OUT_FILE_R = "validation/results/results_R_PFAS.csv"

# Definition of validation scenarios and accompanying comparisons
VALIDATION_SCENARIOS = [{
        'id': 'PFOS',
        'model_path': './model/PBK_PFAS_LT.ant',
        'param_file': './validation/params_LT_validation_PFOS.csv',
        'comparisons': [
            {
                'id': 'CArt',
                'label': 'Concentration arterial plasma PFOS (ug/L)',
                'output_id_r': 'Cart_PFOS',
                'output_id_ant': '[AArt_Plas]'
            },
            {
                'id': 'CVen',
                'label': 'Concentration venous plasma PFOS (ug/L)',
                'output_id_r': 'Cven_PFOS',
                'output_id_ant': '[AVen_Plas]'
            },
            {
                'id': 'BW',
                'label': 'Bodyweight',
                'output_id_r': 'BW',
                'output_id_ant': 'BW'
            }
        ]
    },
    {
        'id': 'PFOA',
        'model_path': './model/PBK_PFAS_LT.ant',
        'param_file': './validation/params_LT_validation_PFOA.csv',
        'comparisons': [
            {
                'id': 'CArt',
                'label': 'Concentration arterial plasma PFOA (ug/L)',
                'output_id_r': 'Cart_PFOA',
                'output_id_ant': '[AArt_Plas]'
            },
            {
                'id': 'CVen',
                'label': 'Concentration venous plasma PFOA (ug/L)',
                'output_id_r': 'Cven_PFOA',
                'output_id_ant': '[AVen_Plas]'
            }
        ]
    }
]

def main():
    parser = argparse.ArgumentParser(description="Demo of force_recompute argument")
    parser.add_argument(
        '-f',
        '--force_recompute',
        action='store_true',
        help="Force re-calculation of validation results instead of using cached results."
    )

    args = parser.parse_args()
    force_recompute = args.force_recompute

    out_dir = 'validation/results'
    run_antimony_validation_scenarios(out_dir, force_recompute)
    run_r_validation_scenarios(force_recompute)
    for scenario in VALIDATION_SCENARIOS:
        plot_scenario_results(out_dir, scenario)

def run_antimony_validation_scenarios(
    out_dir: str,
    force_recompute: bool
):
    for scenario in VALIDATION_SCENARIOS:
        run_validaton_scenario_antimony(
            scenario = scenario,
            out_file = os.path.join(out_dir, f'results_ant_{scenario['id']}.csv'),
            force_recompute = force_recompute
        )

def run_r_validation_scenarios(force_recompute: bool):
    # Check if out file exists or whether we want to force recalculation
    if os.path.exists(OUT_FILE_R) or not force_recompute:
        print("Skipping R validation scenarios: results already available")
        return

    # Run R validation scenarios
    print("Running R validation scenarios")
    subprocess.run(['Rscript', PATH_VALIDATION_R_SCRIPT])

def run_validaton_scenario_antimony(
  scenario,
  out_file: str,
  force_recompute: bool
):
    # Skip if output already available and no forced recalculation
    if os.path.exists(out_file) and not force_recompute:
        print(f"Skipping Antimony validation scenario {scenario['id']}: results already available")
        return

    print(f"Running Antimony validation scenario {scenario['id']}")

    # Define dosing regimen
    daily_intake = 0.001  # ug/day
    days_of_exposure = int(80 * 365.25)
    days_after_exposure = int(0 * 365.25)
    num_days = days_of_exposure + days_after_exposure

    # Load the PBK model
    f_ant = scenario['model_path']
    rr_model = te.loada(f_ant)

    # Make sure input is not constant and does not have boundary conditions
    input_id = 'AGut'
    rr_model.setInitAmount(input_id, 0)
    rr_model.setConstant(input_id, False)
    rr_model.setBoundary(input_id, False)

    # Set chemical parameters
    load_parametrisation(rr_model, scenario['param_file'])

    # Set age and BW parameters
    rr_model.AgeInit = 0
    rr_model.AgeRef = 0
    rr_model.BWRef = rr_model.BWBirth

    # Create a repeating daily oral dosing
    eid = f"oral_dosing"
    rr_model.addEvent(eid, False, f"time % 1 == 0 && time < {days_of_exposure}", False)
    rr_model.addEventAssignment(eid, input_id, f"{input_id} + BW * {daily_intake}", False)
    rr_model.regenerateModel(True, True)

    # Define the output selections
    selections = ['time'] + [output['output_id_ant'] for output in scenario['comparisons']]

    # Simulate the PBPK model
    results = rr_model.simulate(0, num_days, num_days + 1, selections)

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

def plot_scenario_results(out_dir: str, scenario):
    # Load data
    out_file_ant = os.path.join(out_dir, f'results_ant_{scenario['id']}.csv')
    df_out_ant = pd.read_csv(out_file_ant)
    df_out_r = pd.read_csv(OUT_FILE_R)

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
