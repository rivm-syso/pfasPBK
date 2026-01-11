"""
Description:
This script generates a markdown report for the annotated SBML file.

Usage:
Run the script from the command line with the following syntax:
  python scripts/create_report.py

Dependencies:
See requirements.txt (install using `pip install -r requirements.txt`).
"""

import os
from pathlib import Path
import libsbml as ls

from sbmlpbkutils import (
    PbkModelReportGenerator,
    DiagramCreator,
    NamesDisplay,
    RenderMode
)

MODELS_PATH = './model/'
REPORT_PATH = './docs/'

def create_report(models_path: str, report_path: str):

    # Ensure report path
    os.makedirs(report_path, exist_ok=True)

    for file in os.listdir(models_path):
        if file.endswith('.sbml'):
            sbml_file = os.path.join(MODELS_PATH, file)

            # Load SBML file
            document = ls.readSBML(sbml_file)

            # Create report
            sbml_basename = os.path.basename(file)
            report_file = os.path.join(report_path, Path(sbml_basename).with_suffix('.report.md'))
            generator = PbkModelReportGenerator(document)
            generator.create_md_report(report_file, math_render_mode = RenderMode.TEXT)

            # Diagram creator instance
            generator = DiagramCreator()

            # Create diagram without annotations
            output_file = os.path.join(report_path, Path(sbml_basename).with_suffix('.svg'))
            generator.create_diagram(
                document,
                output_file,
                names_display = NamesDisplay.ELEMENT_IDS,
                draw_species = True,
                draw_reaction_ids = False
            )

if __name__ == '__main__':
    create_report(MODELS_PATH, REPORT_PATH)
