"""
Description:
This script creates the annotated SBML file from the Antimony PBK model implementation.

Usage:
Run the script from the command line with the following syntax:
  python models/compile_models.py

Dependencies:
See requirements.txt (install using `pip install -r requirements.txt`).
"""

import os
import uuid
import tellurium as te
import libsbml as ls
import logging
from pathlib import Path

model_path = './model/'

from sbmlpbkutils import PbkModelValidator
from sbmlpbkutils import AnnotationsTemplateGenerator
from sbmlpbkutils import PbkModelAnnotator

def create_file_logger(logfile: str) -> logging.Logger:
    logger = logging.getLogger(uuid.uuid4().hex)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(logfile, 'w+')
    formatter = logging.Formatter('[%(levelname)s] - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

citation_file = './CITATION.cff'

for file in os.listdir(model_path):
    if file.endswith('.ant'):
        ant_file = os.path.join(model_path, file)
        sbml_file = Path(ant_file).with_suffix('.sbml')

        print(f"Creating SBML file [{sbml_file}] from Antimony file [{ant_file}].")
        r = te.loada(ant_file)
        r.exportToSBML(sbml_file)

        document = ls.readSBML(sbml_file)

        annotations_file = Path(sbml_file).with_suffix('.annotations.csv')
        if not os.path.exists(annotations_file):
            # create annotations (csv) file if it does not exist
            print(f"Annotations file not found: creating annotations file [{annotations_file}].")
            model = document.getModel()
            annotations_template_generator = AnnotationsTemplateGenerator()
            annotations = annotations_template_generator.generate(model)
            annotations.to_csv(annotations_file, index=False)

        annotations_log_file = Path(sbml_file).with_suffix('.annotations.log')
        annotated_sbml_file = Path(sbml_file).with_suffix('.sbml')
        annotator = PbkModelAnnotator()
        logger = create_file_logger(annotations_log_file)
        print(f"Creating annotated SBML file [{annotated_sbml_file}] from SBML file [{sbml_file}] with annotations file [{annotations_file}].")
        annotator.annotate(
            document,
            annotations_file,
            citation_file,
            logger
        )
        ls.writeSBML(document, str(annotated_sbml_file))

        validation_log_file = Path(sbml_file).with_suffix('.validation.log')
        validator = PbkModelValidator()
        logger = create_file_logger(validation_log_file)
        validator.validate(annotated_sbml_file, logger)
