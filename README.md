# PFAS PBK model

[![Licence](https://img.shields.io/github/license/rivm-syso/pfasPBK)](https://github.com/rivm-syso/pfasPBK/blob/main/LICENSE)
[![Build](https://img.shields.io/github/actions/workflow/status/rivm-syso/pfasPBK/build.yml?label=build)](https://github.com/rivm-syso/pfasPBK/actions/workflows/build.yml)
[![Test](https://img.shields.io/github/actions/workflow/status/rivm-syso/pfasPBK/test.yml?label=test)](https://github.com/rivm-syso/pfasPBK/actions/workflows/test.yml)
[![Validate](https://img.shields.io/github/actions/workflow/status/rivm-syso/pfasPBK/validate.yml?label=validate)](https://github.com/rivm-syso/pfasPBK/actions/workflows/validate.yml)

This repository contains the source code of a PFAS PBK model implementation based on the PFOA and PFOS PBK model by [Loccisano et al. (2011)](https://www.sciencedirect.com/science/article/abs/pii/S0273230010002242?via%3Dihub]) and the subsequently modified version thereof by [Westerhout et al. (2024)](https://doi.org/10.1093/toxsci/kfae006). The model is implemented in Antimony and available as an annotated SBML file compliant with the FAIR PBK standard.

# Model code

The model code can be found in the model folder. The file [PBK_PFAS.ant](model/PBK_PFAS.ant) contains the Antimony implementation. The file [PBK_PFAS.csv](model/PBK_PFAS.csv) contains the unit specifications and the model element annotations according the FAIR PBK standard, using the harmonized terminology of the [PBPKO ontology](https://github.com/InSilicoVida-Research-Lab/pbpko). The file [PBK_PFAS.sbml](model/PBK_PFAS.sbml) is generated automatically from the model code and annotations file using the [SBML PBK workflow](https://github.com/jwkruisselbrink/sbml-pbk-workflow). This generated SBML file is the interoperable and reusable digital resource of this model implementation.

# Demonstration of use

The Jupyter notebook [test_dosing.ipynb](notebooks/test_dosing.ipynb) demonstrates how this model can be used in simulations.

# Running the notebooks

To run the notebooks, you need Python with Jupyter Notebook and the python packages listed in the [requirements](requirements.txt) file. Install notebook and the required python packages using the command:

```
pip install -r requirements.txt
```
