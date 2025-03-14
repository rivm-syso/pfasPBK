# PFAS PBK model

This repository contains the source code of a PFAS PBK model implementation based on the PFOA and PFOS PBK model by [Loccisano et al.](https://www.sciencedirect.com/science/article/abs/pii/S0273230010002242?via%3Dihub]) and the subsequently modified version thereof by [Westerhout et al.](https://doi.org/10.1093/toxsci/kfae006). The model is implemented in Antimony and available as a FAIR PBK annotated SBML.

# Model implementation

The model implementation files can be found in the model folder. The file [PBK_PFAS.ant](model/PBK_PFAS.ant) contains the Antimony implementation. The file [PBK_PFAS.csv](model/PBK_PFAS.csv) contains the unit specifications and the model element annotations according the FAIR PBK standard, using the harmonized terminology of the [PBPKO ontology](https://github.com/InSilicoVida-Research-Lab/pbpko).

# Running the model 

The Jupyter notebook [test_dosing.ipynb](notebooks/test_dosing.ipynb) demonstrates how this model can be used in simulations. To run this notebook, you need Python with Jupyter Notebook and the python packages listed in the [requirements](requirements.txt) file.

Install the required python packages using the command:

```
pip install -r requirements.txt
```
