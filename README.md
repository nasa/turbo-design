# Turbo Design 3 
This tool is a streamline turbomachinery design tool solving the radial equilibrium equations. It can be used for designing compressors and turbines. The designs can have counter rotating stages, different working fluids, and cooling. The intent of this tool is to enable added flexibility in which loss models are used. Because it's a python, it can connect with custom machine learning based loss models.

# Getting Loss Models working
Loss models need to be built. I have stored set of models on github as .pkl files. They should automatically download but depending on your python version, the pickle binaries may have issues reading. 

Another way to generate them is to navigate to *turbo-design/references/Turbines/AinleyMathieson* (KackerOkapuu, Traupel, CraigCox) and run `python build_dataset.py`. This will create the loss models and save it to the .cache folder for Linux and Mac and for windows it will save to a different .cache folder in your home directory.  

# Examples 
## Turbines
[OptTurb](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/optturb-turbine/optturb.ipynb) OptTurb is part of Paht's PhD work. It's a single stage HPT Turbine designed for Purdue's Experimental aeroThermal LAb (PETAL). It's an excellent candidate for verification because it can be easily modeled using a spreadsheet [OptTurb-SingleStage.xlsx](https://github.com/nasa/turbo-design/blob/main/examples/optturb-turbine/optturb-fixed_pressure_loss2.xlsm) 

[OptTurb-multistage](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/optturb-multistage/optturb-multistage.ipynb) Multi-stage example of OptTurb. This is based off a meanline spreadsheet model [OptTurb-MultiStage.xlsx](https://github.com/nasa/turbo-design/blob/main/examples/optturb-multistage/multistage-fixed_pressure_loss2.xlsx) 

[OptTurb-radial](https://colab.research.google.com/github/nasa/turbo-design/blob/main/examples/optturb-radial-turbine/) Radial Turbine Example. This example is not based on a meanline spreadsheet since radius is constantly changing. It still needs to be tested. It uses radial equilibrium to balance the massflow. 

## Building Turbine Loss Models from Correlations
The loss correlations below were estimated using Axial steam turbines. Correlation figures are extracted and surface fitted. Each of these tutorials shows how to create and save the correlation files. 

[Ainley Mathieson](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/AinleyMathieson/ainley_mathieson.ipynb)

[Craig Cox](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/CraigCox/craig_cox.ipynb)

[Traupel](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/Traupel/traupel.ipynb)

[KackerOkapuu](https://colab.research.google.com/github/nasa/turbo-design/blob/main/references/Turbines/KackerOkapuu/kacker_okapuu.ipynb)

Need to add Dunham-Came, Moustapha-Kacker

## Compressor

Need to add Koch & Smith, Wright & Miller


# Contributors

## Fortran Verson
| Person | Contribution/Role | Dates |
| ------ | ------ | ------ |
| Simon Chen | AXOD | - 2020 |
| Arthur Glassman | TD2 | unknown |
| Paht Juangphanich | Maintainer | 2020-2022 |


## Python Version TD3
| Person | Contribution/Role | Dates | Email |
| ------ | ----------------- | ----- | ------|
| Paht Juangphanich | Turbo Design 3 | - | paht.juangphanich@nasa.gov |
| Andress William | TD2 | Summer 2021 | bill.andress@gmail.com |

