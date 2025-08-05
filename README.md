# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">

Python toolkit for the measurement of neutron spectra and validation of nuclear data with activation foils. Some parts still in development

Dependent on the SANDY package for nuclear data extraction. Also includes some postprocessing and preprocessing codes for FISPACT, OpenMC, MCNP, SPECTRA-UF (but is not dependent on them)

Reference to be provided. Package is unlicenced, as IP belongs to the University of Birmingham: please contact me (LJB841@bham.ac.uk) if you would like to use it

## installing nFoils

to install with dependencies (once you are set up in your chosen environment):
```
git clone https://github.com/louisbutt338/nFoils.git
cd nFoils
pip install .
```

to use the SANDY nuclear data extraction tool, NJOY2016 must also be installed
clone NJOY2016 with git and follow the installation instructions provided on the [NJOY website](https://docs.njoy21.io/install.html).
then set the path to the NJOY2016 executable
```
export NJOY=/path/to/njoy
```
