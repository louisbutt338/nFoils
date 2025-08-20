# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">
*Copyright (C) University of Birmingham - All Rights Reserved*

Python toolkit for the measurement of neutron spectra and validation of nuclear data with activation foils. Still in development - get in touch if anything is misbehaving

Dependent on the SANDY package for nuclear data extraction, actigamma package for radioactive isotope properties, plus numpy/scipy/matplotlib for the usual. Also includes some postprocessing and preprocessing codes for FISPACT, OpenMC, MCNP, SPECTRA-UF

Reference to be provided. Please contact LJB841@bham.ac.uk if you would like to use the package

## installing nFoils

Ensure you are set up in your chosen python virtual environment with pip and setuptools

to install non-editable versions of the python modules (feel free to install an editable version and make your own extra functions):
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

