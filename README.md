# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">

Python toolkit for the characterisation of neutron spectra and validation of nuclear data with activation foils. In development - contact 
LJB841@bham.ac.uk for further information or if anything is misbehaving

Dependent on [SANDY](https://github.com/luca-fiorito-11/sandy/) and [NJOY2016](https://github.com/njoy/NJOY2016) for nuclear data extraction, 
and [actigamma](https://github.com/fispact/actigamma) for radioactive isotope data. All worth checking out if you have got this far

Please reference this repo appropriately in your work - accompanying paper to be provided late 2025

## installing nFoils

first, set up a virtual environment for the package with installation dependencies:

for conda,
```
conda create --name nfoils_env python pip setuptools openblas
```

second, install the package:

for non-editable version,
```
git clone https://github.com/louisbutt338/nFoils.git
cd nFoils
pip install -r requirements.txt .
```
(use ```pip install -r requirements.txt -e .``` to install an editable version)

third, install NJOY2016 to use the reaction.py module for nuclear data extraction:

clone NJOY2016 with git and follow the installation instructions provided on the [NJOY website](https://docs.njoy21.io/install.html)

then set the path to the NJOY2016 executable,
```
export NJOY=/path/to/njoy
```

\*tested for python>=3.11.0 pip>=25.0.0 setuptools>=64.0.0 git>=2.39.5 openblas>=0.3.30

*Development was partly supported by an agreement between the University of Birmingham and UKAEA on a Joint Research Laboratory for 
Fusion Environment Impact on Materials, part-funded by EPSRC (grant number EP/W006839/1). It was also supported 
by a contract (contract number 14455) awarded by UKAEA to Develop a small solid lithium ceramic breeder with in-line tritium detection 
capability for calibrated neutron sources. LB is grateful for a PhD scholarship awarded by the University of Birmingham and UKAEA.*
