# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">

Python package for the measurement of neutron spectra with activation foils. Contributors are welcome!

Contact LJB841@bham.ac.uk for further information or if anything is behaving badly. 
Please reference this code appropriately in your work - accompanying papers to be provided ~2025~2026.

### Neutron spectrum measurement

The main purpose of the package is for the probabilistic unfolding of complex fast neutron spectra 
(such as a fusion spectrum). This requires all input data should be considered as probability distributions, 
e.g. foil response functions from nuclear data libraries and reaction rate calculations from foil 
gamma spectrum measurements. It also requires an unfolding method (used to estimate the neutron spectrum) 
which takes these distributions into account. For this we use Bayesian unfolding, sampling from a posterior 
made up of the prior physics-informed guess for what the spectrum should look like parametrically, and the 
likelihood made up of the input data. 

To do this, the package includes modules for extracting and calculating response functions, 
calibrating a gamma detector, calculating foil activities and reaction rates, and the Bayesian unfolding 
program to measure the spectrum. Uncertainty on all inputs/outputs is carefully considered everywhere.

Tools for performing a flux estimation on a lithium target neutron source, 
postprocessing of ASCII gamma spec files, and validating of a model/unfolded neutron spectrum by comparing 
FISPACT activity predictions with the experimental results are also included, as a special treat.

## Installation

First set up a virtual environment with prerequisites installed \*, 
then install a development version of the package with the following commands:
```
git clone https://github.com/louisbutt338/nFoils.git
cd nFoils
pip install -r requirements.txt -e .
```

To enable nuclear data extraction, install NJOY2016 by following the installation instructions on the 
[NJOY website](https://docs.njoy21.io/install.html), 
then set the path to the NJOY2016 executable with the command:
```
export NJOY=/path/to/njoy
```

\* environment should include python>=3.11.0 pip>=25.0.0 setuptools>=64.0.0 git>=2.39.5 openblas>=0.3.30

## How to use

Key example folders are:

- *extract_nuclear_data* for getting relevant foil response functions and uncertainties, 
using [SANDY](https://github.com/luca-fiorito-11/sandy/) and [NJOY2016](https://github.com/njoy/NJOY2016)
- *calculate_foil_activities* for measuring foil activities from gamma spectrum peak
measurements, using [actigamma](https://github.com/fispact/actigamma)
- *unfold_neutron_spectrum* for Bayesian unfolding of a neutron spectrum, using 
[emcee](https://github.com/dfm/emcee) and [corner](https://github.com/dfm/corner.py)

*Development was partly supported by an agreement between the University of Birmingham and UKAEA on a Joint Research Laboratory for 
Fusion Environment Impact on Materials, part-funded by EPSRC (grant number EP/W006839/1). It was also supported 
by a contract (contract number 14455) awarded by UKAEA to Develop a small solid lithium ceramic breeder with in-line tritium detection 
capability for calibrated neutron sources. LB is grateful for a PhD scholarship awarded by the University of Birmingham and UKAEA.*
