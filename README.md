# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">

[![MIT](http://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/louisbutt338/nfoils/develop/LICENSE)
[![DOI:10.5281/zenodo.19108622](http://img.shields.io/badge/DOI-10.5281/zenodo.19108622-B31B1B.svg)](https://doi.org/10.5281/zenodo.19108622)

Python toolkit for neutron spectrometry with activation foils

Contact LJB841@bham.ac.uk for further information or if anything is behaving badly. Contributors are welcome!
Please reference this code appropriately in your work - accompanying papers to be provided ~2025~ ~2026~ 
when hell freezes over

### Neutron spectrum measurement

nFoils was developed for measuring complex fast neutron spectra with uncertainty quantification.
This requires (a) the calculation of probability distributions for the measurement data and (b) an unfolding method 
which takes these distributions, and estimates a distribution of likely neutron spectra. 
The method is a custom Bayesian parametric unfolding algorithm: this takes the parameterised spectrum 
as a prior and the input data distributions as a likelihood, and then samples from the combined posterior.

The package includes modules for extracting and calculating response functions, 
calibrating a gamma detector, calculating foil activities and reaction rates, and the parametric unfolding 
program to estimate the spectrum. Uncertainty on all inputs/outputs is considered everywhere

Extra tools for performing a flux estimation on a lithium target neutron source, 
postprocessing of ASCII gamma spec files, and validating of a model/unfolded neutron spectrum by comparing 
FISPACT activity predictions with the experimental results are also included, as a special treat just for you

## Installation

First set up a virtual environment with prerequisites installed \*, 
then install a development version of the package \** with the following commands:
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

\** development version encouraged, so you can adapt the code to your heart's content 

## How to use

We suggest testing the example scripts in */examples/* and using them as a starting point for your own measurement. 
Key example folders for a basic foil measurement are:

- *extract_nuclear_data* for getting distributions for the foil response functions from nuclear data, 
using [SANDY](https://github.com/luca-fiorito-11/sandy/) and [NJOY2016](https://github.com/njoy/NJOY2016)
- *calculate_foil_activities* for measuring foil activities and uncertainties from gamma spectrum peak
measurements, using [actigamma](https://github.com/fispact/actigamma)
- *unfold_neutron_spectrum* for probabilistic unfolding of a neutron spectrum, using 
[emcee](https://github.com/dfm/emcee) and [corner](https://github.com/dfm/corner.py)

## Case studies

- Analysis scripts for the initial characterisation of a Birmingham proton-lithium neutron source are available in 
*/materials/proton_lithium_march2024/*, as an appendix to a manuscript currently in review (link to be provided soon)
- Slides on the unfolding algorithm, presented at OSSFE2026, are also available in */materials/*

*Development was partly supported by an agreement between the University of Birmingham and UKAEA on a Joint Research Laboratory for 
Fusion Environment Impact on Materials, part-funded by EPSRC (grant number EP/W006839/1). It was also supported 
by a contract (contract number 14455) awarded by UKAEA to Develop a small solid lithium ceramic breeder with in-line tritium detection 
capability for calibrated neutron sources. LB is grateful for a PhD scholarship awarded by the University of Birmingham and UKAEA.*
