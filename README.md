# nFoils <img src="https://www.cmosc.org/wp-content/uploads/2020/04/tinfoil.jpg" alt="Foils image" width="100" height="60">

Python toolkit for the measurement of neutron spectra and validation of nuclear data with activation foils. Still in development - get in touch if anything is misbehaving

Dependent on the SANDY package for nuclear data extraction, actigamma package for radioactive isotope properties, plus numpy/scipy/matplotlib for the usual. Also includes some postprocessing and preprocessing codes related to FISPACT and SPECTRA-UF

Please reference this repo appropriately in your work - accompanying paper to be provided late 2025. Contact LJB841@bham.ac.uk for further information

## installing nFoils

Ensure you are set up in your chosen environment with python/pip/setuptools/git 

to install an editable version of the package (recommend):
```
git clone https://github.com/louisbutt338/nFoils.git
cd nFoils
pip install -r requirements.txt -e .
```

to use the SANDY nuclear data extraction tool, NJOY2016 must also be installed
clone NJOY2016 with git and follow the installation instructions provided on the [NJOY website](https://docs.njoy21.io/install.html).
then set the path to the NJOY2016 executable
```
export NJOY=/path/to/njoy
```

*Development was partly supported by an agreement between the University of Birmingham and UKAEA on a Joint Research Laboratory for Fusion Environment Impact on Materials, part-funded by the EPSRC Energy Programme [grant number EP/W006839/1]. It was also supported by a contract (Contract Number 14455) awarded by UKAEA to Develop a small solid lithium ceramic breeder with in-line tritium detection capability for calibrated neutron sources. LB is grateful for a PhD scholarship awarded by the University of Birmingham and UKAEA.*
