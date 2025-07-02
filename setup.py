import setuptools

setuptools.setup(
    name="nfoils",
    version="1.0.0",
    author="Louis Butt",
    author_email="LJB841@bham.ac.uk",
    description="Suite of tools for dosimetry foil characterisation of fast neutron sources",
    url="https://github.com/louisbutt338/nFoils",
    packages=setuptools.find_packages(),
    #entry_points= dict(console_scripts=['sum_asciis=scripts.sum_asciis:main']), # Add command line tools here
    install_requires=['numpy', 'matplotlib', 'scipy','sandy', 'actigamma'], # Add dependencies here
    classifiers=["Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
)