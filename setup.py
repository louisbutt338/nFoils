import setuptools

setuptools.setup(
    name="nfoils",
    version="1.0.0",
    author="Louis Butt",
    author_email="LJB841@bham.ac.uk",
    description="Suite of tools for NAA analysis of neutron sources",
    long_description_content_type="text/markdown",
    url="https://github.com/louisbutt338/nFoils",
    packages=setuptools.find_packages(),
    entry_points= dict(console_scripts=['reporttocsv=scripts.ReportToCSV:main',
                                        'plotspectrum=scripts.PlotSpectrum:main',
                                        'gammalines=scripts.GammaLines:main',
                                        'cascadecorrection=scripts.CascadeCorrection:main']), # Add command line tools here
    install_requires=['os','json','math','numpy', 'actigamma', 'datetime', 'matplotlib', 'pandas', 'scipy','dataclasses'], # Add dependencies here
    classifiers=["Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
)