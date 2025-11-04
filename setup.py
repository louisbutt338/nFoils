""" setup the nFoils package
"""
import setuptools

setuptools.setup(
    name="nfoils",
    version="1.0.0",
    author="Louis Butt",
    author_email="LJB841@bham.ac.uk",
    description="nFoils suite for foil neutron measurements",
    url="https://github.com/louisbutt338/nFoils",
    packages=setuptools.find_packages(),
    python_requires='>3.11',
    install_requires=['numpy',
                      'matplotlib',
                      'scipy',
                      'sandy',
                      'actigamma'],
    classifiers=["Development Status :: 4 - Beta",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "Natural Language :: English",
                 "Topic :: Scientific/Engineering",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"]
)
