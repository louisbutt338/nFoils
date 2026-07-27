""" setup the bFoils package
"""
import setuptools

setuptools.setup(
    name="bfoils",
    version="1.1.2",
    author="Louis Butt",
    author_email="LJB841@bham.ac.uk",
    description="bFoils suite for foil energy spectrum measurements",
    url="https://github.com/louisbutt338/bFoils",
    packages=setuptools.find_packages(),
    python_requires='>3.11',
    install_requires=['numpy',
                      'matplotlib',
                      'scipy',
                      'sandy',
                      'actigamma',
                      'emcee',
                      'corner'],
    classifiers=["Development Status :: 4 - Beta",
                 "Intended Audience :: Developers",
                 "Intended Audience :: Science/Research",
                 "Natural Language :: English",
                 "Topic :: Scientific/Engineering",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"]
)
