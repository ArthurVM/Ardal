# setup.py
from setuptools import setup, Extension, find_packages
import pybind11
import sys
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
ext_modules = [
    Extension(
        '_ardal',
        sources=['src/AlleleMatrix.cpp'],
        include_dirs=[
            os.path.abspath(os.path.join(os.path.dirname(__file__), "pybind11/include/")),
            os.path.abspath(os.path.dirname(__file__))
        ],
        language='c++',
        extra_compile_args=['-O3', '-march=native', '-ffast-math'],
        extra_link_args=['-O3']
    )
]

setup(
    name='ardal',
    version='0.1.0',
    author="A. V. Morris",
    long_description=long_description,
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bioinformatics"
    ],
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "pyjson",
        "humanize"
    ],
    ext_modules=ext_modules,
)