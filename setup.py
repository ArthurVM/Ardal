# setup.py
from setuptools import setup, Extension, find_packages
import pybind11
import sys
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

class get_pybind_include:
    def __init__(self, user=False):
        self.user = user
    def __str__(self):
        return os.path.abspath("pybind11/include")
    
ext_modules = [
    Extension(
        '_ardal',
        sources=['src/AlleleMatrix.cpp'],
        include_dirs=[
            get_pybind_include(),
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