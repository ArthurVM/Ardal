# setup.py
from setuptools import setup, Extension, find_packages
import sys
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
ext_modules = [
    Extension(
        'ardal._ardal',
        sources=['src/BitMatrix.cpp',
                 'src/RoaringMatrix.cpp',
                 'src/HybridMatrix.cpp',
                 'src/bindings.cpp',
                 'src/roaring/roaring.c'],
        include_dirs=[
            os.path.abspath("include/"),
            os.path.abspath(os.path.dirname(__file__))
        ],
        language='c++',
        extra_compile_args=['-O3', '-mavx2', '-ffast-math', '-fopenmp'],
        extra_link_args=['-O3', '-fopenmp']
    )
]

setup(
    name='ardal',
    version='0.3.0-alpha',
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
        "humanize",
        "Bio"
    ],
    ext_modules=ext_modules,
    python_requires='>=3.8'
)