# setup.py
from setuptools import setup, Extension, find_packages
from pathlib import Path
import re
import sys
import os

with open("README.md", "r") as fh:
    long_description = fh.read()


def read_version() -> str:
    version_file = Path(__file__).resolve().parent / "ardal" / "_version.py"
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', version_file.read_text(), re.M)
    if not match:
        raise RuntimeError("Could not find __version__ in ardal/_version.py")
    return match.group(1)
    
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
    version=read_version(),
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
        "numpy>=2.3.3",
        "pandas>=2.3.2",
        "scipy>=1.16.2",
        "pyjson>=1.4.1",
        "humanize>=4.13.0",
        "biopython>=1.85"
    ],
    ext_modules=ext_modules,
    python_requires='>=3.8'
)
