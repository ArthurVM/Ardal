# setup.py
from setuptools import setup, Extension, find_packages
from pathlib import Path
import re
import sys
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
def read_version() -> str:
    """Resolve the package version from a single source of truth."""
    env_version = os.environ.get("ARDAL_VERSION")
    if env_version:
        return env_version

    version_pattern = re.compile(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', re.M)
    here = Path(__file__).resolve().parent
    candidates = [
        here / "src" / "_version.py",
        here.parent / "src" / "_version.py",
    ]
    for path in candidates:
        if path.is_file():
            match = version_pattern.search(path.read_text())
            if match:
                return match.group(1)

    return "0.3.0-alpha"

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
        "numpy",
        "pandas",
        "scipy",
        "pyjson",
        "humanize"
    ],
    ext_modules=ext_modules,
    python_requires='>=3.8'
)
