# setup.py
from setuptools import setup, Extension, find_packages
import pybind11
import sys
import os

ext_modules = [
    Extension(
        '_ardal',
        sources=['src/AlleleMatrix.cpp'],
        include_dirs=[
            pybind11.get_include(),
            os.path.abspath(os.path.dirname(__file__))
        ],
        language='c++',
        extra_compile_args=['-O3', '-march=native',  '-ffast-math'],
        extra_link_args=['-O3']
    )
]

setup(
    name='ardal',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "pyjson",
        "humanize"
    ],
    setup_requires=['pybind11'],
    ext_modules=ext_modules,
)