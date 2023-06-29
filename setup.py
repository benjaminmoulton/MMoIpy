"""MMoIpy: A python implementation of a method for calculating the mass, center of gravity (CG), and mass moments of inertia of an aircraft system."""

from setuptools import setup
import os
import sys

setup(name = 'MMoIpy',
    version = '1.0.0',
    description = "MMoIpy: A python implementation of a method for calculating the mass, center of gravity (CG), and mass moments of inertia of an aircraft system.",
    url = 'https://github.com/benjaminmoulton/MMoIpy',
    author = 'benjaminmoulton',
    author_email = 'ben.moulton@usu.edu',
    install_requires = ['numpy>=1.18', 'pytest', 'matplotlib',],
    python_requires ='>=3.6.0',
    license = 'MIT',
    packages = ['mmoipy'],
    zip_safe = False)