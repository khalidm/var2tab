#!/usr/bin/env python

from setuptools import setup

setup(
    name='vcf2tab',
    version='0.1',
    author='Khalid Mahmood',
    author_email='khalid.mahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['thepipeline = src.main:main']
    },
    url='https://github.com/khalidm/vcf2tab',
    license='LICENSE.txt',
    description='vcf2tab is simple tool for converting vcf file to a tsv file.',
    long_description=open('README.md').read(),
    install_requires=[
        "PyVCF == 0.6.7"
    ],
)
