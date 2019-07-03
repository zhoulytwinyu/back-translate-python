#!/usr/bin/env python

from distutils.core import setup

setup(name='backtranslate',
      version='1.0',
      description='efficient peptide sequence back translation algorithm',
      author='Lingyu Zhou',
      author_email='zhoulytwin@gmail.com',
      install_requires=['docopt>=0.6.2'],
      packages=['backtranslate'],
      package_data={'backtranslate':["data/standard_codon.tsv"]},
      scripts=['scripts/back-translate.py'],
      include_package_data=True
     )
