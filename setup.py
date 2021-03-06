#!/usr/bin/env python

from distutils.core import setup

setup(name='ortholotree',
      version='0.5',
      description='search and evaluate gene orthologs using hmm and phylogeny',
      author='Peter Oxley',
      author_email='oxpeter+git@gmail.com',
      url='https://github.com/oxpeter/ortholotree',
      package_data={'ortholotree': ['data/*.txt'],
                    'ortholotree': ['data/*.pep'],
                    'ortholotree': ['data/*.cfg'],
                    'testcode': ['*.txt']},
      packages=['orthomods', 'testcode', 'ortholotree'],
      py_modules=[  'config',
                    'ortholotree',
                  ],
      requires=['argparse',
                'Bio',
                'matplotlib',
                'numpy' ]
     )