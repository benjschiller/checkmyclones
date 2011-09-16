#!/usr/bin/env python
# Setup script
"""Description:
Setup script for checkmyclones
This package is free code.
"""

import os
import os.path
import sys
from distutils.core import setup
try: import py2exe
except ImportError: pass
try: import py2app
except ImportError: pass

def main():
	if not float(sys.version[:3])>=2.7:
		sys.stderr.write("CRITICAL: Python version must greater than or equal to 2.7! python 2.7.2 is recommended!\n")
		sys.exit(1)
	setup(name='checkmyclones',
	      version='0.0.10',
	      description="""Provides tools to check Sanger sequencing results
	      (or any plain-text or FASTQ files) against a set of reference
	      sequences, provided in any reasonable format (including coordinates)""",
	      author='Benjamin Schiller',
	      author_email='benjamin.schiller@ucsf.edu',
	      requires = ['cogent (>=1.5.0)', 'scripter (>2.9.1, <3.0)'],
	      url='http://github.com/benjschiller/checkmyclones',
	      scripts = ['scripts/checkmyclones.py'],
	      packages = ['clonechecker'],
	      package_dir = {'': 'src'},
  	      classifiers = [
				'Development Status :: 3 - Alpha',
				'License :: Freely Distributable',
				'Intended Audience :: Developers',
				'Intended Audience :: Science/Research',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: Microsoft :: Windows',
				'Operating System :: POSIX',
				'Programming Language :: Python :: 2.7',
				'Topic :: Scientific/Engineering :: Bio-Informatics'
				]
	      )
	
if __name__ == '__main__':
	main()
