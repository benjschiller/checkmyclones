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
	if not float(sys.version[:3])>=2.6:
		sys.stderr.write("CRITICAL: Python version must greater than or equal to 2.6! python 2.7.1 is recommended!\n")
		sys.exit(1)
	setup(name='checkmyclones',
	      version='0.0.1',
	      description="""Provides tools to check Sanger sequencing results
	      (or any plain-text or FASTQ files) against a set of reference
	      sequences, provided in any reasonable format (including coordinates)""",
	      author='Benjamin Schiller',
	      author_email='benjamin.schiller@ucsf.edu',
	      install_requires = ['cogent'],
	      packages = ['checkmyclones'],
	      package_dir = {'': 'src'},
  	      classifiers = [
				'Development Status :: 3 - Alpha',
				'License :: Freely Distributable',
				'Intended Audience :: Developers',
				'Intended Audience :: Science/Research',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: Microsoft :: Windows',
				'Operating System :: POSIX',
				'Programming Language :: Python :: 2.6',
				'Programming Language :: Python :: 2.7',
				]
	      )
	
if __name__ == '__main__':
	main()
