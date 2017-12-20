#!/usr/bin/env python
import setuptools
from distutils.core import setup

setup(name='EasyGromacs',
      version='0.0.6',
      description='EasyGromacs Package',
      author='Mustafa Tekpinar',
      author_email='tekpinar@buffalo.edu',
      url='www.tekpinarlab.org',
      license='LGPLv2',
      py_modules=['easyGromacs'],
      install_requires=['wxPython'],
      classifiers=[
              # How mature is this project? Common values are
              #   3 - Alpha
              #   4 - Beta
              #   5 - Production/Stable
              'Development Status :: 3 - Alpha',

              # Indicate who your project is intended for
              'Intended Audience :: Science/Research',
              'Topic :: Scientific/Engineering :: Bio-Informatics',

              # Pick your license as you wish (should match "license" above)
               'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',

              # Specify the Python versions you support here. In particular, ensure
              # that you indicate whether you support Python 2, Python 3 or both.
              'Programming Language :: Python :: 2',
              'Programming Language :: Python :: 2.6',
              'Programming Language :: Python :: 2.7',],
      keywords='molecular dynamics gromacs ',
)
