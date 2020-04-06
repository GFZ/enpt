#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz-potsdam.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz-potsdam.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz-potsdam.de)
#
# This software was developed within the context of the EnMAP project supported
# by the DLR Space Administration with funds of the German Federal Ministry of
# Economic Affairs and Energy (on the basis of a decision by the German Bundestag:
# 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. Please note the following exception: `EnPT` depends on tqdm, which
# is distributed under the Mozilla Public Licence (MPL) v2.0 except for the files
# "tqdm/_tqdm.py", "setup.py", "README.rst", "MANIFEST.in" and ".gitignore".
# Details can be found here: https://github.com/tqdm/tqdm/blob/master/LICENCE.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

__author__ = 'Daniel Scheffler'


with open('README.rst') as readme_file:
    readme = readme_file.read()

version = {}
with open("enpt/version.py", encoding='utf-8') as version_file:
    exec(version_file.read(), version)

requirements = [  # put package requirements here
    'numpy', 'pandas', 'scipy', 'geoarray>=0.8.11', 'py_tools_ds>=0.14.23', 'arosics>=0.9.2', 'sensormapgeo>=0.2.1',
    'cerberus', 'jsmin', 'matplotlib', 'tqdm', 'utm', 'lxml', 'numpy-indexed', 'mvgavg',
    'sicor @ git+https://gitext.gfz-potsdam.de/EnMAP/sicor.git#egg=sicor'
]

setup_requirements = ['setuptools-git']  # needed for package_data version controlled by GIT

test_requirements = ['coverage', 'nose', 'nose-htmloutput', 'rednose']

setup(
    author="Karl Segl, Daniel Scheffler, Niklas Bohn, Stéphane Guillaso",
    author_email="segl@gfz-potsdam.de, danschef@gfz-potsdam.de, nbohn@gfz-potsdam.de, stephane.guillaso@gfz-potsdam.de",
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="EnMAP Processing Tool",
    keywords=['EnPT', 'EnMAP', 'EnMAP-Box', 'hyperspectral', 'remote sensing', 'satellite', 'processing chain'],
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme,
    name='enpt',
    package_dir={'enpt': 'enpt'},
    # NOTE: if the 'package_data' files are not under CVS or Subversion version control, we need setuptools-git here,
    #       otherwise they are not included in the PyPi upload content
    package_data={"enpt": ["resources/**/**/*"]},
    packages=find_packages(exclude=['tests*', 'examples*']),  # does not seems to work, therefore use MANIFEST.in
    scripts=['bin/enpt_cli.py'],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT',
    version=version['__version__'],
    zip_safe=False
)
