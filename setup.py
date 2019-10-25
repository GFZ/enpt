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
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
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

with open('HISTORY.rst') as history_file:
    history = history_file.read()

version = {}
with open("enpt/version.py") as version_file:
    exec(version_file.read(), version)

requirements = [  # put package requirements here
    'numpy', 'pandas', 'scipy', 'geoarray>=0.8.11', 'py_tools_ds>=0.14.23', 'arosics>=0.9.2', 'cerberus', 'jsmin',
    'matplotlib', 'tqdm', 'utm', 'lxml', 'numpy-indexed'
    # 'sicor', # pip install git+https://gitext.gfz-potsdam.de/EnMAP/sicor.git
]

test_requirements = ['coverage', 'nose', 'nose-htmloutput', 'rednose']

setup(
    name='enpt',
    version=version['__version__'],
    description="EnMAP Processing Tools",
    long_description=readme + '\n\n' + history,
    author=["Karl Segl", "Daniel Scheffler", "Niklas Bohn", "Stéphane Guillaso"],
    author_email=['segl@gfz-potsdam.de', 'danschef@gfz-potsdam.de', 'nbohn@gfz-potsdam.de',
                  'stephane.guillaso@gfz-potsdam.de'],
    url='https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT',
    packages=find_packages(exclude=['tests*']),
    package_dir={'enpt':
                 'enpt'},
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    zip_safe=False,
    scripts=['bin/enpt_cli.py'],
    data=[],  # TODO
    keywords='enpt',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
