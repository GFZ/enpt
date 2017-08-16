#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='enpt',
    version='0.1.0',
    description="EnMAP PT",
    long_description=readme + '\n\n' + history,
    author="Karl Segl",
    author_email='segl@gfz-potsdam.de',
    url='https://github.com/audreyr/enpt',
    packages=[
        'enpt',
    ],
    package_dir={'enpt':
                 'enpt'},
    include_package_data=True,
    install_requires=requirements,
    license="GNU General Public License v3",
    zip_safe=False,
    scripts = [], # TODO
    data = [], # TODO
    keywords='enpt',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
