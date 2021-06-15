.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitLab issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitLab issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

EnPT could always use more documentation, whether as part of the
official EnPT docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

Commit Changes
--------------

How to
~~~~~~

Ready to contribute? Here's how to set up `enpt` for local development.

1. Fork the `enpt` repo on GitLab.
2. Clone your fork locally::

    $ git clone git@git.gfz-potsdam.de:your_name_here/enpt.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up
   your fork for local development::

    $ mkvirtualenv enpt
    $ cd enpt/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the codestyle and the software tests, including
   testing other Python versions with tox::

    $ make lint
    $ python -m unittest
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitLab::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a merge request through the GitLab website.


Sign your commits
~~~~~~~~~~~~~~~~~

Please note that our license terms only allow signed commits.
A guideline how to sign your work can be found here: https://git-scm.com/book/en/v2/Git-Tools-Signing-Your-Work

If you are using the PyCharm IDE, the `Commit changes` dialog has an option called `Sign-off commit` to
automatically sign your work.


License header
~~~~~~~~~~~~~~

If you commit new Python files, please note that they have to contain the following license header:

.. code:: bash

    # EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
    #
    # Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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



Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The merge request should include tests.
2. If the merge request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The merge request should work for Python 2.6, 2.7, 3.4, 3.5, 3.6, 3.7 and 3.8. Check
   https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/merge_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::


    $ python -m unittest tests.test_enpt

