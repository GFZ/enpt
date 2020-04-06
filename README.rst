
.. image:: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/img/EnPT_logo_final.svg
   :width: 300px
   :alt: EnPT Logo

============================
EnPT - EnMAP Processing Tool
============================

The EnPT Python package is an automated pre-processing pipeline for the new EnMAP hyperspectral satellite data.
It provides free and open-source features to transform EnMAP Level-1B data to Level-2A. The package has been developed
at the German Research Centre for Geosciences Potsdam (GFZ) as an alternative to the processing chain of the EnMAP
Ground Segment.

Please check the documentation_ for installation and usage instructions and in depth information.

License
-------
Free software: GNU General Public License v3

All images contained in any (sub-)directory of this repository are licensed under the CC0 license which can be found
`here <https://creativecommons.org/publicdomain/zero/1.0/legalcode.txt>`__.

Feature overview
----------------

* read EnMAP Level-1B input data
* radiometric conversion to top-of-atmosphere radiance
* dead pixel correction
* atmospheric correction (based on SICOR_)
* conversion of top-of-atmosphere-radiance to top-of-atmosphere-reflectance
* detection and correction of geometric misregistrations compared to user provided spatial reference (based on AROSICS_)
* orthorectification
* write EnMAP Level-2 output data

Status
------

|badge1| |badge2| |badge3| |badge4| |badge5| |badge6| |badge7|

.. |badge1| image:: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/badges/master/pipeline.svg
    :target: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/pipelines

.. |badge2| image:: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/badges/master/coverage.svg
    :target: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/coverage/

.. |badge3| image:: https://img.shields.io/static/v1?label=Documentation&message=GitLab%20Pages&color=orange
    :target: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/doc/

.. |badge4| image:: https://img.shields.io/pypi/v/enpt.svg
    :target: https://pypi.python.org/pypi/enpt

.. |badge5| image:: https://img.shields.io/pypi/l/enpt.svg
    :target: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/-/blob/master/LICENSE

.. |badge6| image:: https://img.shields.io/pypi/pyversions/enpt.svg
    :target: https://img.shields.io/pypi/pyversions/enpt.svg

.. |badge7| image:: https://img.shields.io/pypi/dm/enpt.svg
    :target: https://pypi.python.org/pypi/enpt

See also the latest coverage_ report and the nosetests_ HTML report.

History / Changelog
-------------------

You can find the protocol of recent changes in the EnPT package
`here <https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/-/blob/master/HISTORY.rst>`__.

Credits
-------

This software was developed within the context of the EnMAP project supported by the DLR Space Administration with
funds of the German Federal Ministry of Economic Affairs and Energy (on the basis of a decision by the German
Bundestag: 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _documentation: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/doc/
.. _coverage: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/coverage/
.. _nosetests: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/nosetests_reports/nosetests.html
.. _SICOR: https://gitext.gfz-potsdam.de/EnMAP/sicor
.. _AROSICS: https://gitext.gfz-potsdam.de/danschef/arosics