.. image:: docs/img/EnPT_Logo_clipped.png
   :width: 100px
   :alt: EnPT Logo

=============================
EnPT - EnMAP Processing Tools
=============================

The EnPT Python package is an automated pre-processing pipeline for the new EnMAP hyperspectral satellite data.
It provides free and open-source features to transform EnMAP Level-1B data to Level-2A. The package has been developed
at the German Research Centre for Geosciences Potsdam (GFZ) as an alternative to the original DLR processing chain.

Please check the documentation_ for usage and in depth information.

License
-------
Free software: GNU General Public License v3

Features
--------

* read EnMAP Level-1B input data
* dead pixel correction
* import of an overlapping digital elevation model
* radiometric conversion to top-of-atmosphere radiance
* atmospheric correction (based on SICOR_)
* detection and correction of geometric misregistrations compared to user provided spatial reference (based on AROSICS_)
* orthorectification
* write EnMAP Level-2 output data

Status
------

|badge1| |badge2|

.. |badge1| image:: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/badges/master/build.svg

.. |badge2| image:: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/badges/master/coverage.svg

See also the latest coverage_ report and the nosetests_ HTML report.


Credits
-------

The development of the gms_preprocessing package was funded by the German Federal Ministry of Education and Research
(BMBF).

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _documentation: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/doc/
.. _coverage: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/coverage/
.. _nosetests: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/nosetests_reports/nosetests.html
.. _SICOR: https://gitext.gfz-potsdam.de/EnMAP/sicor
.. _AROSICS: https://gitext.gfz-potsdam.de/danschef/arosics