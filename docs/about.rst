=====
About
=====

.. image:: img/EnPT_Logo_clipped.png
   :width: 150px
   :alt: EnPT Logo

The EnPT Python package is an automated pre-processing pipeline for the new EnMAP hyperspectral satellite data.
It provides free and open-source features to transform EnMAP Level-1B data to Level-2A. The package has been developed
at the German Research Centre for Geosciences Potsdam (GFZ) as an alternative to the original DLR processing chain.

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

.. _SICOR: https://gitext.gfz-potsdam.de/EnMAP/sicor
.. _AROSICS: https://gitext.gfz-potsdam.de/danschef/arosics
