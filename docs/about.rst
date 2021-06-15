=====
About
=====

.. image:: img/EnPT_Logo_clipped.png
   :width: 150px
   :alt: EnPT Logo

The EnPT Python package is an automated pre-processing pipeline for the new EnMAP hyperspectral satellite data.
It provides free and open-source features for users to transform EnMAP Level-1B data themselves to Level-2A products.
The package has been developed at the German Research Centre for Geosciences Potsdam (GFZ).

Feature overview
----------------

* read EnMAP Level-1B input data
* radiometric conversion to top-of-atmosphere radiance
* dead pixel correction
* atmospheric correction (based on SICOR_ for land and `ACwater Polymer`_ via Polymer_ for water surfaces)
* detection and correction of geometric misregistrations compared to user provided spatial reference (based on AROSICS_)
* orthorectification
* write EnMAP Level-2 output data

.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
.. _AROSICS: https://git.gfz-potsdam.de/danschef/arosics
.. _`ACwater Polymer`: https://gitlab.awi.de/phytooptics/acwater
.. _Polymer: https://forum.hygeos.com
