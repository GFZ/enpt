# -*- coding: utf-8 -*-

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

"""EnPT path generator module for generating file paths for all kinds of EnMAP images."""

from glob import glob
import os
from xml.etree import ElementTree

from ..model.metadata import L1B_product_props

__author__ = 'Daniel Scheffler'


# NOTE:
# paths belonging to providers L1B product are included in the *_header.xml file and read within metadata reader


class PathGenL1BProduct(object):
    """Path generator class for generating file pathes corresponding to the EnMAP L1B product."""

    # TODO update this class

    def __init__(self, root_dir: str, detector_name: str):
        """Get an instance of the EnPT L1B image path generator.

        :param root_dir:
        :param detector_name:
        """
        self.root_dir = root_dir
        assert len(os.listdir(self.root_dir)) > 0, 'Image root directory must contain files.'

        self.detector_name = detector_name
        self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.detector_fn_suffix = L1B_product_props['fn_detector_suffix'][detector_name]
        self.xml = ElementTree.parse(self.get_path_metaxml()).getroot()

    def get_path_metaxml(self):
        """Return the path of the metadata XML file."""
        return glob(os.path.join(self.root_dir, "*_header.xml"))[0]

    def get_path_data(self):
        """Return the path of the image data file."""
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/filename" % self.detector_label))

    def get_path_cloudmask(self):
        """Return the path of the cloud mask file."""
        # FIXME filename currently not included in XML
        return glob(os.path.join(self.root_dir, "*_%s_cloudmask.tif" % self.detector_fn_suffix))[0]

    def get_path_deadpixelmap(self):
        """Return the path of the dead pixel mask file."""
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/dead_pixel_map/filename" % self.detector_label))

    def get_path_quicklook(self):
        """Return the path of the quicklook file."""
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/quicklook/filename" % self.detector_label))

    def _find_in_metaxml(self, expression):
        return self.xml.findall(expression)[0].text.replace("\n", "").strip()


def get_path_ac_options() -> str:
    """Return the path of the options json file needed for atmospheric correction."""
    from sicor import options
    path_ac = os.path.join(os.path.dirname(options.__file__), 'enmap_options.json')
    # FIXME temporarily disabled because not implemented at the moment:
    # path_ac = os.path.join(os.path.dirname(options.__file__), 'sicor_enmap_user_options.json')

    return path_ac
