# -*- coding: utf-8 -*-
"""EnPT path generator module for generating file paths for all kinds of EnMAP images."""

from glob import glob
import os
from xml.etree import ElementTree

from ..model.metadata import L1B_product_props


class PathGenL1BProduct(object):
    """Path generator class for generating file pathes corresponding to the EnMAP L1B product."""

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
    """Returns the path of the options json file needed for atmospheric correction."""
    from sicor import options
    path_ac = os.path.join(os.path.dirname(options.__file__), 'enmap_options.json')

    return path_ac


def get_path_ac_aerosol_table():
    import sicor
    return os.path.join(sicor.__path__[0], 'tables',
                        'linear_atm_functions_ncwv_5_npre_4_ncoz_2_ntmp_2_wvl_350.0_2550.0_1.00_pca.h5')


def get_path_ac_ch4_table():
    import sicor
    return os.path.join(sicor.__path__[0], 'tables',
                        'linear_atm_functions_ncwv_4_npre_2_ncoz_2_ntmp_1_nch4_4_wvl_350.0_2550.0_1.00_pca.h5')
