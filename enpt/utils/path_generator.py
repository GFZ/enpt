# -*- coding: utf-8 -*-

from glob import glob
import os
from xml.etree import ElementTree

from ..model.metadata import L1B_product_props


class PathGenL1BProduct(object):
    def __init__(self, root_dir: str, detector_name: str):
        self.root_dir = root_dir
        self.detector_name = detector_name
        self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.detector_fn_suffix = L1B_product_props['fn_detector_suffix'][detector_name]
        self.xml = ElementTree.parse(self.get_path_metaxml()).getroot()

    def get_path_metaxml(self):
        """Return the path of the metadata XML file."""
        return glob(os.path.join(self.root_dir, "*_header.xml"))[0]

    def get_path_imagedata(self):
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/filename" % self.detector_label))

    def get_path_cloudmask(self):
        # FIXME filename currently not included in XML
        return glob(os.path.join(self.root_dir, "*_%s_cloudmask.tif" % self.detector_fn_suffix))[0]

    def get_path_deadpixelmap(self):
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/dead_pixel_map/filename" % self.detector_label))

    def get_path_quicklook(self):
        return os.path.join(self.root_dir, self._find_in_metaxml("%s/quicklook/filename" % self.detector_label))

    def _find_in_metaxml(self, expression):
        return self.xml.findall(expression)[0].text.replace("\n", "").strip()
