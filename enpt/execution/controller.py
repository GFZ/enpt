# -*- coding: utf-8 -*-
"""EnPT process controller module."""

import os
from datetime import datetime
import tempfile
import zipfile
import shutil
import weakref
import warnings

from ..options.config import EnPTConfig
from ..io.reader import L1B_Reader
from ..model.images import EnMAPL1Product_SensorGeo


class EnPT_Controller(object):
    """Class of EnPT process controller."""

    def __init__(self, config: EnPTConfig=None, **config_kwargs: dict):
        """Initialize the Process Controller.

        :param config:          an instance of the EnPTConfig class (overrides config_kwargs)
        :param config_kwargs:   configuration parameters to be passed to EnPTConfig class
        """
        self.cfg = config or EnPTConfig(**config_kwargs)  # type: EnPTConfig

        # generate temporary directory (must be deleted at the end of the pipeline)
        self.tempDir = tempfile.mkdtemp()

        # setup a finalizer that destroys remaining data (directories, etc.) in case of unexpected exit
        self._finalizer = weakref.finalize(self, self._cleanup, self.tempDir,
                                           warn_message="Implicitly cleaning up {!r}".format(self))

        # read the L1B data
        self.L1_obj = self.read_L1B_data(self.cfg.path_l1b_enmap_image)

    def extract_zip_archive(self, path_zipfile: str) -> str:
        """Extract the given EnMAP image zip archive and return the L1B image root directory path.

        :param path_zipfile:    /path/to/zipfile.zip
        :return:                /tmp/tmpk2qp0yri/rootdir/
        """
        with zipfile.ZipFile(path_zipfile, "r") as zf:
            zf.extractall(self.tempDir)

        return os.path.join(self.tempDir, os.path.basename(path_zipfile).split('.zip')[0])

    def read_L1B_data(self, path_enmap_image: str) -> EnMAPL1Product_SensorGeo:
        """

        :param path_enmap_image:
        :return:
        """
        # input validation
        if not os.path.isdir(path_enmap_image) and \
           not (os.path.exists(path_enmap_image) and path_enmap_image.endswith('.zip')):
            raise ValueError("The parameter 'path_enmap_image' must be a directory or the path to an existing zip "
                             "archive. Received %s." % path_enmap_image)

        # extract L1B image archive if needed
        if path_enmap_image.endswith('.zip'):
            path_enmap_image = self.extract_zip_archive(path_enmap_image)
            if not os.path.isdir(path_enmap_image):
                raise NotADirectoryError(path_enmap_image)

        # run the reader
        RD = L1B_Reader(config=self.cfg)
        L1_obj = RD.read_inputdata(path_enmap_image, observation_time=datetime(2015, 12, 7, 10))

        return L1_obj

    def run_toaRad2toaRef(self):
        """Run conversion from TOA radiance to TOA reflectance."""
        # get a new instance of radiometric transformer
        from ..processors.radiometric_transform import Radiometric_Transformer
        RT = Radiometric_Transformer(self.cfg)

        # run transformation to TOARef
        self.L1_obj = RT.transform_TOARad2TOARef(self.L1_obj)

    def run_geometry_processor(self):
        pass

    def run_atmospheric_correction(self):
        """Run atmospheric correction only."""
        self.L1_obj.run_AC()

    def write_output(self):
        if self.cfg.output_dir:
            self.L1_obj.save(self.cfg.output_dir)

    def run_all_processors(self):
        """Run all processors at once."""
        try:
            if self.cfg.run_deadpix_P:
                self.L1_obj.correct_dead_pixels()
            # self.run_toaRad2toaRef()  # this is only needed for geometry processor but AC expects radiance
            self.run_geometry_processor()
            self.run_atmospheric_correction()
            self.write_output()
        finally:
            self.cleanup()

    @classmethod
    def _cleanup(cls, tempDir, warn_message):
        """Clean up implicitly (not to be called directly)."""
        if tempDir:
            shutil.rmtree(tempDir)
        warnings.warn(warn_message, ResourceWarning)

    def cleanup(self):
        """Clean up (to be called directly)."""
        if self._finalizer.detach():
            if self.tempDir:
                shutil.rmtree(self.tempDir)
