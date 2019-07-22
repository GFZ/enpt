# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tools - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
#
# This software was developed within the context of the EnMAP project supported
# by the DLR Space Administration with funds of the German Federal Ministry of
# Economic Affairs and Energy (on the basis of a decision by the German Bundestag:
# 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""EnPT process controller module."""

import os
import tempfile
import zipfile
import shutil
import weakref
import warnings

from ..options.config import EnPTConfig
from ..io.reader import L1B_Reader
from ..model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo  # noqa: F401


class EnPT_Controller(object):
    """Class of EnPT process controller."""

    def __init__(self, config: EnPTConfig = None, **config_kwargs: dict):
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

        # defaults
        self.L1_obj = None  # type: EnMAPL1Product_SensorGeo
        self.L2_obj = None  # type: EnMAPL2Product_MapGeo

    def extract_zip_archive(self, path_zipfile: str, subdir: str = '') -> str:
        """Extract the given EnMAP image zip archive and return the L1B image root directory path.

        :param path_zipfile:    /path/to/zipfile.zip
        :param subdir:          subdirectory name to be created within temporary directory
        :return:                /tmp/tmpk2qp0yri/rootdir/
        """
        outdir = os.path.join(self.tempDir, subdir) if self.cfg.is_dlr_dataformat else self.tempDir

        with zipfile.ZipFile(path_zipfile, "r") as zf:
            zf.extractall(outdir)

        if not os.path.isdir(outdir):
            raise NotADirectoryError(outdir)

        if self.cfg.is_dlr_dataformat:
            return outdir
        else:
            return os.path.join(self.tempDir, os.path.basename(path_zipfile).split('.zip')[0])

    def read_L1B_data(self) -> None:
        """Read the provider L1B data given in config and return an EnMAP image object."""
        path_enmap_image = self.cfg.path_l1b_enmap_image
        path_enmap_image_gapfill = self.cfg.path_l1b_enmap_image_gapfill

        # input validation
        if not os.path.isdir(path_enmap_image) and \
           not (os.path.exists(path_enmap_image) and path_enmap_image.endswith('.zip')):
            raise ValueError("The parameter 'path_enmap_image' must be a directory or the path to an existing zip "
                             "archive. Received %s." % path_enmap_image)

        # extract L1B image archive if needed
        if path_enmap_image.endswith('.zip'):
            path_enmap_image = self.extract_zip_archive(path_enmap_image, subdir='image_main')

        # extract L1B gap fill image archive if needed
        if path_enmap_image_gapfill and path_enmap_image_gapfill.endswith('.zip'):
            path_enmap_image_gapfill = self.extract_zip_archive(path_enmap_image_gapfill, subdir='image_gapfill')

        # run the reader
        RD = L1B_Reader(config=self.cfg)
        self.L1_obj = RD.read_inputdata(root_dir_main=path_enmap_image, root_dir_ext=path_enmap_image_gapfill,
                                        n_line_ext=self.cfg.n_lines_to_append)

    def run_toaRad2toaRef(self):
        """Run conversion from TOA radiance to TOA reflectance."""
        # get a new instance of radiometric transformer
        from ..processors.radiometric_transform import Radiometric_Transformer
        RT = Radiometric_Transformer(self.cfg)

        # run transformation to TOARef
        self.L1_obj = RT.transform_TOARad2TOARef(self.L1_obj)

    def run_dem_processor(self):
        self.L1_obj.get_preprocessed_dem()

    def run_geometry_processor(self):
        pass

    def run_atmospheric_correction(self):
        """Run atmospheric correction only."""
        self.L1_obj.run_AC()

    def run_orthorectification(self):
        """Run orthorectification to transform L1 sensor geometry image to L2 map geometry."""
        # run transformation from sensor to map geometry
        self.L2_obj = EnMAPL2Product_MapGeo.from_L1B_sensorgeo(config=self.cfg, enmap_ImageL1=self.L1_obj)

    def write_output(self):
        if self.cfg.output_dir:
            try:
                self.L2_obj.save(self.cfg.output_dir)
            except NotImplementedError:
                self.L2_obj.logger.warning('Currently L2A EnMAP images cannot be written to disk. '
                                           'Writing level 1 image instead.')
                self.L1_obj.save(self.cfg.output_dir)

    def run_all_processors(self):
        """Run all processors at once."""
        try:
            self.read_L1B_data()
            if self.cfg.run_deadpix_P:
                self.L1_obj.correct_dead_pixels()
            # self.run_toaRad2toaRef()  # this is only needed for geometry processor but AC expects radiance
            self.run_dem_processor()
            if self.cfg.enable_ac:
                self.L1_obj.logger.info('Skipping atmospheric correction as configured and '
                                        'computing top-of-atmosphere reflectance instead.')
                self.run_atmospheric_correction()
            else:
                self.run_toaRad2toaRef()
            self.run_geometry_processor()
            self.run_orthorectification()
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
