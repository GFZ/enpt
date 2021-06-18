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

"""EnPT process controller module."""

import os
import sys
import tempfile
import zipfile
import shutil
import weakref
import warnings
import pickle
from typing import Optional
from time import time
from datetime import timedelta
from glob import glob

from ..options.config import EnPTConfig
from ..io.reader import L1B_Reader
from ..model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo  # noqa: F401

__author__ = 'Daniel Scheffler'


class EnPT_Controller(object):
    """Class of EnPT process controller."""

    def __init__(self, config: EnPTConfig = None, **config_kwargs: dict):
        """Initialize the Process Controller.

        :param config:          an instance of the EnPTConfig class (overrides config_kwargs)
        :param config_kwargs:   configuration parameters to be passed to EnPTConfig class
        """
        self.cfg: EnPTConfig = config or EnPTConfig(**config_kwargs)

        # generate temporary directory (must be deleted at the end of the pipeline)
        self.tempDir = tempfile.mkdtemp()

        # setup a finalizer that destroys remaining data (directories, etc.) in case of unexpected exit
        self._finalizer = weakref.finalize(self, self._cleanup, self.tempDir,
                                           warn_message="Implicitly cleaning up {!r}".format(self))

        # record startup time
        self._time_startup = time()

        # defaults
        self.L1_obj: Optional[EnMAPL1Product_SensorGeo] = None
        self.L2_obj: Optional[EnMAPL2Product_MapGeo] = None

    def extract_zip_archive(self, path_zipfile: str, subdir: str = '') -> str:
        """Extract the given EnMAP image zip archive and return the L1B image root directory path.

        :param path_zipfile:    /path/to/zipfile.zip
        :param subdir:          subdirectory name to be created within temporary directory
        :return:                /tmp/tmpk2qp0yri/rootdir/
        """
        outdir = os.path.join(self.tempDir, subdir) if not self.cfg.is_dummy_dataformat else self.tempDir

        with zipfile.ZipFile(path_zipfile, "r") as zf:
            zf.extractall(outdir)

        # move the data one level up in case they are within a sub-folder in the zip file
        if not self.cfg.is_dummy_dataformat:
            content = glob(os.path.join(outdir, '*'))

            if len(content) == 1 and os.path.isdir(content[0]):
                for fp in glob(os.path.join(outdir, '**', '*')):
                    shutil.move(fp, outdir)
                shutil.rmtree(content[0])

        if not os.path.isdir(outdir):
            raise NotADirectoryError(outdir)

        if not self.cfg.is_dummy_dataformat:
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

    def run_spatial_optimization(self):
        # get a new instance of radiometric transformer
        from ..processors.spatial_optimization import Spatial_Optimizer
        SpO = Spatial_Optimizer(self.cfg)

        # run optimization
        self.L1_obj = SpO.optimize_geolayer(self.L1_obj)

    def run_atmospheric_correction(self):
        """Run atmospheric correction only."""
        self.L1_obj.run_AC()

    def run_orthorectification(self):
        """Run orthorectification to transform L1 sensor geometry image to L2 map geometry."""
        # run transformation from sensor to map geometry
        self.L2_obj = EnMAPL2Product_MapGeo.from_L1B_sensorgeo(config=self.cfg, enmap_ImageL1=self.L1_obj)

    def write_output(self):
        if self.cfg.output_dir:
            if self.L2_obj is not None:
                self.L2_obj.save(self.cfg.output_dir)
            else:
                self.L1_obj.save(self.cfg.output_dir)

    def run_all_processors(self):
        """Run all processors at once."""
        if os.getenv('IS_ENPT_GUI_TEST') != "1":
            try:
                if os.getenv('IS_ENPT_GUI_CALL') == "1":
                    self._write_to_stdout_stderr()

                if self.cfg.log_level == 'DEBUG':
                    self._print_received_configuration()

                self.read_L1B_data()
                if self.cfg.run_deadpix_P:
                    self.L1_obj.correct_dead_pixels()
                if self.cfg.enable_absolute_coreg:
                    # self.run_toaRad2toaRef()  # this is only needed for geometry processor but AC expects radiance
                    self.run_spatial_optimization()
                self.run_dem_processor()
                if self.cfg.enable_ac:
                    self.run_atmospheric_correction()

                    # re-apply dead pixel correction
                    self.L1_obj.logger.info(
                        'Re-applying dead pixel correction to correct for spectral spikes due to fringe effect.')
                    self.L1_obj.correct_dead_pixels()
                else:
                    self.L1_obj.logger.info('Skipping atmospheric correction as configured and '
                                            'computing top-of-atmosphere reflectance instead.')
                    self.run_toaRad2toaRef()
                # self.run_spatial_optimization()
                self.run_orthorectification()
                self.write_output()

                self.L1_obj.logger.info('Total runtime of the processing chain: %s'
                                        % timedelta(seconds=time() - self._time_startup))
            finally:
                self.cleanup()

        else:
            self._print_received_configuration()

            if not os.path.isdir(self.cfg.output_dir):
                raise NotADirectoryError(self.cfg.output_dir)

            with open(os.path.join(self.cfg.output_dir, 'received_args_kwargs.pkl'), 'wb') as outF:
                pickle.dump(
                    dict(
                        json_config=self.cfg.json_config,
                        kwargs=self.cfg.kwargs
                    ),
                    outF)

    @staticmethod
    def _write_to_stdout_stderr():
        """Write to STDOUT and STDERR to reveal Queue.name of these streams to enpt_enmapboxapp."""
        sys.stdout.write('Connecting to EnPT STDOUT stream.\n')
        sys.stdout.flush()
        sys.stderr.write('Connecting to EnPT STDERR stream.\n')
        sys.stderr.flush()

    def _print_received_configuration(self):
        print('EnPT Controller received the following configuration:')
        print(repr(self.cfg))

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
