# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz.de)
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
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""Module to download the files needed for ISOFIT which cannot be included in the package."""

import urllib.request
import os
import hashlib
from logging import Logger


def download_isofit_resources(dir_output: str, logger: Logger = None):
    """
    Download surface spectra and LUT file from remote EnPT repository for running ISOFIT.

    :param dir_output:  directory where to store the downloaded files
    :param logger:      logger to use
    """
    url_enpt = "https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/"
    url_spec = url_enpt + "-/raw/main/enpt/resources/isofit/isofit_surface_spectra.zip"
    url_lut = url_enpt + "-/raw/main/enpt/resources/isofit/lut.zip"

    checksums = {
        "isofit_surface_spectra.zip": "d3b3ea65d2931240cb0adde0bb75cb75",
        "lut.zip": "1ac45c0056128d34651be0790fc16841",
    }

    for url, desc in zip([url_spec, url_lut], ["ISOFIT surface spectra", "LUT"]):
        fn = os.path.basename(url)
        p_out = os.path.join(dir_output, fn)

        download = not os.path.isfile(p_out) or md5(p_out) != checksums[fn]

        if download:
            os.makedirs(dir_output, exist_ok=True)

            md5sum = ''
            for i in range(3):
                urllib.request.urlretrieve(url, p_out)

                if os.path.isfile(p_out):
                    # verify checksum
                    md5sum = md5(p_out)
                    if md5sum == checksums[fn]:
                        break

            if not os.path.isfile(p_out):
                raise FileNotFoundError(f"Download of {desc} zipfile failed after 3 attempts. Please download it "
                                        f"manually from {url} and store it at {dir_output} directory. "
                                        "Otherwise, the ISOFIT AC will not work.")
            elif md5sum != checksums[fn]:
                raise ValueError(f"Downloaded {desc} zipfile is corrupted. "
                                 f"Please download it manually from {url} and store it at {dir_output} directory. "
                                 "Otherwise, the ISOFIT AC will not work.")

            if logger is not None:
                logger.info(f"{desc} zipfile successfully downloaded.")

        else:
            if logger is not None:
                logger.info(f"{desc} zipfile already downloaded. Proceeding.")


def md5(file_path):
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

    return file_hash.hexdigest()
