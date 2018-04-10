# -*- coding: utf-8 -*-
"""Reader module for reading all kinds of EnMAP images."""

from datetime import datetime
import logging
import tempfile
import zipfile
import numpy as np
from scipy.interpolate import interp1d

from ..model.images import EnMAPL1Product_SensorGeo
from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo
from ..options.config import EnPTConfig


class L1B_Reader(object):
    """Reader for EnMAP Level-1B products."""

    def __init__(self, config: EnPTConfig, logger: logging.Logger=None):
        """Get an instance of L1B_Reader."""
        self.cfg = config
        self.logger = logger or logging.getLogger(__name__)

    def read_inputdata(self, root_dir, observation_time: datetime, lon_lat_smpl: tuple=(15, 15), nsmile_coef: int=5,
                       compute_snr: bool=True):
        # TODO move to init?
        """Read L1B, DEM and spatial reference data.

        :param root_dir: Root directory of EnMAP Level-1B product
        :param observation_time: datetime of observation time (currently missing in metadata)
        :param lon_lat_smpl: number if sampling points in lon, lat fields
        :param nsmile_coef: number of polynomial coefficients for smile
        :return: instance of EnMAPL1Product_MapGeo
        """
        # get a new instance of EnMAPL1Product_MapGeo
        L1_obj = EnMAPL1Product_SensorGeo(root_dir, config=self.cfg)

        # read metadata
        L1_obj.meta = EnMAP_Metadata_L1B_SensorGeo(L1_obj.paths.metaxml, config=self.cfg)
        L1_obj.meta.read_metadata(observation_time=observation_time, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)

        # read VNIR data
        # call L1_obj.vnir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
        L1_obj.vnir.data = L1_obj.paths.vnir.data
        L1_obj.vnir.mask_clouds = L1_obj.paths.vnir.mask_clouds
        L1_obj.vnir.deadpixelmap = L1_obj.paths.vnir.deadpixelmap
        L1_obj.vnir.detector_meta = L1_obj.meta.vnir

        # read SWIR data
        # call L1_obj.swir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
        L1_obj.swir.data = L1_obj.paths.swir.data
        L1_obj.swir.mask_clouds = L1_obj.paths.swir.mask_clouds
        L1_obj.swir.deadpixelmap = L1_obj.paths.swir.deadpixelmap
        L1_obj.swir.detector_meta = L1_obj.meta.swir

        # compute radiance
        L1_obj.DN2TOARadiance()

        # compute SNR
        if compute_snr:
            with tempfile.TemporaryDirectory() as tmpDir, zipfile.ZipFile(self.cfg.path_l1b_snr_model, "r") as zf:
                zf.extractall(tmpDir)

                if L1_obj.meta.vnir.unitcode == 'TOARad':
                    L1_obj.vnir.detector_meta.calc_snr_from_radiance(rad_data=L1_obj.vnir.data, dir_snr_models=tmpDir)
                if L1_obj.meta.swir.unitcode == 'TOARad':
                    L1_obj.swir.detector_meta.calc_snr_from_radiance(rad_data=L1_obj.swir.data, dir_snr_models=tmpDir)

        return L1_obj

    def validate_input(self):
        """Validate user inputs."""
        pass

    def validate_output(self):
        """Validate outputs of L1B_Reader."""
        pass


def Solar_Irradiance_reader(path_solar_irr_model: str, resol_nm: float=None, wvl_min_nm: float=None,
                            wvl_max_nm: float=None) -> np.ndarray:
    """Read the given solar irradiance file and return an array of irradiances.

    :param path_solar_irr_model:    file path to solar irradiance model
                                    -> must be arranged like that:
                                        col0 = Wavelength[nm]; col1 = Solar Irradiance [W/m2/µm])
    :param resol_nm:                spectral resolution for returned irradiances [nanometers]
    :param wvl_min_nm:              minimum wavelength of returned irradiances [nanometers]
    :param wvl_max_nm:              maximum wavelength of returned irradiances [nanometers]
    :return:
    """
    sol_irr = np.loadtxt(path_solar_irr_model, skiprows=1)
    if resol_nm is not None and isinstance(resol_nm, (float, int)):
        wvl_min = (np.min(sol_irr[:, 0]) if wvl_min_nm is None else wvl_min_nm)
        wvl_max = (np.max(sol_irr[:, 0]) if wvl_max_nm is None else wvl_max_nm)
        wvl_rsp = np.arange(wvl_min, wvl_max, resol_nm)
        sol_irr = interp1d(sol_irr[:, 0], sol_irr[:, 1], kind='linear')(wvl_rsp)
    return sol_irr
