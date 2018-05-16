# -*- coding: utf-8 -*-
"""Reader module for reading all kinds of EnMAP images."""

# from datetime import datetime
import logging
import tempfile
import zipfile
import numpy as np
from scipy.interpolate import interp1d

from ..model.images import EnMAPL1Product_SensorGeo
# from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo
from ..options.config import EnPTConfig


class L1B_Reader(object):
    """Reader for EnMAP Level-1B products."""

    def __init__(self,
                 config: EnPTConfig,
                 logger: logging.Logger=None,
                 root_dir_main: str=None,
                 root_dir_ext: str=None,
                 n_line_ext: int=None):
        # Add option to init as suggested.
        """Get an instance of L1B_Reader.

        :param config:  instance of EnPTConfig class
        :param logger:  instance of logging.Logger (NOTE: This logger is only used to log messages within L1B_Reader.
                                                          It is not appended to the read L1B EnMAP object).
        :param root_dir_main: Root directory of EnMAP Level-1B product (the main image)
        :param root_dir_ext:  Root directory of EnMAP Level-1B product (to extend the main image)
        :param n_line_ext:    [Optional] add number of line to be added to the main image from the extended image
        """
        self.cfg = config
        self.logger = logger or logging.getLogger(__name__)

        # read data if root_dir_main is given or not
        if root_dir_main is not None:
            self.read_inputdata(root_dir_main, root_dir_ext, n_line_ext)

    def read_inputdata(self,
                       root_dir_main,
                       root_dir_ext: str=None,
                       n_line_ext: int=None,
                       lon_lat_smpl: tuple=(15,15),
                       compute_snr: bool=True):
        # All information are read from data itself now
        # In case of multiple files, temporary files are created to store them.
        """
        Read L1B EnMAP data. Extend the image by adding a second image [entire, partial]
        :param root_dir_main: Root directory of the main EnMAP Level-1B product
        :param root_dir_ext:  Root directory of the extended EnMAP Level-1B product [optional]
        :param n_line_ext:    Number of line to be added to the main image [if None, use the whole image]
        :param lon_lat_smpl:  number if sampling points in lon, lat fields
        :param compute_snr:   whether to compute SNR or not (default: True)
        :return: instance of EnMAPL1Product_SensorGeo
        """
        self.logger.info("Reading Input Data")

        # Get a new instance of the EnMAPL1Product_SensorGeo for the main image
        l1b_main_obj = EnMAPL1Product_SensorGeo(root_dir_main, config=self.cfg, logger=self.logger,
                                                lon_lat_smpl=lon_lat_smpl)

        # load data from the main object
        l1b_main_obj.vnir.data = l1b_main_obj.paths.vnir.data
        l1b_main_obj.vnir.mask_clouds = l1b_main_obj.paths.vnir.mask_clouds
        l1b_main_obj.vnir.deadpixelmap = l1b_main_obj.paths.vnir.deadpixelmap
        l1b_main_obj.swir.data = l1b_main_obj.paths.swir.data
        l1b_main_obj.swir.mask_clouds = l1b_main_obj.paths.swir.mask_clouds
        l1b_main_obj.swir.deadpixelmap = l1b_main_obj.paths.swir.deadpixelmap
        l1b_main_obj.DN2TOARadiance()

        # in case of a second file, we create new files that will be temporary save into a temporary directory
        # and their path will be stored into a the paths of l1b_main_obj
        # NOTE: We do the following hypothesis:
        #         - The dead pixel map will not change when acquiring 2 adjacent images.
        if root_dir_ext is not None:
            l1b_ext_obj = EnMAPL1Product_SensorGeo(root_dir_ext, config=self.cfg, lon_lat_smpl=lon_lat_smpl)
            l1b_main_obj.append_new_image(l1b_ext_obj, n_line_ext)

        # compute SNR
        if compute_snr:
            with tempfile.TemporaryDirectory(dir=self.cfg.working_dir) as tmpDir,\
                    zipfile.ZipFile(self.cfg.path_l1b_snr_model, "r") as zf:
                zf.extractall(tmpDir)

                if l1b_main_obj.meta.vnir.unitcode == 'TOARad':
                    l1b_main_obj.vnir.detector_meta.calc_snr_from_radiance(rad_data=l1b_main_obj.vnir.data,
                                                                           dir_snr_models=tmpDir)
                if l1b_main_obj.meta.swir.unitcode == 'TOARad':
                    l1b_main_obj.swir.detector_meta.calc_snr_from_radiance(rad_data=l1b_main_obj.swir.data,
                                                                           dir_snr_models=tmpDir)


        # Return the l1b_main_obj
        return l1b_main_obj

    # def read_inputdata_old(self, root_dir, observation_time: datetime, lon_lat_smpl: tuple=(15, 15), nsmile_coef: int=5,
    #                    compute_snr: bool=True):
    #     # TODO move to init? --> This has been added in the init phase (will call the new read_inputdata method
    #     """Read L1B, DEM and spatial reference data.
    #
    #     :param root_dir:            Root directory of EnMAP Level-1B product
    #     :param observation_time:    datetime of observation time (currently missing in metadata)
    #     :param lon_lat_smpl:        number if sampling points in lon, lat fields
    #     :param nsmile_coef:         number of polynomial coefficients for smile
    #     :param compute_snr:         whether to compute SNR or not (default: True)
    #     :return:    instance of EnMAPL1Product_MapGeo
    #     """
    #     # get a new instance of EnMAPL1Product_MapGeo
    #     L1_obj = EnMAPL1Product_SensorGeo(root_dir, config=self.cfg)
    #
    #     # read metadata
    #     L1_obj.meta = EnMAP_Metadata_L1B_SensorGeo(L1_obj.paths.metaxml, config=self.cfg, logger=L1_obj.logger)
    #     L1_obj.meta.read_metadata(observation_time=observation_time, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)
    #
    #     # read VNIR data
    #     # call L1_obj.vnir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
    #     L1_obj.vnir.data = L1_obj.paths.vnir.data
    #     L1_obj.vnir.mask_clouds = L1_obj.paths.vnir.mask_clouds
    #     L1_obj.vnir.deadpixelmap = L1_obj.paths.vnir.deadpixelmap
    #     L1_obj.vnir.detector_meta = L1_obj.meta.vnir
    #
    #     # read SWIR data
    #     # call L1_obj.swir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
    #     L1_obj.swir.data = L1_obj.paths.swir.data
    #     L1_obj.swir.mask_clouds = L1_obj.paths.swir.mask_clouds
    #     L1_obj.swir.deadpixelmap = L1_obj.paths.swir.deadpixelmap
    #     L1_obj.swir.detector_meta = L1_obj.meta.swir
    #
    #     # compute radiance
    #     L1_obj.DN2TOARadiance()
    #
    #     # compute SNR
    #     if compute_snr:
    #         with tempfile.TemporaryDirectory(dir=self.cfg.working_dir) as tmpDir,\
    #                 zipfile.ZipFile(self.cfg.path_l1b_snr_model, "r") as zf:
    #             zf.extractall(tmpDir)
    #
    #             if L1_obj.meta.vnir.unitcode == 'TOARad':
    #                 L1_obj.vnir.detector_meta.calc_snr_from_radiance(rad_data=L1_obj.vnir.data, dir_snr_models=tmpDir)
    #             if L1_obj.meta.swir.unitcode == 'TOARad':
    #                 L1_obj.swir.detector_meta.calc_snr_from_radiance(rad_data=L1_obj.swir.data, dir_snr_models=tmpDir)
    #
    #     return L1_obj

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
