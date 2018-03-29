# -*- coding: utf-8 -*-
"""Reader module for reading all kinds of EnMAP images."""

from datetime import datetime
import logging

from ..model.images import EnMAPL1Product_SensorGeo
from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo


class L1B_Reader(object):
    """Reader for EnMAP Level-1B products."""

    def __init__(self, logger=None, **user_inputs):
        """Get an instance of L1B_Reader."""
        self.logger = logger or logging.getLogger(__name__)
        self.cfg = user_inputs

    def read_inputdata(self, root_dir, observation_time: datetime, lon_lat_smpl: tuple=(15, 15), nsmile_coef: int=5,
                       snr_vnir: str=None, snr_swir: str=None):
        # TODO move to init?
        """Read L1B, DEM and spatial reference data.

        :param root_dir: Root directory of EnMAP Level-1B product
        :param observation_time: datetime of observation time (currently missing in metadata)
        :param lon_lat_smpl: number if sampling points in lon, lat fields
        :param nsmile_coef: number of polynomial coefficients for smile
        :return: instance of EnMAPL1Product_MapGeo
        """
        # get a new instance of EnMAPL1Product_MapGeo
        L1_obj = EnMAPL1Product_SensorGeo(root_dir)

        # read metadata
        L1_obj.meta = EnMAP_Metadata_L1B_SensorGeo(L1_obj.paths.metaxml)
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

        # compute radiance and calculate snr
        L1_obj.DN2TOARadiance()
        if snr_vnir is not None and L1_obj.meta.vnir.unit == 'mW m^-2 sr^-1 nm^-1':
            self.logger.info("Compute SNR for vnir: %s" % snr_vnir)
            L1_obj.vnir.detector_meta.calc_snr_vnir(detector=L1_obj.vnir, snr_data_fn=snr_vnir)
        if snr_swir is not None and L1_obj.meta.vnir.unit == 'mW m^-2 sr^-1 nm^-1':
            self.logger.info("Compute SNR for swir: %s" % snr_swir)
            L1_obj.swir.detector_meta.calc_snr_swir(detector=L1_obj.swir, snr_data_fn=snr_swir)

        return L1_obj

    def validate_input(self):
        """Validate user inputs."""
        pass

    def validate_output(self):
        """Validate outputs of L1B_Reader."""
        pass
