# -*- coding: utf-8 -*-

from datetime import datetime
import logging

from ..model.images import EnMAPL1Product_ImGeo
from ..model.metadata import EnMAP_Metadata_ImGeo


class L1B_Reader(object):
    """Reader for EnMAP Level-1B products."""

    def __init__(self, logger=None, **user_inputs):
        self.logger = logger or logging.getLogger(__name__)
        self.cfg = user_inputs

    @staticmethod
    def read_inputdata(root_dir, observation_time: datetime, lon_lat_smpl=(15, 15), nsmile_coef=4):
        # TODO move to init?
        """Read L1B, DEM and spatial reference data

        :param root_dir: Root directory of EnMAP Level-1B product
        :param observation_time: datetime of observation time (currently missing in metadata)
        :param lon_lat_smpl: number if sampling points in lon, lat fields
        :param nsmile_coef: number of polynomial coefficients for smile
        :return: instance of EnMAPL1Product_MapGeo
        """
        # get a new instance of EnMAPL1Product_MapGeo
        L1_obj = EnMAPL1Product_ImGeo(root_dir)

        # read metadata
        L1_obj.meta = EnMAP_Metadata_ImGeo(L1_obj.paths.metaxml)
        L1_obj.meta.read_metadata(observation_time=observation_time, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)

        # read VNIR data
        # call L1_obj.vnir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
        L1_obj.vnir.data = L1_obj.paths.vnir.imagedata
        L1_obj.vnir.mask_clouds = L1_obj.paths.vnir.mask_clouds
        L1_obj.vnir.deadpixelmap = L1_obj.paths.vnir.deadpixelmap
        L1_obj.vnir.meta = L1_obj.meta.vnir

        # read SWIR data
        # call L1_obj.swir.arr.setter which sets L1_obj.swir.arr to an instance of GeoArray class
        L1_obj.swir.data = L1_obj.paths.swir.imagedata
        L1_obj.swir.mask_clouds = L1_obj.paths.swir.mask_clouds
        L1_obj.swir.deadpixelmap = L1_obj.paths.swir.deadpixelmap
        L1_obj.swir.meta = L1_obj.meta.swir

        # compute radiance and calculate snr
        L1_obj.DN2TOARadiance()
        L1_obj.vnir.meta.calc_snr(data=L1_obj.vnir.data)
        L1_obj.swir.meta.calc_snr(data=L1_obj.swir.data)

        return L1_obj

    def validate_input(self):
        pass

    def validate_output(self):
        pass
