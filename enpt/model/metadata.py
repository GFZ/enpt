# -*- coding: utf-8 -*-

from datetime import datetime
from xml.etree import ElementTree
import logging
import numpy as np
from scipy.interpolate import interp2d

L1B_product_props = dict(
    xml_detector_label=dict(
        VNIR='detector1',
        SWIR='detector2'),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2')
)


##############
# BASE CLASSES
##############


class _EnMAP_Metadata_Detector_ImGeo(object):
    """Base class for all EnMAP metadata associated with a single EnMAP detector in image geometry.

    NOTE:
        - All metadata that have VNIR and SWIR detector in image geometry in common should be included here.
        - All classes representing VNIR or SWIR specific metadata should inherit from _EnMAP_Metadata_Detector_ImGeo.

    """
    def __init__(self, detector_name: str, logger: logging.Logger=None):
        """

        :param detector_name: Name of the detector (VNIR or SWIR)
        :param logger:
        """

        self.detector_name = detector_name  # type: str
        self.logger = logger or logging.getLogger()

        self.fwhm = None  # type: np.ndarray  # Full width half maximmum for each EnMAP band
        self.wvl_center = None  # type: np.ndarray  # Center wavelengths for each EnMAP band
        self.nwvl = None  # type: int  # Number of wave bands
        self.nrows = None  # type: int  # number of rows
        self.ncols = None  # type: int  # number of columns
        self.smile_coef = None  # type: np.ndarray  # smile coefficients needed for smile computation
        self.nsmile_coef = None  # type: int  # number of smile coefficients
        self.smile = None  # type: np.ndarray  # smile for each EnMAP image column
        self.l_min = None  # type: np.ndarray
        self.l_max = None  # type: np.ndarray
        self.geom_view_zenith = None  # type: float  # viewing zenith angle
        self.geom_view_azimuth = None  # type: float  # viewing azimuth angle
        self.geom_sun_zenith = None  # type: float  # sun zenith angle
        self.geom_sun_azimuth = None  # type: float  # sun azimuth angle
        self.mu_sun = None  # type: float  # earth-sun distance # TODO doc correct?
        self.lat_UL_UR_LL_LR = None  # type:  list  # latitude coordinates for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR = None  # type:  list  # longitude coordinates for UL, UR, LL, LR
        self.lats = None  # type: np.ndarray  # 2D array of latitude coordinates according to given lon/lat sampling
        self.lons = None  # type: np.ndarray  # 2D array of longitude coordinates according to given lon/lat sampling
        self.unit = ''  # type: str  # radiometric unit of pixel values
        self.snr = None  # type: np.ndarray  # Signal to noise ratio as computed from radiance data

    def _read_metadata(self, path_xml, detector_label_xml, lon_lat_smpl=(15, 15), nsmile_coef=4):
        """Read the metadata of a specific EnMAP detector in image geometry.

        :param path_xml:  file path of the metadata XML file
        :param detector_label_xml:
        :param lon_lat_smpl:  number if sampling points in lon, lat fields
        :param nsmile_coef:  number of polynomial coefficients for smile
        """
        xml = ElementTree.parse(path_xml).getroot()
        lbl = detector_label_xml
        self.logger.info("Load data for: %s" % lbl)

        self.fwhm = np.array(xml.findall("%s/fwhm" % lbl)[0].text.replace("\n", "").split(), dtype=np.float)
        self.wvl_center = np.array(
            xml.findall("%s/centre_wavelength" % lbl)[0].text.replace("\n", "").split(), dtype=np.float)
        self.nwvl = len(self.wvl_center)
        self.nrows = np.int(xml.findall("%s/rows" % lbl)[0].text)
        self.ncols = np.int(xml.findall("%s/columns" % lbl)[0].text)
        self.smile_coef = np.array(xml.findall("%s/smile" % lbl)[0].text.replace("\n", "").split(), dtype=np.float) \
                            .reshape((-1, nsmile_coef + 1))[:, 1:]
        self.nsmile_coef = nsmile_coef
        self.smile = self.calc_smile()
        self.l_min = np.array(xml.findall("%s/L_min" % lbl)[0].text.split(), dtype=np.float)
        self.l_max = np.array(xml.findall("%s/L_max" % lbl)[0].text.split(), dtype=np.float)
        self.geom_view_zenith = np.float(
            xml.findall("%s/observation_geometry/zenith_angle" % lbl)[0].text.split()[0])
        self.geom_view_azimuth = np.float(
            xml.findall("%s/observation_geometry/azimuth_angle" % lbl)[0].text.split()[0])
        self.geom_sun_zenith = np.float(
            xml.findall("%s/illumination_geometry/zenith_angle" % lbl)[0].text.split()[0])
        self.geom_sun_azimuth = np.float(
            xml.findall("%s/illumination_geometry/azimuth_angle" % lbl)[0].text.split()[0])
        self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))
        self.lat_UL_UR_LL_LR = \
            [float(xml.findall("%s/geometry/bounding_box/%s_northing" % (lbl, corner))[0].text.split()[0])
             for corner in ("UL", "UR", "LL", "LR")]
        self.lon_UL_UR_LL_LR = \
            [float(xml.findall("%s/geometry/bounding_box/%s_easting" % (lbl, corner))[0].text.split()[0])
             for corner in ("UL", "UR", "LL", "LR")]
        self.lats = self.interpolate_corners(*self.lat_UL_UR_LL_LR, *lon_lat_smpl)
        self.lons = self.interpolate_corners(*self.lon_UL_UR_LL_LR, *lon_lat_smpl)
        self.unit = 'none'  # '" ".join(xml.findall("%s/radiance_unit" % lbl)[0].text.split())
        self.unitcode = 'DN'

    def calc_smile(self):
        """Compute smile for each EnMAP column.
        The sum in (1) is expressed as inner product of two arrays with inner dimension as the polynomial smile
        coefficients shape = (ncols, nsmile_coef) of polynomial x

        :return:
        """
        # smile(icol, iwvl) = sum_(p=0)^(nsmile_coef-1) smile_coef[iwvl, p] * icol**p (1)
        return np.inner(
            np.array([[icol ** p for p in np.arange(self.nsmile_coef)] for icol in np.arange(self.ncols)]),
            self.smile_coef  # shape = (nwvl, nsmile_coef)
        )  # shape = (ncols, nwvl)

    def calc_snr(self, data: np.ndarray):
        """Compute EnMAP snr from radiance data.

        :param data: Numpy array with radiance for scene
        """
        self.logger.info("Compute snr for: %s" % self.detector_name)
        self.logger.warning("SNR model missing -> const. value of 500 is returned")
        return 500 * np.ones(data.shape, dtype=np.float)

    @staticmethod
    def interpolate_corners(ul: float, ur: float, ll: float, lr: float, nx: int, ny: int):
        """Compute interpolated field from corner values of a scalar field given at: ul, ur, ll, lr.

        :param nx, ny: final shape
        """
        ff = interp2d(x=[0, 1], y=[0, 1], z=[[ul, ur], [ll, lr]], kind='linear')
        rr = np.zeros((nx, ny), dtype=np.float)
        for i, x in enumerate(np.linspace(0, 1, nx)):
            for j, y in enumerate(np.linspace(0, 1, ny)):
                rr[i, j] = ff(x, y)
        return rr


########################################################
# EnPT metadata objects for EnMAP data in image geometry
########################################################


class EnMAP_Metadata_ImGeo(object):
    """EnMAP Metadata class holding the metadata of the complete EnMAP product in image geometry
    including VNIR and SWIR detector.

    Attributes:
        - logger(logging.Logger):  None or logging instance
        - observation_datetime(datetime.datetime):  datetime of observation time (currently missing in metadata)
        - vnir(EnMAP_Metadata_VNIR_ImGeo)
        - swir(EnMAP_Metadata_SWIR_ImGeo)
    """
    def __init__(self, path_metaxml, logger=None):
        self.logger = logger or logging.getLogger()
        self._path_xml = path_metaxml

        # defaults
        self.observation_datetime = None  # type: datetime  # Date and Time of image observation
        self.vnir = None  # type: EnMAP_Metadata_VNIR_ImGeo # metadata of VNIR only
        self.swir = None  # type: EnMAP_Metadata_SWIR_ImGeo # metadata of SWIR only

    def read_common_meta(self, observation_time: datetime=None):
        # FIXME observation time is currently missing in the XML
        self.observation_datetime = observation_time

    def read_metadata(self, observation_time: datetime=None, lon_lat_smpl=(15, 15), nsmile_coef=4):
        self.read_common_meta(observation_time)
        self.vnir = EnMAP_Metadata_VNIR_ImGeo(self._path_xml)
        self.vnir.read_metadata(lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)
        self.swir = EnMAP_Metadata_SWIR_ImGeo(self._path_xml)
        self.swir.read_metadata(lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)


class EnMAP_Metadata_VNIR_ImGeo(_EnMAP_Metadata_Detector_ImGeo):
    """EnMAP Metadata class holding the metadata of the VNIR detector in image geometry.

    NOTE:
        - inherits all attributes from base class _EnMAP_Metadata_Detector_ImGeo.
    """
    def __init__(self, path_metaxml, logger=None):
        # get all attributes from base class '_EnMAP_Metadata_Detector_ImGeo'
        super(EnMAP_Metadata_VNIR_ImGeo, self).__init__('VNIR', logger=logger)
        self._path_xml = path_metaxml
        self.detector_label = L1B_product_props['xml_detector_label']['VNIR']

    def read_metadata(self, lon_lat_smpl=(15, 15), nsmile_coef=4):
        super(EnMAP_Metadata_VNIR_ImGeo, self)\
            ._read_metadata(self._path_xml, self.detector_label, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)


class EnMAP_Metadata_SWIR_ImGeo(_EnMAP_Metadata_Detector_ImGeo):
    """EnMAP Metadata class holding the metadata of the SWIR detector in image geometry.

    NOTE:
        - inherits all attributes from base class _EnMAP_Metadata_Detector_ImGeo.
    """
    def __init__(self, path_metaxml, logger=None):
        # get all attributes from base class '_EnMAP_Metadata_Detector_ImGeo'
        super(EnMAP_Metadata_SWIR_ImGeo, self).__init__('SWIR', logger=logger)
        self._path_xml = path_metaxml
        self.detector_label = L1B_product_props['xml_detector_label']['SWIR']

    def read_metadata(self, lon_lat_smpl=(15, 15), nsmile_coef=4):
        """Read the metadata of the SWIR detector in image geometry.

        :param lon_lat_smpl: number if sampling points in lon, lat fields
        :param nsmile_coef: number of polynomial coefficients for smile
        :return:
        """
        super(EnMAP_Metadata_SWIR_ImGeo, self)\
            ._read_metadata(self._path_xml, self.detector_label, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)
