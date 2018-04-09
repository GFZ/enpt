# -*- coding: utf-8 -*-
"""EnPT metadata modules. All object and functions regarding EnMAP metadata are implemented here."""

from datetime import datetime
from xml.etree import ElementTree
import logging
import os
from typing import Union
import numpy as np
from scipy.interpolate import interp2d
import spectral as sp
from geoarray import GeoArray

from ..options.config import EnPTConfig
from .srf import SRF


L1B_product_props = dict(
    xml_detector_label=dict(
        VNIR='detector1',
        SWIR='detector2'),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2')
)


#########################################################
# EnPT metadata objects for EnMAP data in sensor geometry
#########################################################


class EnMAP_Metadata_L1B_Detector_SensorGeo(object):
    """Class for all EnMAP metadata associated with a single EnMAP detector in sensor geometry.

    NOTE:
        - All metadata that have VNIR and SWIR detector in sensor geometry in common should be included here.

    """

    def __init__(self, detector_name: str, config: EnPTConfig, logger: logging.Logger=None):
        """Get a metadata object containing the metadata of a single EnMAP detector in sensor geometry.

        :param detector_name: Name of the detector (VNIR or SWIR)
        :param logger:
        """
        self.cfg = config
        self.detector_name = detector_name  # type: str
        self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.logger = logger or logging.getLogger()

        self.fwhm = None  # type: np.ndarray  # Full width half maximmum for each EnMAP band
        self.wvl_center = None  # type: np.ndarray  # Center wavelengths for each EnMAP band
        self.srf = None  # type: SRF  # SRF object holding the spectral response functions for each EnMAP band
        self.solar_irrad = None  # type: dict  # solar irradiance in [W/m2/nm] for eac band
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
        self.lat_UL_UR_LL_LR = None  # type:  list  # latitude coordinates for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR = None  # type:  list  # longitude coordinates for UL, UR, LL, LR
        self.lats = None  # type: np.ndarray  # 2D array of latitude coordinates according to given lon/lat sampling
        self.lons = None  # type: np.ndarray  # 2D array of longitude coordinates according to given lon/lat sampling
        self.unit = ''  # type: str  # radiometric unit of pixel values
        self.unitcode = ''  # type: str  # code of radiometric unit
        self.snr = None  # type: np.ndarray  # Signal to noise ratio as computed from radiance data

    def read_metadata(self, path_xml, lon_lat_smpl, nsmile_coef):
        """Read the metadata of a specific EnMAP detector in sensor geometry.

        :param path_xml:  file path of the metadata XML file
        :param lon_lat_smpl:  number if sampling points in lon, lat fields
        :param nsmile_coef:  number of polynomial coefficients for smile
        """
        xml = ElementTree.parse(path_xml).getroot()
        lbl = self.detector_label
        self.logger.info("Load data for: %s" % lbl)

        self.fwhm = np.array(xml.findall("%s/fwhm" % lbl)[0].text.replace("\n", "").split(), dtype=np.float)
        self.wvl_center = np.array(
            xml.findall("%s/centre_wavelength" % lbl)[0].text.replace("\n", "").split(), dtype=np.float)
        self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
        self.solar_irrad = self.calc_solar_irradiance_CWL_FWHM_per_band()
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
        self.lat_UL_UR_LL_LR = \
            [float(xml.findall("%s/geometry/bounding_box/%s_northing" % (lbl, corner))[0].text.split()[0])
             for corner in ("UL", "UR", "LL", "LR")]
        self.lon_UL_UR_LL_LR = \
            [float(xml.findall("%s/geometry/bounding_box/%s_easting" % (lbl, corner))[0].text.split()[0])
             for corner in ("UL", "UR", "LL", "LR")]
        self.lats = self.interpolate_corners(*self.lat_UL_UR_LL_LR, *lon_lat_smpl)
        self.lons = self.interpolate_corners(*self.lon_UL_UR_LL_LR, *lon_lat_smpl)

        try:
            self.unitcode = xml.findall("%s/unitcode" % lbl)[0].text
            # '" ".join(xml.findall("%s/radiance_unit" % lbl)[0].text.split())
            self.unit = xml.findall("%s/unit" % lbl)[0].text
        except IndexError:
            self.unitcode = 'DN'
            self.unit = 'none'
        except Exception:
            raise

        self.snr = None

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

    def calc_snr_from_radiance(self, rad_data: Union[GeoArray, np.ndarray], dir_snr_models: str):
        """Compute EnMAP SNR from radiance data for the given detector.

        SNR equation:    SNR = p0 + p1 * LTOA + p2 * LTOA ^ 2   [W/(m^2 sr nm)].
        NOTE:   The SNR model files (SNR_D1.bsq/SNR_D2.bsq) contain polynomial coefficients needed to compute SNR.
                SNR_D1.bsq: SNR model for EnMAP SWIR detector (contains high gain and low gain model coefficients)
                            - 1000 columns (for 1000 EnMAP columns)
                            - 88 bands for 88 EnMAP VNIR bands
                            - 7 lines:  - 3 lines: high gain coefficients
                                        - 3 lines: low gain coefficients
                                        - 1 line: threshold needed to decide about high gain or low gain
                SNR_D2.bsq: SNR model for EnMAP SWIR detector
                            - 1000 columns (for 1000 EnMAP columns)
                            - x bands for x EnMAP SWIR bands
                            - 3 lines for 3 coefficients

        :param rad_data:        image radiance data of EnMAP_Detector_SensorGeo
        :param dir_snr_models:  root directory where SNR model data is stored (must contain SNR_D1.bsq/SNR_D2.bsq)
        """
        path_snr_model = os.path.join(dir_snr_models, "SNR_D1.hdr" if self.detector_name == 'VNIR' else "SNR_D2.hdr")
        rad_data = np.array(rad_data)

        assert self.unitcode == 'TOARad'
        self.logger.info("Computing SNR for: %s using: %s" % (self.detector_name, path_snr_model))

        if self.detector_name == 'VNIR':
            gA = sp.open_image(path_snr_model)
            coeffs_highgain = gA[0:3, :, :]
            coeffs_lowgain = gA[3:6, :, :]
            gain_threshold = np.squeeze(gA[6, :, :])

            self.snr = np.zeros(rad_data.shape)
            for irow in range(self.nrows):
                highgain_mask = rad_data[irow, :, :] < gain_threshold  # a single row
                rad_highgain = rad_data[irow, highgain_mask]
                self.snr[irow, :, :][highgain_mask] = \
                    coeffs_highgain[0, highgain_mask] + \
                    coeffs_highgain[1, highgain_mask] * rad_highgain + \
                    coeffs_highgain[2, highgain_mask] * rad_highgain ** 2

                lowgain_mask = rad_data[irow, :, :] >= gain_threshold
                rad_lowgain = rad_data[irow, lowgain_mask]
                self.snr[irow, :, :][lowgain_mask] = \
                    coeffs_lowgain[0, lowgain_mask] + \
                    coeffs_lowgain[1, lowgain_mask] * rad_lowgain + \
                    coeffs_lowgain[2, lowgain_mask] * rad_lowgain ** 2

        else:
            coeffs = sp.open_image(path_snr_model)[:, :, :]
            self.snr = coeffs[0, :, :] + coeffs[1, :, :] * rad_data[:, :, :] + coeffs[2, :, :] * rad_data[:, :, :] ** 2

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

    def calc_solar_irradiance_CWL_FWHM_per_band(self) -> dict:
        from ..io.reader import Solar_Irradiance_reader

        self.logger.debug('Calculating solar irradiance...')

        sol_irr = Solar_Irradiance_reader(path_solar_irr_model=self.cfg.path_solar_irr, wvl_min_nm=350, wvl_max_nm=2500)

        irr_bands = {}
        for band in self.srf.bands:
            if self.srf[band] is None:
                irr_bands[band] = None
            else:
                WVL_band = self.srf.srfs_wvl if self.srf.wvl_unit == 'nanometers' else self.srf.srfs_wvl * 1000
                RSP_band = self.srf.srfs_norm01[band]
                sol_irr_at_WVL = np.interp(WVL_band, sol_irr[:, 0], sol_irr[:, 1], left=0, right=0)

                irr_bands[band] = round(np.sum(sol_irr_at_WVL * RSP_band) / np.sum(RSP_band), 2)

        return irr_bands


class EnMAP_Metadata_L1B_SensorGeo(object):
    """EnMAP Metadata class holding the metadata of the complete EnMAP product in sensor geometry incl. VNIR and SWIR.

    Attributes:
        - logger(logging.Logger):  None or logging instance
        - observation_datetime(datetime.datetime):  datetime of observation time (currently missing in metadata)
        - vnir(EnMAP_Metadata_VNIR_SensorGeo)
        - swir(EnMAP_Metadata_SWIR_SensorGeo)

    """

    def __init__(self, path_metaxml, config: EnPTConfig, logger=None):
        """Get a metadata object instance for the given EnMAP L1B product in sensor geometry.

        :param path_metaxml:  file path of the EnMAP L1B metadata XML file
        :param logger:  instance of logging.logger or subclassed
        """
        self.cfg = config
        self.logger = logger or logging.getLogger()
        self._path_xml = path_metaxml

        # defaults
        self.observation_datetime = None  # type: datetime  # Date and Time of image observation
        self.earthSunDist = None  # type: float  # earth-sun distance # TODO doc correct?
        self.vnir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of VNIR only
        self.swir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of SWIR only
        self.detector_attrNames = ['vnir', 'swir']

    def read_common_meta(self, observation_time: datetime=None):
        """Read the metadata belonging to both, the VNIR and SWIR detector of the EnMAP L1B product in sensor geometry.

        :param observation_time:  date and time of image observation (datetime.datetime)
        """
        # FIXME observation time is currently missing in the XML
        self.observation_datetime = observation_time
        self.earthSunDist = self.get_earth_sun_distance(self.observation_datetime)

    def get_earth_sun_distance(self, acqDate: datetime):
        """Get earth sun distance (requires file of pre calculated earth sun distance per day)

        :param acqDate:
        """

        if not os.path.exists(self.cfg.path_earthSunDist):
            self.logger.warning("\n\t WARNING: Earth Sun Distance is assumed to be "
                                "1.0 because no database can be found at %s.""" % self.cfg.path_earthSunDist)
            return 1.0
        if not acqDate:
            self.logger.warning("\n\t WARNING: Earth Sun Distance is assumed to be 1.0 because "
                                "acquisition date could not be read from metadata.")
            return 1.0

        with open(self.cfg.path_earthSunDist, "r") as EA_dist_f:
            EA_dist_dict = {}
            for line in EA_dist_f:
                date, EA = [item.strip() for item in line.split(",")]
                EA_dist_dict[date] = EA

        return float(EA_dist_dict[acqDate.strftime('%Y-%m-%d')])

    def read_metadata(self, observation_time: datetime, lon_lat_smpl, nsmile_coef):
        """Read the metadata of the whole EnMAP L1B product in sensor geometry.

        :param observation_time:  date and time of image observation (datetime.datetime)
        :param lon_lat_smpl:  number if sampling points in lon, lat fields
        :param nsmile_coef:  number of polynomial coefficients for smile
        """
        self.read_common_meta(observation_time)
        self.vnir = EnMAP_Metadata_L1B_Detector_SensorGeo('VNIR', config=self.cfg, logger=self.logger)
        self.vnir.read_metadata(self._path_xml, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)
        self.swir = EnMAP_Metadata_L1B_Detector_SensorGeo('SWIR', config=self.cfg, logger=self.logger)
        self.swir.read_metadata(self._path_xml, lon_lat_smpl=lon_lat_smpl, nsmile_coef=nsmile_coef)
