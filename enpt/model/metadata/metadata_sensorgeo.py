# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2026 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
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

"""EnPT metadata objects for EnMAP data in sensor geometry."""

from datetime import datetime
from lxml import etree as ElementTree
import logging
import os
import fnmatch
from typing import List  # noqa: F401
from collections import OrderedDict
from packaging.version import parse as parse_version
import numpy as np
from py_tools_ds.geo.vector.topology import Polygon, get_footprint_polygon  # noqa: F401  # flake8 issue
from geoarray import GeoArray

from ...options.config import EnPTConfig
from ..srf import SRF
from ...processors.spatial_transform import RPC_3D_Geolayer_Generator

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class EnMAP_Metadata_L1B_Detector_SensorGeo(object):
    """Class for all EnMAP metadata associated with a single EnMAP detector in sensor geometry.

    NOTE:
        - All metadata that have VNIR and SWIR detector in sensor geometry in common should be included here.

    """

    def __init__(self, detector_name: str, config: EnPTConfig, logger: logging.Logger = None):
        """Get a metadata object containing the metadata of a single EnMAP detector in sensor geometry.

        :param detector_name:   Name of the detector (VNIR or SWIR)
        :param config:          EnPT configuration object
        :param logger:          instance of logging.logger or subclassed
        """
        self.cfg = config
        self.detector_name: str = detector_name
        self.detector_label = detector_name.lower()
        self.logger = logger or logging.getLogger()

        # These lines are used to load path information
        self.filename_data: str = ''  # detector data filename
        self.scene_basename: str = ''  # basename of the EnMAP image
        self.filename_deadpixelmap: str = ''  # filename of the dead pixel file
        self.filename_quicklook: str = ''  # filename of the quicklook file
        self.filename_testflags: str = ''  # filename of the testflags file
        # FIXME masks of BOTH detectors
        self.filename_mask_landwater: str = ''  # filename of the land/water mask file
        self.filename_mask_snow: str = ''  # filename of the snow mask file
        self.filename_mask_cloudshadow: str = ''  # filename of the cloud shadow mask file
        self.filename_mask_clouds: str = ''  # filename of the cloud mask file
        self.filename_mask_haze: str = ''  # filename of the haze mask file
        self.filename_mask_cirrus: str = ''  # filename of the cirrus mask file

        self.wvl_center: np.ndarray | None = None  # Center wavelengths for each EnMAP band
        self.fwhm: np.ndarray | None = None  # Full width half maximum for each EnMAP band
        self.srf: SRF | None = None  # SRF object holding the spectral response functions for each EnMAP band
        self.solar_irrad: np.ndarray | None = None  # solar irradiance in [W/m2/nm] for each band
        self.nwvl: int | None = None  # Number of wave bands
        self.nrows: int | None = None  # number of rows
        self.ncols: int | None = None  # number of columns
        self.smile_coef: np.ndarray | None = None  # smile coefficients needed for smile computation
        self.nsmile_coef: int | None = None  # number of smile coefficients
        self.smile: np.ndarray | None = None  # smile for each EnMAP image column
        self.gains: np.ndarray | None = None  # band-wise gains for computing radiance from DNs
        self.offsets: np.ndarray | None = None   # band-wise offsets for computing radiance from DNs
        self.l_min: np.ndarray | None = None  # band-wise l-min for computing radiance from DNs
        self.l_max: np.ndarray | None = None  # band-wise l-max for computing radiance from DNs
        self.goodbands_inds: list | None = None  # list of band indices included in the processing (all other bands are removed)  # noqa
        self.lat_UL_UR_LL_LR: List[float, float, float, float] | None = None  # latitude coords for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR: List[float, float, float, float] | None = None  # longitude coords for UL, UR, LL, LR
        self.epsg_ortho: int | None = None  # EPSG code of the orthorectified image
        self.rpc_coeffs: OrderedDict = OrderedDict()  # RPC coefficients for geolayer computation
        self.ll_mapPoly: Polygon | None = None  # footprint polygon in longitude/latitude map coordinates
        self.lats: np.ndarray | None = None  # 2D/3D array of latitude coordinates
        self.lons: np.ndarray | None = None  # 2D/3D array of longitude coordinates
        self.geolayer_has_keystone: bool | None = None  # indicates if lon/lat geolayer considers keystone (3D array)
        self.unit: str = ''  # radiometric unit of pixel values
        self.unitcode: str = ''  # code of radiometric unit
        self.preview_wvls: list[float] | None = None  # wavelengths to be used for quicklook images
        self.preview_bands: list[int] | None = None  # band indices to be used for quicklook images
        self.snr: np.ndarray | None = None   # Signal-to-noise ratio as computed from radiance data

    def read_metadata(self, path_xml):
        """Read the metadata of a specific EnMAP detector in sensor geometry.

        :param path_xml: file path of the metadata file
        :return: None
        """
        xml = ElementTree.parse(path_xml).getroot()

        lbl = self.detector_label
        self.logger.info("Reading metadata for %s detector..." % self.detector_name)

        # read data filenames
        all_filenames = [ele.text for ele in xml.findall("product/productFileInformation/file/name")]

        def get_filename(matching_exp: str):
            matches = []
            for ext in ['', '.TIF', '.GEOTIFF', '.BSQ', '.BIL', '.BIP', 'JPEG2000', '.JP2', '.jp2']:
                matches.extend(fnmatch.filter(all_filenames, f'{matching_exp}{ext}'))

                if matches:
                    break

            if not matches:
                raise FileNotFoundError("Could not find a file that matches the expression '%s'." % matching_exp)

            return matches[0]

        self.filename_data = xml.find("product/image/%s/name" % lbl).text
        self.scene_basename = self.filename_data.split('-SPECTRAL_IMAGE')[0]
        self.filename_quicklook = xml.find("product/quicklook/%s/name" % lbl).text
        self.filename_deadpixelmap = get_filename('*QL_PIXELMASK_%s' % self.detector_name)
        self.filename_mask_landwater = get_filename('*QL_QUALITY_CLASSES')
        self.filename_mask_snow = get_filename('*QL_QUALITY_SNOW')
        self.filename_mask_cloudshadow = get_filename('*QL_QUALITY_CLOUDSHADOW')
        self.filename_mask_clouds = get_filename('*-QL_QUALITY_CLOUD')
        self.filename_mask_haze = get_filename('*QL_QUALITY_HAZE')
        self.filename_mask_cirrus = get_filename('*QL_QUALITY_CIRRUS')
        self.filename_testflags = get_filename('*QL_QUALITY_TESTFLAGS_%s' % self.detector_name)

        # FIXME combine different cloud masks?

        # read some basic information concerning the detector
        self.nrows = int(xml.find("product/image/%s/dimension/rows" % lbl).text)
        self.ncols = int(xml.find("product/image/%s/dimension/columns" % lbl).text)
        self.unitcode = 'DN'
        self.unit = ''

        # Read image coordinates
        # FIXME why do we have the same corner coordinates for VNIR and SWIR?
        points = xml.findall("base/spatialCoverage/boundingPolygon/point")
        coords_xy = {p.find('frame').text: (float(p.find('longitude').text),
                                            float(p.find('latitude').text))
                     for p in points}
        coord_tags = ['upper_left', 'upper_right', 'lower_left', 'lower_right']
        self.lon_UL_UR_LL_LR = [coords_xy[ct][0] for ct in coord_tags]
        self.lat_UL_UR_LL_LR = [coords_xy[ct][1] for ct in coord_tags]

        # read the band related information: wavelength, fwhm
        self.nwvl = int(xml.find("product/image/%s/channels" % lbl).text)
        # FIXME hardcoded - DLR does not provide any smile information
        #   => smile coefficients are set to zero until we have valid ones
        self.nsmile_coef = 5
        self.smile_coef = np.zeros((self.nwvl, self.nsmile_coef), dtype=float)

        startband = 0 if self.detector_name == 'VNIR' else int(xml.find("product/image/vnir/channels").text)
        subset = slice(startband, startband + self.nwvl)
        bi = "specific/bandCharacterisation/bandID/"
        self.wvl_center = np.array([float(ele.text) for ele in xml.findall(bi + "wavelengthCenterOfBand")[subset]])
        self.fwhm = np.array([float(ele.text) for ele in xml.findall(bi + "FWHMOfBand")[subset]])
        self.gains = np.array([float(ele.text) for ele in xml.findall(bi + "GainOfBand")[subset]])
        self.offsets = np.array([float(ele.text) for ele in xml.findall(bi + "OffsetOfBand")[subset]])

        # read preview bands
        wvl_red = float(xml.find("product/image/%s/qlChannels/red" % lbl).text)
        wvl_green = float(xml.find("product/image/%s/qlChannels/green" % lbl).text)
        wvl_blue = float(xml.find("product/image/%s/qlChannels/blue" % lbl).text)
        self.preview_wvls = [wvl_red, wvl_green, wvl_blue]
        self.preview_bands = np.array([np.argmin(np.abs(self.wvl_center - wvl)) for wvl in self.preview_wvls])

        # read RPC coefficients
        for bID in xml.findall("product/navigation/RPC/bandID")[subset]:
            bN = 'band_%d' % np.int64(bID.attrib['number'])

            keys2combine = ('row_num', 'row_den', 'col_num', 'col_den')

            tmp = OrderedDict([(ele.tag.lower(), float(ele.text)) for ele in bID.findall('./')])
            self.rpc_coeffs[bN] = {k: v for k, v in tmp.items() if not k.startswith(keys2combine)}

            for n in keys2combine:
                self.rpc_coeffs[bN]['%s_coeffs' % n.lower()] = \
                    np.array([v for k, v in tmp.items() if k.startswith(n)])

        # drop bad bands from metadata if desired
        if self.cfg.drop_bad_bands:
            wvls = list(self.wvl_center)
            self.goodbands_inds = gbl = [wvls.index(wvl) for wvl in wvls
                                         if not (1358 < wvl < 1453 or
                                                 1814 < wvl < 1961)]

            if len(gbl) < len(wvls):
                for attrN in ['wvl_center', 'fwhm', 'offsets', 'gains', 'l_min', 'l_max']:
                    if getattr(self, attrN) is not None:
                        setattr(self, attrN, getattr(self, attrN)[gbl])

                self.nwvl = len(gbl)
                self.smile_coef = self.smile_coef[gbl, :]
                self.rpc_coeffs = OrderedDict([(k, v) for i, (k, v) in enumerate(self.rpc_coeffs.items())
                                               if i in gbl])

        # compute metadata derived from read data
        self.smile = self.calc_smile()
        self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
        self.solar_irrad = self.calc_solar_irradiance_CWL_FWHM_per_band()
        self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                          self.lat_UL_UR_LL_LR)), fix_invalid=True)
        from ...processors.spatial_transform import get_UTMEPSG_from_LonLat_cornersXY
        # NOTE:   self.cfg.target_epsg is set if user-provided or in case of Lon/Lat coordinates
        self.epsg_ortho = \
            self.cfg.target_epsg or \
            get_UTMEPSG_from_LonLat_cornersXY(lons=self.lon_UL_UR_LL_LR, lats=self.lat_UL_UR_LL_LR)

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

    def calc_snr_from_radiance(self, rad_data: GeoArray | np.ndarray, dir_snr_models: str):
        """Compute EnMAP SNR from radiance data for the given detector.

        SNR equation:    SNR = p0 + p1 * LTOA + p2 * LTOA ^ 2   [W/(m^2 sr nm)].

        NOTE:   The SNR model files (SNR_D1.bsq/SNR_D2.bsq) contain polynomial coefficients needed to compute SNR.

                SNR_D1.bsq: SNR model for EnMAP VNIR detector (contains high gain and low gain model coefficients)

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
        rad_data = np.array(rad_data)
        assert self.unitcode == 'TOARad'
        self.logger.info(f"Computing SNR from {self.detector_name} TOA radiance.")

        if self.detector_name == 'VNIR':
            gA = self._get_snr_model(dir_snr_models)

            coeffs_highgain = gA[0:3, :, :]  # [3 x ncols x nbands]
            coeffs_lowgain = gA[3:6, :, :]  # [3 x ncols x nbands]
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
            gA = self._get_snr_model(dir_snr_models)
            coeffs = gA[:]
            self.snr = coeffs[0, :, :] + coeffs[1, :, :] * rad_data[:, :, :] + coeffs[2, :, :] * rad_data[:, :, :] ** 2

    def _get_snr_model(self, dir_snr_models: str) -> GeoArray:
        """Get the SNR model coefficients for the current detector.

        NOTE: Missing bands are linearly interpolated.

        :param dir_snr_models:  directory containing the SNR models
        """
        path_snr_model = os.path.join(dir_snr_models, "SNR_D1.bsq" if self.detector_name == 'VNIR' else "SNR_D2.bsq")
        gA = GeoArray(path_snr_model)

        if gA.bands == self.nwvl:
            return gA

        else:
            from scipy.interpolate import make_interp_spline
            # equivalent to legacy scipy.interpolate.interp1d
            # interp1d(wvl, gA[:], axis=2, kind='linear', fill_value="extrapolate", bounds_error=False)(self.wvl_center)
            coeffs_interp = make_interp_spline(gA.meta.band_meta['wavelength'], gA[:], k=1, axis=2)(self.wvl_center)
            gA_ = GeoArray(coeffs_interp)
            gA_.meta.band_meta['wavelength'] = self.wvl_center

            return gA_

    def compute_geolayer_for_cube(self, elevation: GeoArray | float):
        """
        Compute the geolayer for the current detector.

        :param elevation:  elevation to be used
                           (DEM in map geometry or single value in meter above sea level)
        """
        self.logger.info('Computing %s geolayer...' % self.detector_name)
        GeolayerGen = \
            RPC_3D_Geolayer_Generator(
                rpc_coeffs_per_band=self.rpc_coeffs,
                elevation=elevation,
                enmapIm_cornerCoords=tuple(zip(self.lon_UL_UR_LL_LR, self.lat_UL_UR_LL_LR)),
                enmapIm_dims_sensorgeo=(self.nrows, self.ncols),
                CPUs=self.cfg.CPUs
            )
        lons, lats = GeolayerGen.compute_geolayer()

        self.geolayer_has_keystone = len(GeolayerGen.bandgroups_with_unique_rpc_coeffs) > 1

        return lons, lats

    def calc_solar_irradiance_CWL_FWHM_per_band(self) -> np.ndarray:
        from ...io.reader import read_solar_irradiance

        self.logger.debug('Calculating solar irradiance...')

        sol_irr = read_solar_irradiance(path_solar_irr_model=self.cfg.path_solar_irr, wvl_min_nm=350, wvl_max_nm=2500)

        irr_bands = []
        for band in self.srf.bands:
            WVL_band = self.srf.srfs_wvl if self.srf.wvl_unit == 'nanometers' else self.srf.srfs_wvl * 1000
            RSP_band = self.srf.srfs_norm01[band]
            sol_irr_at_WVL = np.interp(WVL_band, sol_irr[:, 0], sol_irr[:, 1], left=0, right=0)

            irr_bands.append(np.round(np.sum(sol_irr_at_WVL * RSP_band) / np.sum(RSP_band), 2))

        return np.array(irr_bands)


class EnMAP_Metadata_L1B_SensorGeo(object):
    """EnMAP Metadata class for the metadata of the complete EnMAP L1B product in sensor geometry incl. VNIR and SWIR.

    Attributes:
        - logger(logging.Logger):  None or logging instance
        - observation_datetime(datetime.datetime):  datetime of observation time
        - geom_view_zenith: viewing zenith angle
        - geom_view_azimuth: viewing azimuth angle
        - geom_sun_zenith: sun zenith angle
        - geom_sun_azimuth: sun azimuth angle
        - avg_elevation: average elevation of the scene in meters above sea level
        - mu_sun: needed by SICOR for TOARad > TOARef conversion
        - vnir(EnMAP_Metadata_VNIR_SensorGeo)
        - swir(EnMAP_Metadata_SWIR_SensorGeo)
        - detector_attrNames: attribute names of the detector objects
    """

    def __init__(self, path_metaxml, config: EnPTConfig, logger=None):
        """Get a metadata object instance for the given EnMAP L1B product in sensor geometry.

        :param path_metaxml:    file path of the EnMAP L1B metadata XML file
        :para, config:          EnPT configuration object
        :param logger:          instance of logging.logger or subclassed
        """
        self.cfg = config
        self.logger = logger or logging.getLogger()
        self.path_xml = path_metaxml
        self.rootdir = os.path.dirname(path_metaxml)

        # defaults - Common
        self.proc_level: str | None = None   # Dataset processing level
        self.version_provider: str = ''  # version of ground segment processing system
        self.observation_datetime: datetime | None = None  # Date and Time of image observation
        self.geom_view_zenith: float | None = None  # viewing zenith angle
        self.geom_view_azimuth: float | None = None  # viewing azimuth angle
        self.geom_sun_zenith: float | None = None  # sun zenith angle
        self.geom_sun_azimuth: float | None = None   # sun azimuth angle
        self.geom_angles_all: dict | None = None  # all view and sun angles available
        self.avg_elevation: float | None = None  # average elevation of the scene in meters above sea level
        self.mu_sun: float | None = None   # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist: float | None = None  # earth-sun distance
        self.aot: float | None = None  # scene aerosol optical thickness
        self.water_vapour: float | None = None  # scene water vapor [cm]
        self.vnir: EnMAP_Metadata_L1B_Detector_SensorGeo | None = None  # metadata of VNIR only
        self.swir: EnMAP_Metadata_L1B_Detector_SensorGeo | None = None  # metadata of SWIR only
        self.detector_attrNames: list = ['vnir', 'swir']  # attribute names of the detector objects
        self.filename_metaxml: str | None = None  # filename of XML metadata file

        self._scene_basename: str | None = None  # basename of the EnMAP image

    @property
    def scene_basename(self):
        if self.vnir:
            self._scene_basename = self.vnir.scene_basename

        return self._scene_basename

    def read_common_meta(self, path_xml):
        """Read the common metadata, principally stored in General Info.

        - the acquisition time
        - the geometrical observation and illumination

        :param path_xml: path to the main XML file
        :return: None
        """
        # load the metadata xml file
        xml = ElementTree.parse(path_xml).getroot()

        self.filename_metaxml = os.path.basename(path_xml)

        # read processing level
        self.proc_level = xml.find("base/level").text
        if self.proc_level != 'L1B':
            raise RuntimeError(self.proc_level, "Unexpected input data processing level. Expected 'L1B'.")

        # read version of ground segment processing system
        self.version_provider = xml.find("base/revision").text

        # raise a warning in case of old processing version (de-striping was implemented in version 01.02.00)
        if parse_version(self.version_provider) < parse_version('01.02.00'):
            self.logger.warning(
                f"The input EnMAP Level-1B image was processed with an old version of the ground segment "
                f"processing system (version {self.version_provider}), which, e.g. did not include de-striping. "
                f"It is highly recommended to re-download the dataset in the latest processing version from the "
                f"archive via the EOWEB GeoPortal (www.eoweb.dlr.de) before passing it to EnPT.")

        # read the acquisition time
        self.observation_datetime = \
            datetime.strptime(xml.find("base/temporalCoverage/startTime").text, '%Y-%m-%dT%H:%M:%S.%fZ')

        # get the distance earth sun from the acquisition date
        self.earthSunDist = self.get_earth_sun_distance(self.observation_datetime)

        # read Geometry (observation/illumination) angle
        # NOTE: EnMAP metadata provide also the angles for the image corners
        #       -> would allow even more precise computation (e.g., specific/sunElevationAngle/upper_left)
        # NOTE: alongOffNadirAngle is always near 0 and therefore ignored here (not relevant for AC)
        # FIXME VZA may be negative in DLR L1B data -> correct to always use the absolute value for SICOR?
        self.geom_view_zenith = np.abs(float(xml.find("specific/acrossOffNadirAngle/center").text))
        # FIXME correct to directly use sceneAzimuthAngle (14.24 (DLR) vs. 101.1 (AlpineTest)
        self.geom_view_azimuth = float(xml.find("specific/sceneAzimuthAngle/center").text)
        self.geom_sun_zenith = 90 - float(xml.find("specific/sunElevationAngle/center").text)
        self.geom_sun_azimuth = float(xml.find("specific/sunAzimuthAngle/center").text)
        self.avg_elevation = int(float(xml.find("specific/meanGroundElevation").text))
        self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))
        self.aot = float(xml.find("specific/qualityFlag/sceneAOT").text) / 1000  # scale factor is 1000
        self.water_vapour = float(xml.find("specific/qualityFlag/sceneWV").text) / 1000  # scale factor is 1000

        # TODO: revise this later to get rid of the duplicates with self.geom_xxx
        self.geom_angles_all = dict(
            view_zenith={e.tag: abs(float(e.text)) for e in xml.findall("specific/acrossOffNadirAngle/")},
            view_azimuth={e.tag: float(e.text) for e in xml.findall("specific/sceneAzimuthAngle/")},
            sun_zenith={e.tag: 90 - float(e.text)for e in xml.findall("specific/sunElevationAngle/")},
            sun_azimuth={e.tag: float(e.text) for e in xml.findall("specific/sunAzimuthAngle/")}
        )

    def get_earth_sun_distance(self, acqDate: datetime):
        """Get earth-sun distance (requires file of pre-calculated earth sun distance per day).

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

    def read_metadata(self):
        """Read the metadata of the entire EnMAP L1B product in sensor geometry."""
        # first read common metadata
        self.read_common_meta(self.path_xml)

        # define and read the VNIR metadata
        self.vnir = EnMAP_Metadata_L1B_Detector_SensorGeo('VNIR', config=self.cfg, logger=self.logger)
        self.vnir.read_metadata(self.path_xml)

        # define and read the SWIR metadata
        self.swir = EnMAP_Metadata_L1B_Detector_SensorGeo('SWIR', config=self.cfg, logger=self.logger)
        self.swir.read_metadata(self.path_xml)

    def to_XML(self) -> str:
        """Generate an XML metadata string from the L1B metadata."""
        xml = ElementTree.parse(self.path_xml).getroot()

        for lbl, detMeta in zip(['vnir', 'swir'], [self.vnir, self.swir]):
            xml.find("product/image/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
            xml.find("product/image/%s/dimension/columns" % lbl).text = str(detMeta.ncols)
            xml.find("product/quicklook/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
            xml.find("product/quicklook/%s/dimension/columns" % lbl).text = str(detMeta.ncols)

        xml_string = ElementTree.tostring(xml, encoding='unicode', pretty_print=True)

        return xml_string
