"""Handling of EnMAP Level-1B products."""

from xml.etree import ElementTree
from types import SimpleNamespace
from datetime import datetime
from os import path
import numpy as np
import spectral
from scipy.interpolate import interp2d


class EnMAPL1BProduct(object):
    """Reader for EnMAP Level-1B products.

    Attributes:
        - vnir
            - ...
        - swir
            - same as vor [vnir]
        - root_dir: path to root dir of EnMAP product

    """

    def __init__(self, header_fn: str, observation_time: datetime, lon_lat_smpl=(15, 15), nsmile_coef=4):
        """Level-1B product object for EnMAP data.

        :param header_fn: Filename of EnMAP Level-1B product XML file
        :param observation_time: datetime of observation time (currently missing in metadata)
        :param lon_lat_smpl: number if sampling points in lon, lat fields fields
        :param nsmile_coef: number of polynomial coefficients for smile
        """
        self.vnir = SimpleNamespace()
        self.swir = SimpleNamespace()
        self.header_fn = header_fn
        self.root_dir = path.dirname(self.header_fn)

        xml = ElementTree.parse(self.header_fn).getroot()

        for detector, detector_label in zip((self.vnir, self.swir), ("detector1", "detector2")):
            detector.fwhm = np.array(xml.findall("%s/fwhm" % detector_label)[0].text.replace("\n", "").split(),
                                     dtype=np.float)
            detector.wvl_center = np.array(
                xml.findall("%s/centre_wavelength" % detector_label)[0].text.replace("\n", "").split(), dtype=np.float)
            detector.nwvl = len(detector.wvl_center)
            detector.nrows = np.int(xml.findall("%s/rows" % detector_label)[0].text)
            detector.ncols = np.int(xml.findall("%s/columns" % detector_label)[0].text)
            detector.smile_coef = np.array(xml.findall("%s/smile" % detector_label)[0].text.replace("\n", "").split(),
                                           dtype=np.float
                                           ).reshape((-1, nsmile_coef + 1))[:, 1:]
            detector.nsmile_coef = nsmile_coef
            detector.l_min = np.array(xml.findall("%s/L_min" % detector_label)[0].text.split(), dtype=np.float)
            detector.l_max = np.array(xml.findall("%s/L_max" % detector_label)[0].text.split(), dtype=np.float)
            detector.geom_view_zenith = np.float(
                xml.findall("%s/observation_geometry/zenith_angle" % detector_label)[0].text.split()[0])
            detector.geom_view_azimuth = np.float(
                xml.findall("%s/observation_geometry/zenith_angle" % detector_label)[0].text.split()[0])
            detector.geom_illu_zenith = np.float(
                xml.findall("%s/illumination_geometry/zenith_angle" % detector_label)[0].text.split()[0])
            detector.geom_illu_azimuth = np.float(
                xml.findall("%s/illumination_geometry/zenith_angle" % detector_label)[0].text.split()[0])
            detector.mu_sun = np.cos(np.deg2rad(detector.geom_illu_zenith))
            detector.observation_time = observation_time

            # smile(icol, iwvl) = sum_(p=0)^(nsmile_coef-1) smile_coef[iwvl, p] * icol**p (1)
            detector.smile = np.inner(
                # the sum in (1) is expressed as inner product of two arrays with
                # inner dimension as the polynomial smile coefficients
                # shape = (ncols, nsmile_coef) of polynomial x
                np.array([[icol ** p for p in np.arange(detector.nsmile_coef)] for icol in np.arange(detector.ncols)]),
                detector.smile_coef  # shape = (nwvl, nsmile_coef)
            )  # shape = (ncols, nwvl)

            detector.lat_UL_UR_LL_LR = (
                float(xml.findall("%s/geometry/bounding_box/%s_northing" %
                                  (detector_label, corner))[0].text.split()[0]) for corner in ("UL", "UR", "LL", "LR"))
            detector.lon_UL_UR_LL_LR = (
                float(xml.findall("%s/geometry/bounding_box/%s_easting" %
                                  (detector_label, corner))[0].text.split()[0]) for corner in ("UL", "UR", "LL", "LR"))

            detector.lats = self.interpolate_corners(*detector.lat_UL_UR_LL_LR, *lon_lat_smpl)
            detector.lons = self.interpolate_corners(*detector.lon_UL_UR_LL_LR, *lon_lat_smpl)

            detector.data_fn = xml.findall("%s/filename" % detector_label)[0].text.split()[0]
            detector.data = spectral.open_image(path.join(
                self.root_dir, detector.data_fn.replace(".bsq", ".hdr")))[:, :, :]  # unit of [dn]

    @staticmethod
    def interpolate_corners(ul: float, ur: float, ll: float, lr: float, nx: int, ny: int):
        """Compute interpolated field from corner values of a scalar field given at: ul, ur, ll, lr.

        :param nx, ny: final shape
        """
        ff = interp2d(x=[0, 1], y=[0, 1], z=[
            [ul, ur],
            [ll, lr]], kind='linear')
        rr = np.zeros((nx, ny), dtype=np.float)
        for i, x in enumerate(np.linspace(0, 1, nx)):
            for j, y in enumerate(np.linspace(0, 1, ny)):
                rr[i, j] = ff(x, y)
        return rr


class NoMatchError(Exception):
    """Raised of given path in xml could not be matched."""

    pass


class MultipleMatchError(Exception):
    """Raised if more than one match in xml is found."""

    pass
