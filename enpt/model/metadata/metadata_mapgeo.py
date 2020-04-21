# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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

"""EnPT metadata objects for EnMAP data in map geometry."""

from datetime import datetime
from lxml import etree as ElementTree
import logging
import os
import fnmatch
from typing import Union, List, Tuple, Optional  # noqa: F401
from collections import OrderedDict
import numpy as np
from py_tools_ds.geo.vector.topology import Polygon, get_footprint_polygon  # noqa: F401  # flake8 issue
from geoarray import GeoArray

from .metadata_sensorgeo import EnMAP_Metadata_L1B_SensorGeo
from ...options.config import EnPTConfig, enmap_xres
from ..srf import SRF

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class EnMAP_Metadata_L2A_MapGeo(object):
    def __init__(self,
                 config: EnPTConfig,
                 meta_l1b: EnMAP_Metadata_L1B_SensorGeo,
                 wvls_l2a: Union[List, np.ndarray],
                 dims_mapgeo: Tuple[int, int, int],
                 logger=None):
        """EnMAP Metadata class for the metadata of the complete EnMAP L2A product in map geometry incl. VNIR and SWIR.

        :param config:              EnPT configuration object
        :param meta_l1b:            metadata object of the L1B dataset in sensor geometry
        :param wvls_l2a:            list of center wavelengths included in the L2A product
        :param dims_mapgeo:         dimensions of the EnMAP raster data in map geometry, e.g., (1024, 1000, 218)
        :param logger:              instance of logging.logger or subclassed
        """
        self.cfg = config
        self._meta_l1b = meta_l1b
        self.logger = logger or logging.getLogger()

        # defaults
        self.band_means: Optional[np.ndarray] = None  # band-wise means in unscaled values (percent for reflectance)
        self.band_stds: Optional[np.ndarray] = None  # band-wise standard deviations in unscaled values
        self.fileinfos: list = []  # file informations for each file beloning to the EnMAP L2A product

        self.proc_level = 'L2A'
        self.observation_datetime: datetime = meta_l1b.observation_datetime  # Date and Time of observation
        # FIXME VZA may be negative in DLR data
        self.geom_view_zenith: float = meta_l1b.geom_view_zenith  # viewing zenith angle
        self.geom_view_azimuth: float = meta_l1b.geom_view_azimuth  # viewing azimuth angle
        self.geom_sun_zenith: float = meta_l1b.geom_sun_zenith  # sun zenith angle
        self.geom_sun_azimuth: float = meta_l1b.geom_sun_azimuth  # sun azimuth angle
        self.mu_sun: float = meta_l1b.mu_sun  # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist: float = meta_l1b.earthSunDist  # earth-sun distance

        # generate file names for L2A output
        if not self.cfg.is_dummy_dataformat:
            self.scene_basename = meta_l1b.vnir.filename_data.split('-SPECTRAL_IMAGE')[0].replace('L1B-', 'L2A-')
        else:
            self.scene_basename = os.path.splitext(meta_l1b.vnir.filename_data)[0]
        self.filename_data = meta_l1b.vnir.filename_data.replace('L1B-', 'L2A-').replace('_VNIR', '')
        self.filename_deadpixelmap_vnir = meta_l1b.vnir.filename_deadpixelmap.replace('L1B-', 'L2A-')
        self.filename_deadpixelmap_swir = meta_l1b.swir.filename_deadpixelmap.replace('L1B-', 'L2A-')
        self.filename_quicklook_vnir = meta_l1b.vnir.filename_quicklook.replace('L1B-', 'L2A-')
        self.filename_quicklook_swir = meta_l1b.swir.filename_quicklook.replace('L1B-', 'L2A-')
        self.filename_mask_landwater = meta_l1b.vnir.filename_mask_landwater.replace('L1B-', 'L2A-')
        self.filename_mask_clouds = meta_l1b.vnir.filename_mask_clouds.replace('L1B-', 'L2A-')
        self.filename_mask_cloudshadow = meta_l1b.vnir.filename_mask_cloudshadow.replace('L1B-', 'L2A-')
        self.filename_mask_haze = meta_l1b.vnir.filename_mask_haze.replace('L1B-', 'L2A-')
        self.filename_mask_snow = meta_l1b.vnir.filename_mask_snow.replace('L1B-', 'L2A-')
        self.filename_mask_cirrus = meta_l1b.vnir.filename_mask_cirrus.replace('L1B-', 'L2A-')
        self.filename_metaxml = meta_l1b.filename_metaxml.replace('L1B-', 'L2A-')

        # fuse band-wise metadata (sort all band-wise metadata by wavelengths but band number keeps as it is)
        # get band index order
        wvls_sorted = np.array(sorted(np.hstack([self._meta_l1b.vnir.wvl_center,
                                                 self._meta_l1b.swir.wvl_center])))
        bandidx_order = np.array([np.argmin(np.abs(wvls_sorted - cwl)) for cwl in wvls_l2a])

        self.wvl_center = np.hstack([meta_l1b.vnir.wvl_center, meta_l1b.swir.wvl_center])[bandidx_order]
        self.fwhm = np.hstack([meta_l1b.vnir.fwhm, meta_l1b.swir.fwhm])[bandidx_order]
        self.gains = np.full((dims_mapgeo[2],), 100)  # implies reflectance scaled between 0 and 10000
        self.offsets = np.zeros((dims_mapgeo[2],))
        self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
        self.solar_irrad = np.hstack([meta_l1b.vnir.solar_irrad, meta_l1b.swir.solar_irrad])[bandidx_order]

        if not meta_l1b.vnir.nsmile_coef == meta_l1b.swir.nsmile_coef:
            raise ValueError('Expected equal number of smile coefficients for VNIR and SWIR. Received %d/%s.'
                             % (meta_l1b.vnir.nsmile_coef, meta_l1b.swir.nsmile_coef))

        self.nsmile_coef = meta_l1b.vnir.nsmile_coef
        self.smile_coef = np.vstack([meta_l1b.vnir.smile_coef, meta_l1b.swir.smile_coef])[bandidx_order, :]
        self.smile = np.hstack([meta_l1b.vnir.smile, meta_l1b.swir.smile])[:, bandidx_order]

        if not self.cfg.is_dummy_dataformat:
            self.rpc_coeffs = OrderedDict(zip(
                ['band_%d' % (i + 1) for i in range(dims_mapgeo[2])],
                [meta_l1b.vnir.rpc_coeffs['band_%d' % (i + 1)] if 'band_%d' % (i + 1) in meta_l1b.vnir.rpc_coeffs else
                 meta_l1b.swir.rpc_coeffs['band_%d' % (i + 1)] for i in bandidx_order]))
        else:
            self.rpc_coeffs = OrderedDict()

        self.nrows = dims_mapgeo[0]
        self.ncols = dims_mapgeo[1]
        self.nwvl = dims_mapgeo[2]
        common_UL_UR_LL_LR = self.get_common_UL_UR_LL_LR()
        self.lon_UL_UR_LL_LR = [lon for lon, lat in common_UL_UR_LL_LR]
        self.lat_UL_UR_LL_LR = [lat for lon, lat in common_UL_UR_LL_LR]
        self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                          self.lat_UL_UR_LL_LR)), fix_invalid=True)

        if meta_l1b.vnir.unit != meta_l1b.swir.unit or meta_l1b.vnir.unitcode != meta_l1b.swir.unitcode:
            raise RuntimeError('L2A data should have the same radiometric unit for VNIR and SWIR. '
                               'Received %s in %s for VNIR and %s in %s for SWIR.'
                               % (meta_l1b.vnir.unitcode, meta_l1b.vnir.unit,
                                  meta_l1b.swir.unitcode, meta_l1b.swir.unit))

        self.unit = meta_l1b.vnir.unit
        self.unitcode = meta_l1b.vnir.unitcode
        self.preview_wvls_vnir = meta_l1b.vnir.preview_wvls
        self.preview_wvls_swir = meta_l1b.swir.preview_wvls
        self.preview_bands_vnir = meta_l1b.vnir.preview_bands
        self.preview_bands_swir = np.array([np.argmin(np.abs(self.wvl_center - wvl))
                                            for wvl in meta_l1b.swir.preview_wvls])  # must index from VNIR band 0

        self.snr = None
        if meta_l1b.vnir.snr is not None:
            assert meta_l1b.swir.snr is not None
            self.snr = np.dstack([meta_l1b.vnir.snr, meta_l1b.swir.snr])[:, :, bandidx_order]

    def get_common_UL_UR_LL_LR(self):
        vnir_ulx, vnir_urx, vnir_llx, vnir_lrx = self._meta_l1b.vnir.lon_UL_UR_LL_LR
        vnir_uly, vnir_ury, vnir_lly, vnir_lry = self._meta_l1b.vnir.lat_UL_UR_LL_LR
        swir_ulx, swir_urx, swir_llx, swir_lrx = self._meta_l1b.swir.lon_UL_UR_LL_LR
        swir_uly, swir_ury, swir_lly, swir_lry = self._meta_l1b.swir.lat_UL_UR_LL_LR

        # use OUTER coordinates
        return ((min([vnir_ulx, swir_ulx]), max([vnir_uly, swir_uly])),
                (max([vnir_urx, swir_urx]), max([vnir_ury, swir_ury])),
                (min([vnir_llx, swir_llx]), min([vnir_lly, swir_lly])),
                (max([vnir_lrx, swir_lrx]), min([vnir_lry, swir_lry])))

    def add_band_statistics(self, datastack_vnir_swir: Union[np.ndarray, GeoArray]):
        R, C, B = datastack_vnir_swir.shape
        # NOTE:  DEVIDE by gains to reflectance in percent
        self.band_means = np.mean(datastack_vnir_swir.reshape(1, R * C, B), axis=1) / self.gains
        self.band_stds = np.mean(datastack_vnir_swir.reshape(1, R * C, B), axis=1) / self.gains

    def add_product_fileinformation(self, filepaths: List[str], sizes: List[int] = None, versions: List[str] = None):
        self.fileinfos = []

        for i, fp in enumerate(filepaths):
            ismeta = fp.endswith('METADATA.XML') or fp.endswith('_header.xml')  # FIXME
            if not os.path.exists(fp):
                if ismeta:
                    pass  # does not yet exist
                else:
                    raise FileNotFoundError(fp)

            ext = os.path.splitext(fp)[1]
            fmt = 'binary' if ext in ['.GEOTIFF',
                                      '.TIF',
                                      '.TIFF',
                                      '.GTIFF',
                                      '.BSQ',
                                      '.BIL',
                                      '.BIP',
                                      '.JPEG2000'] \
                else 'xml' if ext == '.XML' \
                else 'NA'
            fileinfo_dict = dict(
                name=os.path.basename(fp),
                size=sizes[i] if sizes else int(os.path.getsize(fp) / 1024) if not ismeta else '',
                version=versions[i] if versions else '',
                format=fmt
            )

            self.fileinfos.append(fileinfo_dict)

    def to_XML(self) -> str:
        """
        Generate an XML metadata string from the L2A metadata.
        """
        # use an XML parser that creates properly indented XML files even if new SubElements have been added
        parser = ElementTree.XMLParser(remove_blank_text=True)

        # parse (use L1B metadata as template)
        xml = ElementTree.parse(self._meta_l1b.path_xml, parser).getroot()

        if self.cfg.is_dummy_dataformat:
            self.logger.warning('No XML metadata conversion implemented for datasets different to the DLR format.'
                                'Metadata XML file will be empty.')
            return ''

        self.logger.warning('Currently, the L2A metadata XML file does not contain all relevant keys and contains '
                            'not updated values!')  # FIXME

        ############
        # metadata #
        ############

        xml.find("metadata/schema/processingLevel").text = self.proc_level
        xml.find("metadata/name").text = self.filename_metaxml
        # xml.find("metadata/comment").text = 'EnMAP Level 0 Product of datatake 987'  # FIXME hardcoded

        ##############
        # processing #
        ##############

        xml.find("processing/terrainCorrection").text = 'Yes'  # FIXME hardcoded {Yes, No}
        xml.find("processing/ozoneValue").text = 'NA'  # FIXME {[200-500], NA}
        xml.find("processing/season").text = 'NA'  # FIXME {summer, winter, NA}
        xml.find("processing/productFormat").text = 'GeoTIFF+Metadata'  # FIXME hardcoded
        # {BSQ+Metadata, BIL+Metadata, BIP+Metadata, JPEG2000+Metadata, GeoTiff+Metadata}
        xml.find("processing/mapProjection").text = 'UTM_Zone_of_Scene_Center'  # FIXME hardcoded
        # {UTM_Zone_of_Scene_Center, UTM_Zone_of_Scene_Center(-1), UTM_Zone_of_Scene_Center(+1),
        #  UTM_Zone_of_Datatake_Center, Geographic, European_Projection_LAEA, NA}
        xml.find("processing/DEMDBVersion").text = 'SRTM-C_v4'  # FIXME hardcoded
        # {SRTM-C-X_vv.rr, best-of-DEM_vv.rr, DEM-derivedfrom-Tandem-X_vv.rr, ASTER-GDEM_vv.rr, NA}
        xml.find("processing/correctionType").text = 'NA'  # FIXME hardcoded {Combined, Land_Mode, Water_Mode, NA}
        xml.find("processing/cirrusHazeRemoval").text = 'NA'  # FIXME hardcoded {Yes, No}
        xml.find("processing/bandInterpolation").text = 'NA'  # FIXME hardcoded {Yes, No}
        xml.find("processing/waterType").text = 'NA'  # FIXME hardcoded {Clear, Turbid, Highly_Turbid, NA}

        ########
        # base #
        ########

        # TODO update corner coordinates? DLR just uses the same coordinates like in L1B
        # xml.find("base/spatialCoverage" % lbl).text =
        xml.find("base/format").text = 'ENMAP_%s' % self.proc_level
        xml.find("base/level").text = self.proc_level
        xml.find("base/size").text = 'NA'  # FIXME Size of product. Attribute unit {byte, Kbyte, Mbyte, Gbyte}

        ############
        # specific #
        ############

        xml.find("specific/code").text = self.proc_level
        bi = "specific/bandCharacterisation/bandID/"
        for ele, gain in zip(xml.findall(bi + "GainOfBand"), self.gains):
            ele.text = str(gain)
        for ele, offset in zip(xml.findall(bi + "OffsetOfBand"), self.offsets):
            ele.text = str(offset)

        ###########
        # product #
        ###########

        if not self.fileinfos:
            raise ValueError('Product file informations must be added before writing metadata. '
                             'Call add_product_fileinformation() before!')

        from . import L2A_product_props_DLR
        for detName, detMetaL1B in zip(['VNIR', 'SWIR'], [self._meta_l1b.vnir, self._meta_l1b.swir]):
            lbl = L2A_product_props_DLR['xml_detector_label'][detName]
            # FIXME DLR uses L0 filenames for VNIR/SWIR separately?!
            xml.find("product/image/%s/name" % lbl).text = detMetaL1B.filename_data
            # FIXME this is the size of the VNIR/SWIR stack
            size = [F['size'] for F in self.fileinfos if os.path.splitext(F['name'])[0].endswith('-SPECTRAL_IMAGE')][0]
            xml.find("product/image/%s/size" % lbl).text = str(size)
            # FIXME DLR data dimensions equal neither L2A data nor L1B data
            xml.find("product/image/%s/channels" % lbl).text = str(detMetaL1B.nwvl)
            xml.find("product/image/%s/dimension/rows" % lbl).text = str(self.nrows)
            xml.find("product/image/%s/dimension/columns" % lbl).text = str(self.ncols)
            # xml.find("product/image/%s/dimension/dimensionGeographic/longitude" % lbl).text = 'NA'  # TODO
            # xml.find("product/image/%s/dimension/dimensionGeographic/latitude" % lbl).text = 'NA'

            fN_quicklook = self.filename_quicklook_vnir if detName == 'VNIR' else self.filename_quicklook_swir
            size_quicklook = [F['size'] for F in self.fileinfos
                              if os.path.splitext(F['name'])[0].endswith('-QL_%s' % detName)][0]
            xml.find("product/quicklook/%s/name" % lbl).text = fN_quicklook
            xml.find("product/quicklook/%s/size" % lbl).text = str(size_quicklook)
            xml.find("product/quicklook/%s/dimension/rows" % lbl).text = str(self.nrows)
            xml.find("product/quicklook/%s/dimension/columns" % lbl).text = str(self.ncols)
            # xml.find("product/quicklook/%s/dimension/dimensionGeographic/longitude" % lbl).text = 'NA'
            # xml.find("product/quicklook/%s/dimension/dimensionGeographic/latitude" % lbl).text = 'NA'

        # productFileInformation
        ########################

        # get L1B product file information
        l1b_fileinfos = xmlSubtree2dict(xml, 'product/productFileInformation/')

        # clear old L1B file information in XML
        pFI_root = xml.findall('product/productFileInformation')[0]
        pFI_root.clear()

        # recreate sub-elements for productFileInformation according to L2A file information
        for i, fileInfo in enumerate(self.fileinfos):
            fn_l1b_exp = fileInfo['name'].replace('L2A', '*').replace('-SPECTRAL_IMAGE', '-SPECTRAL_IMAGE_VNIR')
            l1b_fileInfo = [fI for fI in l1b_fileinfos.values() if fnmatch.fnmatch(fI['name'], fn_l1b_exp)]

            if l1b_fileInfo:
                # TODO update file size of METADATA.XML (has not been written yet)
                fileInfo['size'] = fileInfo['size'] or l1b_fileInfo[0]['size']
                fileInfo['version'] = fileInfo['version'] or l1b_fileInfo[0]['version']
            else:
                # FIXME if no L1B equivalent is found for the file to be written, the file version will be empty ('')
                pass

            sub = ElementTree.SubElement(pFI_root, 'file', number=str(i))

            for k, kw in zip(['name', 'size', 'version', 'format'], [{}, {'unit': 'kbytes'}, {}, {}]):
                ele = ElementTree.SubElement(sub, k, **kw)
                ele.text = str(fileInfo[k])

        # TODO update product/ortho/projection
        #      {UTM_ZoneX_North, UTM_ZoneX_South (where X in {1..60}), Geographic, LAEA-ETRS89, NA}
        xml.find('product/ortho/resolution').text = str(enmap_xres)
        xml.find('product/ortho/resampling').text = self.cfg.ortho_resampAlg

        # band statistics
        #################

        if self.band_means is None or self.band_stds is None:
            raise ValueError('Band statistics have not yet been computed. Compute them first by calling '
                             'add_band_statistics()!')

        bs = "specific/bandStatistics/bandID/"
        for ele, mean in zip(xml.findall(bs + "meanReflectance"), self.band_means):
            ele.text = str(mean)
        for ele, std in zip(xml.findall(bs + "stdDeviation"), self.band_stds):
            ele.text = str(std)

        xml_string = ElementTree.tostring(xml, encoding='unicode', pretty_print=True)

        return xml_string


def xmlSubtree2dict(xml_root, path_subtree) -> OrderedDict:
    outDict = OrderedDict()
    allEle = xml_root.findall(path_subtree)

    for ele in allEle:
        eleKey = '%s_%s' % (ele.tag, ele.get('number'))
        outDict[eleKey] = dict()
        for subele in ele:
            outDict[eleKey][subele.tag] = subele.text

    return outDict
