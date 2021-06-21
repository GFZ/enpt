# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
from ...options.config import EnPTConfig
from ..srf import SRF
from ...version import __version__

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class EnMAP_Metadata_L2A_MapGeo(object):
    def __init__(self,
                 config: EnPTConfig,
                 meta_l1b: EnMAP_Metadata_L1B_SensorGeo,
                 wvls_l2a: Union[List, np.ndarray],
                 dims_mapgeo: Tuple[int, int, int],
                 grid_res_l2a: Tuple[float, float],
                 logger=None):
        """EnMAP Metadata class for the metadata of the complete EnMAP L2A product in map geometry incl. VNIR and SWIR.

        :param config:              EnPT configuration object
        :param meta_l1b:            metadata object of the L1B dataset in sensor geometry
        :param wvls_l2a:            list of center wavelengths included in the L2A product
        :param dims_mapgeo:         dimensions of the EnMAP raster data in map geometry, e.g., (1024, 1000, 218)
        :param grid_res_l2a:        Coordinate grid resolution of the L2A product (x, y)
        :param logger:              instance of logging.logger or subclassed
        """
        self.cfg = config
        self._meta_l1b = meta_l1b
        self.grid_res = grid_res_l2a
        self.logger = logger or logging.getLogger()

        # defaults
        self.band_means: Optional[np.ndarray] = None  # band-wise means in unscaled values (percent for reflectance)
        self.band_stds: Optional[np.ndarray] = None  # band-wise standard deviations in unscaled values
        self.fileinfos: list = []  # file informations for each file beloning to the EnMAP L2A product

        # input validation
        if not len(wvls_l2a) == dims_mapgeo[2]:
            raise ValueError("The list of center wavelength must be qual to the number of bands in the map geometry "
                             "dimensions.")

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
        file_ext_l1b = os.path.splitext(meta_l1b.vnir.filename_data)[1]
        file_ext_l2a = \
            '.TIF' if self.cfg.output_format == 'GTiff' else \
            '.bsq' if self.cfg.output_interleave == 'band' else \
            '.bil' if self.cfg.output_interleave == 'line' else \
            '.bip'

        def convert_fn(fn):
            return fn.replace('L1B-', 'L2A-').replace(file_ext_l1b, file_ext_l2a)

        if not self.cfg.is_dummy_dataformat:
            self.scene_basename = convert_fn(meta_l1b.vnir.filename_data.split('-SPECTRAL_IMAGE')[0])
        else:
            self.scene_basename = os.path.splitext(meta_l1b.vnir.filename_data)[0]
        self.filename_data = convert_fn(meta_l1b.vnir.filename_data).replace('_VNIR', '')
        self.filename_deadpixelmap_vnir = convert_fn(meta_l1b.vnir.filename_deadpixelmap)
        self.filename_deadpixelmap_swir = convert_fn(meta_l1b.swir.filename_deadpixelmap)
        self.filename_quicklook_vnir = convert_fn(meta_l1b.vnir.filename_quicklook)
        self.filename_quicklook_swir = convert_fn(meta_l1b.swir.filename_quicklook)
        self.filename_mask_landwater = convert_fn(meta_l1b.vnir.filename_mask_landwater)
        self.filename_mask_clouds = convert_fn(meta_l1b.vnir.filename_mask_clouds)
        self.filename_mask_cloudshadow = convert_fn(meta_l1b.vnir.filename_mask_cloudshadow)
        self.filename_mask_haze = convert_fn(meta_l1b.vnir.filename_mask_haze)
        self.filename_mask_snow = convert_fn(meta_l1b.vnir.filename_mask_snow)
        self.filename_mask_cirrus = convert_fn(meta_l1b.vnir.filename_mask_cirrus)
        self.filename_metaxml = convert_fn(meta_l1b.filename_metaxml)

        # fuse band-wise metadata (sort all band-wise metadata by wavelengths but band number keeps as it is)
        # get band index order
        wvls_vswir = np.hstack([self._meta_l1b.vnir.wvl_center,
                                self._meta_l1b.swir.wvl_center])
        bandidx_order = np.array([np.argmin(np.abs(wvls_vswir - cwl)) for cwl in wvls_l2a])

        self.wvl_center = np.array(wvls_l2a)
        self.fwhm = np.hstack([meta_l1b.vnir.fwhm, meta_l1b.swir.fwhm])[bandidx_order]
        self.gains = np.full((dims_mapgeo[2],), 1. / self.cfg.scale_factor_boa_ref)  # => value range 0-1
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
            all_rpc_coeffs = OrderedDict(list(meta_l1b.vnir.rpc_coeffs.items()) +
                                         list(meta_l1b.swir.rpc_coeffs.items()))
            self.rpc_coeffs = OrderedDict([(k, v) for i, (k, v) in enumerate(all_rpc_coeffs.items())
                                           if i in bandidx_order])
        else:
            self.rpc_coeffs = OrderedDict()

        self.nrows = dims_mapgeo[0]
        self.ncols = dims_mapgeo[1]
        self.nwvl = dims_mapgeo[2]
        assert self.nwvl == len(self.wvl_center)

        common_UL_UR_LL_LR = self.get_common_UL_UR_LL_LR()
        self.lon_UL_UR_LL_LR = [lon for lon, lat in common_UL_UR_LL_LR]
        self.lat_UL_UR_LL_LR = [lat for lon, lat in common_UL_UR_LL_LR]
        self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                          self.lat_UL_UR_LL_LR)), fix_invalid=True)
        self.epsg = self._meta_l1b.vnir.epsg_ortho

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
        # NOTE:  Multiply by gains to get reflectance in the range 0-1
        data = datastack_vnir_swir[datastack_vnir_swir.mask_nodata[:]]
        self.band_means = np.mean(data, axis=0) * self.gains
        self.band_stds = np.std(data, axis=0) * self.gains

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
                                      '.BSQ', '.bsq',
                                      '.BIL', '.bil',
                                      '.BIP', '.bip',
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
        """Generate an XML metadata string from the L2A metadata."""
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

        def uk(path, value):
            xml.find(path).text = str(value)

        def create_SubElement(_parent, _tag, attrib={}, _text=None, nsmap=None, **_extra):
            result = ElementTree.SubElement(_parent, _tag, attrib, nsmap, **_extra)
            result.text = _text
            return result

        ############
        # metadata #
        ############

        uk("metadata/schema/processingLevel", self.proc_level)
        uk("metadata/name", self.filename_metaxml)
        # uk("metadata/comment").text = 'EnMAP Level 0 Product of datatake 987'  # FIXME hardcoded

        ##############
        # processing #
        ##############

        uk("processing/terrainCorrection", 'Yes')  # FIXME hardcoded {Yes, No}
        uk("processing/ozoneValue", 'NA')  # FIXME {[200-500], NA})
        uk("processing/season", 'NA')  # FIXME {summer, winter, NA}
        uk("processing/productFormat", 'GeoTIFF+Metadata')  # FIXME hardcoded
        # {BSQ+Metadata, BIL+Metadata, BIP+Metadata, JPEG2000+Metadata, GeoTiff+Metadata}
        uk("processing/mapProjection", 'UTM_Zone_of_Scene_Center')  # FIXME hardcoded
        # {UTM_Zone_of_Scene_Center, UTM_Zone_of_Scene_Center(-1), UTM_Zone_of_Scene_Center(+1),
        #  UTM_Zone_of_Datatake_Center, Geographic, European_Projection_LAEA, NA}
        uk("processing/DEMDBVersion", 'SRTM-C_v4')  # FIXME hardcoded
        # {SRTM-C-X_vv.rr, best-of-DEM_vv.rr, DEM-derivedfrom-Tandem-X_vv.rr, ASTER-GDEM_vv.rr, NA}
        uk("processing/correctionType", 'NA')  # FIXME hardcoded {Combined, Land_Mode, Water_Mode, NA}
        uk("processing/cirrusHazeRemoval", 'NA')  # FIXME hardcoded {Yes, No}
        uk("processing/bandInterpolation", 'NA')  # FIXME hardcoded {Yes, No}
        uk("processing/waterType", 'NA')  # FIXME hardcoded {Clear, Turbid, Highly_Turbid, NA}

        ########
        # base #
        ########

        # TODO update corner coordinates? DLR just uses the same coordinates like in L1B
        # xml.find("base/spatialCoverage" % lbl).text =
        uk("base/format", 'ENMAP_%s' % self.proc_level)
        uk("base/level", self.proc_level)
        uk("base/size", 'NA')  # FIXME Size of product. Attribute unit {byte, Kbyte, Mbyte, Gbyte}

        ############
        # specific #
        ############

        uk("specific/code", self.proc_level)

        # delete old band characterisation (as it contains a different set of bands)
        bChar_root = xml.find('specific/bandCharacterisation')
        bChar_root.clear()

        # recreate sub-elements for bandCharacterisation with respect to the current set of bands
        for i in range(self.nwvl):
            sub = ElementTree.SubElement(bChar_root, 'bandID', number=str(i + 1))
            for k, val in zip(['wavelengthCenterOfBand', 'FWHMOfBand', 'GainOfBand', 'OffsetOfBand'],
                              [self.wvl_center[i], self.fwhm[i], self.gains[i], self.offsets[i]]):
                ele = ElementTree.SubElement(sub, k)
                ele.text = str(val)

        ###########
        # product #
        ###########

        if not self.fileinfos:
            raise ValueError('Product file informations must be added before writing metadata. '
                             'Call add_product_fileinformation() before!')

        # image #
        #########

        # recreate sub-elements for image (L1B XML contains sub-elements for VNIR and SWIR here, L2A is merged)
        im_root = xml.find('product/image')
        im_root.clear()

        merge = ElementTree.SubElement(im_root, 'merge')
        for subEleArgs in [
            # tag, attribute dictionary, text
            ('channels', {}, str(self.nwvl)),
            ('name', {}, self.filename_data),
            ('size', {'unit': "Kbyte"}, str([i for i in self.fileinfos if i['name'] == self.filename_data][0]['size'])),
            ('version', {}, __version__),  # we use the EnPT version here
            ('format', {}, 'binary'),
            ('dimension',),
            ('dimensionGeographic',),
            ('qlChannelsSWIR',),
            ('qlChannelsVNIR',),
        ]:
            create_SubElement(merge, *subEleArgs)

        dim_root = xml.find('product/image/merge/dimension')
        for subEleArgs in [
            # tag, attribute dictionary, text
            ('columns', {}, str(self.ncols)),
            ('rows', {}, str(self.nrows)),
        ]:
            create_SubElement(dim_root, *subEleArgs)

        dimgeo_root = xml.find('product/image/merge/dimensionGeographic')
        for subEleArgs in [
            # tag, attribute dictionary, text
            ('longitude', {'unit': 'DEG'}, 'NA'),  # FIXME
            ('latitude', {'unit': 'DEG'}, 'NA'),  # FIXME 0.3314405
        ]:
            create_SubElement(dimgeo_root, *subEleArgs)

        qlswir_root = xml.find('product/image/merge/qlChannelsSWIR')
        for subEleArgs in [
            # tag, attribute dictionary, text
            ('red', {}, str(self.preview_wvls_vnir[0])),
            ('green', {}, str(self.preview_wvls_vnir[1])),
            ('blue', {}, str(self.preview_wvls_vnir[2])),
        ]:
            create_SubElement(qlswir_root, *subEleArgs)

        qlvnir_root = xml.find('product/image/merge/qlChannelsVNIR')
        for subEleArgs in [
            # tag, attribute dictionary, text
            ('red', {}, str(self.preview_wvls_vnir[0])),
            ('green', {}, str(self.preview_wvls_vnir[1])),
            ('blue', {}, str(self.preview_wvls_vnir[2])),
        ]:
            create_SubElement(qlvnir_root, *subEleArgs)

        # quicklook #
        #############

        from . import L2A_product_props_DLR
        for detName, detMetaL1B in zip(['VNIR', 'SWIR'], [self._meta_l1b.vnir, self._meta_l1b.swir]):
            lbl = L2A_product_props_DLR['xml_detector_label'][detName]

            fN_quicklook = self.filename_quicklook_vnir if detName == 'VNIR' else self.filename_quicklook_swir
            size_quicklook = [F['size'] for F in self.fileinfos
                              if os.path.splitext(F['name'])[0].endswith('-QL_%s' % detName)][0]
            uk("product/quicklook/%s/name" % lbl, fN_quicklook)
            uk("product/quicklook/%s/size" % lbl, str(size_quicklook))
            uk("product/quicklook/%s/dimension/rows" % lbl, str(self.nrows))
            uk("product/quicklook/%s/dimension/columns" % lbl, str(self.ncols))
            # uk("product/quicklook/%s/dimension/dimensionGeographic/longitude" % lbl, 'NA')  # FIXME
            # uk("product/quicklook/%s/dimension/dimensionGeographic/latitude" % lbl, 'NA')  # FIXME

        # productFileInformation #
        ##########################

        # get L1B product file information
        l1b_fileinfos = xmlSubtree2dict(xml, 'product/productFileInformation/')

        # clear old L1B file information in XML
        pFI_root = xml.find('product/productFileInformation')
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

        # ortho #
        #########

        if self.epsg == 4326:
            proj_str = 'Geographic'
        elif self.epsg == 3035:
            proj_str = 'LAEA-ETRS89'
        elif len(str(self.epsg)) == 5 and str(self.epsg)[:3] in ['326', '327']:
            zone = int(str(self.epsg)[-2:])
            if not 0 <= zone <= 60:
                raise RuntimeError('Invalid L2A UTM zone: %d.' % zone)
            proj_str = 'UTM_Zone%d_North' % zone if str(self.epsg).startswith('326') else 'UTM_Zone%d_South' % zone
        else:
            proj_str = 'NA'

        uk('product/ortho/projection', proj_str)  # {UTM_ZoneX_North, UTM_ZoneX_South (where X in {1..60}), Geographic, LAEA-ETRS89, NA}  # noqa
        uk('product/ortho/resolution', self.grid_res[0])
        uk('product/ortho/resampling', self.cfg.ortho_resampAlg)

        # band statistics
        #################

        if self.band_means is None or self.band_stds is None:
            raise ValueError('Band statistics have not yet been computed. Compute them first by calling '
                             'add_band_statistics()!')

        # delete old L1B band statistics
        bStats_root = xml.find('product/bandStatistics')
        bStats_root.clear()

        # recreate sub-elements for bandStatistics with respect to the current set of bands
        for i in range(self.nwvl):
            sub = ElementTree.SubElement(bStats_root, 'bandID', number=str(i + 1))
            for k, val in zip(['wavelength', 'mean', 'stdDeviation'],
                              [self.wvl_center[i], self.band_means[i], self.band_stds[i]]):
                ele = ElementTree.SubElement(sub, k)
                ele.text = str(val)

        #######################
        # generate XML string #
        #######################

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
