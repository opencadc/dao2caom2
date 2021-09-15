# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  : 4 $
#
# ***********************************************************************
#

import logging
import re

from caom2 import CalibrationLevel, DataProductType, TargetType, ProductType
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from dao2caom2.dao_name import  DAOName


__all__ = ['factory', 'get_current']


class Telescope:
    """
    This class is the Spectrum implementation.
    """
    def __init__(self):
        self._uri = None
        self._logger = logging.getLogger(self.__class__.__name__)

    def _get_wavelength(self, header):
        return mc.to_float(header.get('WAVELENG'))

    def get_artifact_product_type(self, header):
        obs_type = header.get('OBSTYPE')
        if obs_type == 'object':
            product_type = ProductType.SCIENCE
        else:
            product_type = ProductType.CALIBRATION
        return product_type

    def get_calibration_level(self):
        return CalibrationLevel.RAW_STANDARD

    def get_data_product_type(self):
        return DataProductType.SPECTRUM

    def get_energy_axis_function_delta(self, header):
        wavelength = self._get_wavelength(header)
        cdelt = header.get('DELTA_WL')
        if wavelength is None:
            cdelt = None
        else:
            if cdelt is None:
                dispersion = header.get('DISPERSI')
                dispaxis = self._get_dispaxis(header)
                if dispaxis == 1:
                    xbin = mc.to_float(header.get('XBIN'))
                else:
                    xbin = mc.to_float(header.get('YBIN'))
                cdelt = dispersion * 15.0 * xbin / 1000.0
        return cdelt

    def get_energy_axis_function_refcoord_pix(self, header):
        wavelength = self._get_wavelength(header)
        if wavelength is None:
            crpix = None
        else:
            crpix = header.get('REFPIXEL')
            if crpix is None:
                temp = header.get('DATASEC')
                if temp is not None:
                    datasec = re.sub(
                        r'(\[)(\d+:\d+,\d+:\d+)(\])', r'\g<2>', temp
                    )
                    (dx, dy) = datasec.split(',')
                    (xl, xh) = dx.split(':')
                    (yl, yh) = dy.split(':')
                    dispaxis = self._get_dispaxis(header)
                    if dispaxis == 1:
                        crpix = (int(xh) - int(xl)) / 2.0 + int(xl)
                    else:
                        crpix = (int(yh) - int(yl)) / 2.0 + int(yl)
        return crpix

    def get_energy_axis_function_naxis(self, header):
        return 1

    def get_geo(self):
        # DB 10-09-20
        # Google Maps to give you latitude/longitude if desired. 48.519497
        # and -123.416502.  Not sure of the elevation.
        return ac.get_location(48.519497, -123.416502, 210.0)

    def get_position_function_cd11(self, header):
        return (-1) * self.get_position_function_cd22(header)

    def get_position_function_cd12(self, header):
        return self.get_position_function_cd21(header)

    def get_position_function_cd21(self, header):
        obs_type = header.get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            result = 0.0
        return result

    def get_position_function_cd22(self, header):
        obs_type = header.get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = 0.001388
        return result

    def get_position_function_coord1_pix(self, header):
        return self.get_position_function_coord2_pix(header)

    def get_position_function_coord2_pix(self, header):
        obs_type = header.get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            result = 1.0
        return result

    def get_position_function_dimension_naxis1(self, ignore_header):
        return 1

    def get_position_function_dimension_naxis2(self, ignore_header):
        return 1

    def get_telescope_name(self):
        return None

    def get_target_type(self):
        return TargetType.OBJECT

    def get_time_axis_val(self, header):
        return ac.get_datetime(header.get('DATE-OBS'))

    @property
    def uri(self):
        return self._uri

    @uri.setter
    def uri(self, value):
        self._uri = value


class Dao12Metre(Telescope):

    def __init__(self):
        super().__init__()

    def get_geo(self):
        return ac.get_location(48.52092, -123.42006, 225.0)

    def get_telescope_name(self):
        return 'DAO 1.2-m'


class Dao18Metre(Telescope):

    def __init__(self):
        super().__init__()

    def get_geo(self):
        return ac.get_location(48.51967, -123.41833, 232.0)

    def get_telescope_name(self):
        return 'DAO 1.8-m'


class SkyCam(Telescope):

    def __init__(self):
        super().__init__()

    def get_telescope_name(self):
        return 'DAO Skycam'

    def get_time_axis_val(self, header):
        return ac.get_datetime(header.get('CLOCKVAL'))


class Imaging(Telescope):

    def __init__(self):
        super().__init__()

    def _get_position_by_scale_size_bin(self, header):
        result = None
        platescale = mc.to_float(header.get('PLTSCALE'))
        pixsize = mc.to_float(header.get('PIXSIZE'))
        xbin = mc.to_float(header.get('XBIN'))
        if platescale is not None and pixsize is not None and xbin is not None:
            result = platescale * pixsize * xbin / 3600000.0
        return result

    def get_data_product_type(self):
        return DataProductType.IMAGE

    def get_energy_axis_function_delta(self, header):
        return header.get('BANDPASS')

    def get_energy_axis_function_refcoord_pix(self, ignore_header):
        return 1.0

    def get_position_function_cd11(self, header):
        return self._get_position_by_scale_size_bin(header)

    def get_position_function_cd22(self, header):
        return self._get_position_by_scale_size_bin(header)

    def get_position_function_coord1_pix(self, header):
        return header.get('NAXIS1') / 2.0

    def get_position_function_coord2_pix(self, header):
        return header.get('NAXIS2') / 2.0

    def get_position_function_dimension_naxis1(self, header):
        return header.get('NAXIS1')

    def get_position_function_dimension_naxis2(self, header):
        return header.get('NAXIS2')

    def get_target_type(self):
        return TargetType.FIELD


class Processed(Telescope):

    def __init__(self):
        super().__init__()

    def get_calibration_level(self):
        return CalibrationLevel.CALIBRATED


class Dao12MetreImage(Dao12Metre, Imaging):

    def __init__(self):
        super().__init__()


class Dao12MetreProcessedImage(Dao12MetreImage, Processed):

    def __init__(self):
        super().__init__()

    def get_position_function_cd11(self, header):
        return header.get('CD1_1')

    def get_position_function_cd12(self, header):
        return header.get('CD1_2')

    def get_position_function_cd21(self, header):
        return header.get('CD2_1')

    def get_position_function_cd22(self, header):
        return header.get('CD2_2')

    def get_position_function_coord1_pix(self, header):
        return header.get('CRPIX1')

    def get_position_function_coord2_pix(self, header):
        return header.get('CRPIX2')


class Dao12MetreProcessedSpectrum(Dao12Metre, Processed):

    def __init__(self):
        super().__init__()

    def get_energy_axis_function_naxis(self, header):
        return header.get('NAXIS1')


class Dao12MetreSpectrum(Dao12Metre):

    def __init__(self):
        super().__init__()

    def _get_dispaxis(self, header):
        return header.get('DISPAXIS', 2)

    def get_energy_axis_function_naxis(self, header):
        dispaxis = self._get_dispaxis(header)
        return header.get(f'NAXIS{dispaxis}')


class Dao18MetreImage(Dao18Metre, Imaging):

    def __init__(self):
        super().__init__()


class Dao18MetreProcessedImage(Dao18MetreImage, Processed):

    def __init__(self):
        super().__init__()

    def get_position_function_cd11(self, header):
        return header.get('CD1_1')

    def get_position_function_cd12(self, header):
        return header.get('CD1_2')

    def get_position_function_cd21(self, header):
        return header.get('CD2_1')

    def get_position_function_cd22(self, header):
        return header.get('CD2_2')

    def get_position_function_coord1_pix(self, header):
        return header.get('CRPIX1')

    def get_position_function_coord2_pix(self, header):
        return header.get('CRPIX2')


class Dao18MetreProcessedSpectrum(Dao18Metre, Processed):

    def __init__(self):
        super().__init__()

    def get_energy_axis_function_naxis(self, header):
        return header.get('NAXIS1')


class Dao18MetreSpectrum(Dao18Metre):

    def __init__(self):
        super().__init__()

    def _get_dispaxis(self, header):
        return header.get('DISPAXIS', 1)

    def get_energy_axis_function_naxis(self, header):
        dispaxis = self._get_dispaxis(header)
        return header.get(f'NAXIS{dispaxis}')


# globally accessible pointer, placed here so they survive the
# importlib.import_module call from fits2caom2

# latest instances of Telescope
current = {}
# single instance of DefiningMetadataFinder
defining_metadata_finder = None


def factory(uri):
    # at this point, decompose the classes based on what can be determined
    # from the file name only
    ignore_scheme, ignore_path, f_name = mc.decompose_uri(uri)
    if f_name.startswith('a'):
        logging.error('skycam')
        result = SkyCam()
    else:
        defining_metadata = defining_metadata_finder.get(uri)
        logging.error(defining_metadata)
        f_id = DAOName.remove_extensions(f_name)
        if defining_metadata.data_product_type == DataProductType.IMAGE:
            if (
                re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122')
                ):
                    logging.error('1.2m processed image')
                    result = Dao12MetreProcessedImage()
                else:
                    logging.error('1.8m processed image')
                    result = Dao18MetreProcessedImage()
            elif (
                f_name.startswith('dao_c122') or f_name.startswith('dao_r122')
            ):
                logging.error('1.2m image')
                result = Dao12MetreImage()
            else:
                logging.error('1.8m image')
                result = Dao18MetreImage()
        else:
            if (
                re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122')
                ):
                    logging.error('1.2m processed spectrum')
                    result = Dao12MetreProcessedSpectrum()
                else:
                    logging.error('1.8m processed spectrum')
                    result = Dao18MetreProcessedSpectrum()
            elif (
                f_name.startswith('dao_c122') or f_name.startswith('dao_r122')
            ):
                logging.error('1.2m spectrum')
                result = Dao12MetreSpectrum()
            else:
                logging.error('1.8m spectrum')
                result = Dao18MetreSpectrum()

    result.uri = uri
    global current
    current[uri] = result


def get_current(uri):
    return current.get(uri)
