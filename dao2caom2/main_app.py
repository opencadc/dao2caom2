# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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
#  $Revision: 4 $
#
# ***********************************************************************
#

import importlib
import logging
import os
import re
import sys
import traceback

from astropy.coordinates import SkyCoord, FK5
import astropy.units as u

from caom2 import Observation, TargetType, DataProductType, ProductType
from caom2 import ObservationIntentType, CalibrationLevel
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc


__all__ = ['dao_main_app', 'update', 'DAOName', 'COLLECTION', 'APPLICATION',
           'to_caom2']


APPLICATION = 'dao2caom2'
COLLECTION = 'DAO'


class DAOName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed raw files and compressed product files in storage
    - uncompressed product files have an extension added to the input
      names - lower case if there's a single input, upper case if there's
      multiple inputs, as far as I can tell right now.
    """

    DAO_NAME_PATTERN = '*'

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None):
        self.fname_in_ad = file_name
        if obs_id is None:
            obs_id = DAOName.get_obs_id(file_name)
            file_id = mc.StorageName.remove_extensions(file_name)
        else:
            file_id = obs_id
        super(DAOName, self).__init__(
            obs_id, COLLECTION, DAOName.DAO_NAME_PATTERN, fname_on_disk)
        self._file_id = file_id
        self._logger = logging.getLogger(__name__)
        self._logger.error(self)

    def __str__(self):
        return f'obs id {self._obs_id} file name {self.file_name} ' \
               f'file id {self._file_id}'

    def is_valid(self):
        return True

    @staticmethod
    def is_processed(entry):
        # the entry is a uri
        file_id = mc.CaomName(entry).file_id
        result = False
        if (re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', file_id) or
                re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', file_id)):
            result = True
        return result

    @staticmethod
    def is_unprocessed_reticon(entry):
        # the entry is a uri
        file_id = mc.CaomName(entry).file_id
        result = False
        if re.match('dao_[r]\\d{3}_\\d{4}_\\d{6}', file_id):
            result = True
        return result

    @staticmethod
    def get_obs_id(file_name):
        # observation ID differs from the file ID for processed data, except
        # for composite processed observations (master biases and flats)
        file_id = mc.StorageName.remove_extensions(file_name)
        logging.error(f'file_id is {file_id}')
        if re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aev]', file_id):
            obsID = file_id[0:-2]
        else:
            obsID = file_id
        return obsID


def get_artifact_product_type(header):
    obs_type = _get_obs_type(header)
    if obs_type == 'object':
        product_type = ProductType.SCIENCE
    else:
        product_type = ProductType.CALIBRATION
    return product_type


def get_calibration_level(parameters):
    uri = parameters.get('uri')
    result = CalibrationLevel.RAW_STANDARD
    if DAOName.is_processed(uri):
        result = CalibrationLevel.CALIBRATED
    return result


def get_data_product_type(header):
    obs_mode = header.get('OBSMODE')
    if obs_mode is not None and obs_mode.strip() == 'Imaging':
        data_product_type = DataProductType.IMAGE
    else:
        data_product_type = DataProductType.SPECTRUM
    return data_product_type


def get_energy_function_naxis(parameters):
    uri = parameters.get('uri')
    header = parameters.get('header')
    naxis = 1
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        if DAOName.is_processed(uri):
            naxis = _get_naxis1(header)
        else:
            dispaxis = _get_dispaxis(header)
            logging.error(f'dispaxis in get energy function {dispaxis}')
            if dispaxis == 1:
                naxis = _get_naxis1(header)
            elif dispaxis == 2:
                naxis = _get_naxis2(header)
            else:
                raise mc.CadcException(
                    f'Could not find dispaxis for \'TODO\'')
    return naxis


def get_energy_function_delta(header):
    cdelt = 1.0
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        wavelength = _get_wavelength(header)
        cdelt = header.get('DELTA_WL')
        if wavelength is None:
            pass
        else:
            if cdelt is None:
                logging.error(f'hows the wavelength? xbin {header.get("XBIN")} '
                              f'ybin {header.get("YBIN")}')
                dispersion = header.get('DISPERSI')
                dispaxis = _get_dispaxis(header)
                if dispaxis == 1:
                    xbin = mc.to_float(header.get('XBIN'))
                else:
                    xbin = mc.to_float(header.get('YBIN'))
                cdelt = dispersion * 15.0 * xbin / 1000.0
                logging.error(f'cdelt is {cdelt} dispersion is '
                              f'{dispersion} xbin is {xbin}')
    else:
        cdelt = header.get('BANDPASS')
    return cdelt


def get_energy_function_pix(header):
    crpix = 1.0
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        wavelength = _get_wavelength(header)
        if wavelength is None:
            pass
        else:
            crpix = header.get('REFPIXEL')
            if crpix is None:
                temp = header.get('DATASEC')
                if temp is not None:
                    datasec = re.sub(r'(\[)(\d+:\d+,\d+:\d+)(\])', r'\g<2>', temp)
                    (dx, dy) = datasec.split(',')
                    (xl, xh) = dx.split(':')
                    (yl, yh) = dy.split(':')
                    dispaxis = _get_dispaxis(header)
                    if dispaxis == 1:
                        crpix = (int(xh) - int(xl)) / 2.0 + int(xl)
                    else:
                        crpix = (int(yh) - int(yl)) / 2.0 + int(yl)
    return crpix


def get_energy_resolving_power(header):
    resolving_power = None
    wavelength = _get_wavelength(header)
    if wavelength is not None:
        band_pass = mc.to_float(header.get('BANDPASS'))
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.IMAGE:
            resolving_power = wavelength / band_pass
        else:
            # assume 2.5 pixel wide resolution element
            if band_pass is None:
                cdelt = get_energy_function_delta(header)
                logging.error(f'no bandpadd cdelt {cdelt} wavln {wavelength}')
                resolving_power = wavelength / (2.5 * cdelt)
            else:
                logging.error(f'bandpass {band_pass} wavln {wavelength}')
                resolving_power = wavelength / (2.5 * band_pass)
    return resolving_power


def get_geo_x(header):
    x, ignore_y, ignore_z = _get_geo(header)
    return x


def get_geo_y(header):
    ignore_x, y, ignore_z = _get_geo(header)
    return y


def get_geo_z(header):
    ignore_x, ignore_y, z = _get_geo(header)
    return z


def get_members(header):
    pass


def get_obs_intent(header):
    obs_type = _get_obs_type(header)
    if obs_type == 'object':
        intent = ObservationIntentType.SCIENCE
    else:
        intent = ObservationIntentType.CALIBRATION
    return intent


def get_position_function_coord1_pix(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            result = _get_naxis1(header) / 2.0
    else:
        if data_product_type == DataProductType.IMAGE:
            result = _get_naxis1(header) / 2.0
        else:
            obs_type = header.get('OBSTYPE')
            if obs_type == 'dark':
                result = 1.0
    return result


def get_position_function_coord2_pix(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            result = _get_naxis2(header) / 2.0
    else:
        if data_product_type == DataProductType.IMAGE:
            result = _get_naxis2(header) / 2.0
        else:
            obs_type = header.get('OBSTYPE')
            if obs_type == 'dark':
                result = 1.0
    return result


def get_position_function_coord1_val(header):
    ra, ignore_dec = _get_position(header)
    return ra


def get_position_function_coord2_val(header):
    ignore_ra, dec = _get_position(header)
    return dec


def _pattern(header, science_spectrum, science_image):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = science_spectrum()
        else:
            result = science_image()
    return result


def get_position_function_cd11(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = -0.001388
        else:
            platescale = mc.to_float(header.get('PLTSCALE'))
            pixsize = mc.to_float(header.get('PIXSIZE'))
            xbin = mc.to_float(header.get('XBIN'))
            if (platescale is not None and
                    pixsize is not None and
                    xbin is not None):
                logging.error('am i here?')
                result = platescale * pixsize * xbin / 3600000.0
    else:
        obs_type = header.get('OBSTYPE')
        if data_product_type == DataProductType.IMAGE and obs_type == 'dark':
            platescale = mc.to_float(header.get('PLTSCALE', 0.0))
            pixsize = mc.to_float(header.get('PIXSIZE', 0.0))
            xbin = mc.to_float(header.get('XBIN', 0.0))
            logging.error('must be here then')
            result = platescale * pixsize * xbin / 3600000.0
        else:
            if obs_type == 'dark':
                result = -0.001388
    return result


def get_position_function_cd12(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 0.0
            # naxis1 = _get_naxis1(header)
            # if naxis1 is None:
            #     result = 0.0
            # else:
            #     result = mc.to_float(naxis1) / 2.0
        else:
            result = 0.0
    else:
        obs_type = header.get('OBSTYPE')
        if obs_type == 'dark':
            result = 0.0
    return result


def get_position_function_cd21(header):
    logging.error('called')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 0.0
            # naxis2 = _get_naxis2(header)
            # if naxis2 is None:
            #     result = 0.0
            # else:
            #     result = mc.to_float(naxis2) / 2.0
        else:
            result = 0.0
    else:
        obs_type = header.get('OBSTYPE')
        if obs_type == 'dark':
            result = 0.0
    return result


def get_position_function_cd22(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = 0.001388
        else:
            platescale = mc.to_float(header.get('PLTSCALE'))
            pixsize = mc.to_float(header.get('PIXSIZE'))
            xbin = mc.to_float(header.get('XBIN'))
            if (platescale is not None and
                    pixsize is not None and
                    xbin is not None):
                result = platescale * pixsize * xbin / 3600000.0
    else:
        obs_type = header.get('OBSTYPE')
        if data_product_type == DataProductType.IMAGE and obs_type == 'dark':
            platescale = mc.to_float(header.get('PLTSCALE', 0.0))
            pixsize = mc.to_float(header.get('PIXSIZE', 0.0))
            xbin = mc.to_float(header.get('XBIN', 0.0))
            logging.error('must also be here then')
            result = platescale * pixsize * xbin / 3600000.0
        else:
            if obs_type == 'dark':
                result = 0.001388
    return result


def get_position_function_dimension_naxis1(header):
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        result = 1
    else:
        result = header.get('NAXIS1')
    return result


def get_position_function_dimension_naxis2(header):
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        result = 1
    else:
        result = header.get('NAXIS2')
    return result


def get_target_type(header):
    target_type = TargetType.FIELD
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        target_type = TargetType.OBJECT
    return target_type


def get_telescope_name(header):
    telescope = _get_telescope(header)
    observatory = header.get('OBSERVAT')
    return f'{observatory} {telescope}'


def get_time_axis_delta(header):
    exptime = get_time_exposure(header)
    return exptime / (24.0 * 3600.0)


def get_time_axis_val(header):
    return ac.get_datetime(header.get('DATE-OBS'))


def get_time_exposure(header):
    exptime = header.get('EXPTIME')
    ncombine = mc.to_float(header.get('NCOMBINE'))
    if ncombine is not None:
        exptime = exptime * ncombine
    return exptime


def get_time_resolution(header):
    exptime = header.get('EXPTIME')
    ncombine = mc.to_float(header.get('NCOMBINE'))
    if ncombine is None:
        ncombine = 1
    else:
        exptime = exptime * ncombine
    return exptime / ncombine


def _get_dispaxis(header):
    dispaxis = None
    if get_data_product_type(header) == DataProductType.SPECTRUM:
        dispaxis = header.get('DISPAXIS')
        logging.error(f'dispaxis from header is {dispaxis}')
        if dispaxis is None:
            telescope = get_telescope_name(header)
            if telescope == 'DAO 1.2-m':
                logging.error('am I here? for dispaxis?')
                dispaxis = 2
            else:
                dispaxis = 1
    if dispaxis is None:
        raise mc.CadcException('Could not determine dispaxis for TODO.')
    return dispaxis


def _get_geo(header):
    telescope = _get_telescope(header)
    if telescope == '1.2-m':
        return ac.get_location(48.52092, -123.42006, 225.0)
    elif telescope == '1.8-m':
        return ac.get_location(48.51967, -123.41833, 232.0)
    else:
        raise mc.CadcException(f'Unexpected telescope value of {telescope} for '
                               f'{header.get("DAOPRGID")}')


def _get_naxis1(header):
    return header.get('NAXIS1')


def _get_naxis2(header):
    return header.get('NAXIS2')


def _get_obs_type(header):
    return header.get('OBSTYPE')


def _get_position(header):
    if header.get('EQUINOX') is None:
        return None, None
    else:
        ra = header.get('RA', 0)
        dec = header.get('DEC', 0)
        equinox = f'J{header.get("EQUINOX")}'
        fk5 = FK5(equinox=equinox)
        coord = SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg), frame=fk5)
        j2000 = FK5(equinox='J2000')
        result = coord.transform_to(j2000)
        return result.ra.degree, result.dec.degree


def _get_telescope(header):
    return header.get('TELESCOP')


def _get_wavelength(header):
    return mc.to_float(header.get('WAVELENG'))


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)
    bp.configure_energy_axis(4)
    bp.configure_observable_axis(5)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.clear('Observation.metaRelease')
    bp.add_fits_attribute('Observation.metaRelease',  'RELEASE')
    bp.add_fits_attribute('Observation.metaRelease',  'DATE-OBS')

    bp.set('Observation.telescope.name',  'get_telescope_name(header)')
    bp.set('Observation.telescope.geoLocationX',  'get_geo_x(header)')
    bp.set('Observation.telescope.geoLocationY',  'get_geo_y(header)')
    bp.set('Observation.telescope.geoLocationZ',  'get_geo_z(header)')

    bp.set('Observation.target.type',  'get_target_type(header)')
    bp.clear('Observation.target.moving')
    bp.set_default('Observation.target.moving',  'false')
    bp.clear('Observation.target.standard')
    bp.set_default('Observation.target.standard',  'false')

    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'DAOPRGID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PINAME')

    bp.clear('Observation.environment.humidity')
    bp.add_fits_attribute('Observation.environment.humidity', 'REL_HUMI')
    bp.clear('Observation.environment.photometric')
    bp.set_default('Observation.environment.photometric', 'false')

    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(parameters)')
    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease',  'DATE-OBS')

    bp.set('Plane.provenance.project', 'DAO Science Archive')
    bp.clear('Plane.provenance.name')
    bp.add_fits_attribute('Plane.provenance.name', 'PROCNAME')
    bp.set_default('Plane.provenance.name', 'DAO unprocessed data')
    # bp.set('Plane.provenance.name', 'DAO unprocessed data')
    bp.set('Plane.provenance.producer', 'NRC Herzberg')
    bp.set('Plane.provenance.reference',
           'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/dao/')
    bp.clear('Plane.provenance.version')
    bp.add_fits_attribute('Plane.provenance.version', 'PROCVERS')
    bp.set_default('Plane.provenance.version', None)

    bp.set('Artifact.productType', 'get_artifact_product_type(header)')

    bp.set('Chunk.time.axis.axis.ctype', 'TIME')
    bp.set('Chunk.time.axis.axis.cunit', 'd')
    bp.set('Chunk.time.axis.function.naxis', '1')
    bp.set('Chunk.time.axis.function.delta', 'get_time_axis_delta(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
    bp.set(
        'Chunk.time.axis.function.refCoord.val', 'get_time_axis_val(header)')
    bp.set('Chunk.time.exposure', 'get_time_exposure(header)')
    bp.set('Chunk.time.resolution', 'get_time_resolution(header)')

    bp.set('Chunk.observable.axis.axis.ctype', 'FLUX')
    bp.set('Chunk.observable.axis.axis.cunit', 'COUNTS')
    bp.set('Chunk.observable.axis.function.refCoord.pix', 1)

    bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    bp.set('Chunk.energy.axis.axis.cunit', 'Angstrom')
    bp.set('Chunk.energy.axis.function.delta',
           'get_energy_function_delta(header)')
    bp.set('Chunk.energy.axis.function.naxis',
           'get_energy_function_naxis(parameters)')
    bp.set('Chunk.energy.axis.function.refCoord.pix',
           'get_energy_function_pix(header)')
    bp.clear('Chunk.energy.axis.function.refCoord.val')
    bp.add_fits_attribute('Chunk.energy.axis.function.refCoord.val',
                          'WAVELENG')
    bp.set('Chunk.energy.specsys', 'TOPOCENT')
    bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
    bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
    bp.set('Chunk.energy.resolvingPower', 'get_energy_resolving_power(header)')
    bp.clear('Chunk.energy.bandpassName')
    bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')

    bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
    bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
    bp.set('Chunk.position.axis.axis1.cunit', 'deg')
    bp.set('Chunk.position.axis.axis2.cunit', 'deg')
    bp.set('Chunk.position.axis.function.dimension.naxis1',
           'get_position_function_dimension_naxis1(header)')
    bp.set('Chunk.position.axis.function.dimension.naxis2',
           'get_position_function_dimension_naxis2(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord1.pix',
           'get_position_function_coord1_pix(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord1.val',
           'get_position_function_coord1_val(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix',
           'get_position_function_coord2_pix(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord2.val',
           'get_position_function_coord2_val(header)')
    bp.set('Chunk.position.axis.function.cd11',
           'get_position_function_cd11(header)')
    bp.set('Chunk.position.axis.function.cd22',
           'get_position_function_cd22(header)')
    bp.set('Chunk.position.axis.function.cd12',
           'get_position_function_cd12(header)')
    bp.set('Chunk.position.axis.function.cd21',
           'get_position_function_cd21(header)')

    # composites
    if re.match(r'ad:DAO/dao_[c]\d{3}_\d{4}_\d{6}_[aBF].\w', uri):
        bp.set('CompositeObservation.members', 'get_members(header)')
    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes (an n:n
    relationship between TDM attributes and CAOM attributes). Must have this
    signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = None
    if 'headers' in kwargs:
        headers = kwargs['headers']
    fqn = None
    if 'fqn' in kwargs:
        fqn = kwargs['fqn']

    # correct the *_axis values
    for plane in observation.planes.values():
        if plane.data_product_type == DataProductType.SPECTRUM:
            for artifact in plane.artifacts.values():

                # provenance inputs
                _update_provenance(plane, artifact, headers)
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        chunk.observable_axis = 2
                        chunk.time_axis = 5
                        chunk.energy_axis = 1
                        if DAOName.is_unprocessed_reticon(artifact.uri):
                            cc.reset_energy(chunk)
                        if artifact.product_type == ProductType.SCIENCE:
                            logging.error('science')
                            chunk.position_axis_1 = 3
                            chunk.position_axis_2 = 4
                        else:
                            logging.error(f'not science {artifact.product_type}'
                                          f' type {observation.type}')
                            if observation.type == 'dark':
                                chunk.position_axis_1 = 3
                                chunk.position_axis_2 = 4
                            else:
                                cc.reset_position(chunk)
                            # no energy for calibration?
                            if observation.type not in ['flat', 'comparison', 'dark']:
                                cc.reset_energy(chunk)
        else:  # DataProductType.IMAGE
            logging.error('image')
            for artifact in plane.artifacts.values():
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        # no observable axis when image
                        cc.reset_observable(chunk)
                        if artifact.product_type == ProductType.CALIBRATION:
                            if observation.type != 'dark':
                                cc.reset_position(chunk)
                            if observation.type not in ['flat', 'dark']:
                                cc.reset_energy(chunk)

    logging.debug('Done update.')
    return observation


def _update_provenance(plane, artifact, headers):
    pass
    # if DAOName.is_processed(artifact.uri):
    #     cc.update_plane_provenance(plane, headers,)
# def update_chunk_position(chunk, headers):
#     naxis1 = headers[0].get('NAXIS1')
#     naxis2 = headers[0].get('NAXIS2')
#     ra = headers[0].get('RA')
#     dec = headers[0].get('DEC')
#     equinox = headers[0].get('EQUINOX')
#     platescale = mc.to_float(headers[0].get('PLTSCALE'))
#     pixsize = mc.to_float(headers[0].get('PIXSIZE'))
#     xbin = mc.to_float(headers[0].get('XBIN'))
#     cd11 = platescale * pixsize * xbin / 3600000.0
#     cd22 = cd11
#     cd12 = mc.to_float(naxis1) / 2.0
#     cd21 = cd12
#     chunk.position_axis_1 = 1
#     chunk.position_axis_2 = 2
#     from caom2 import SpatialWCS
#     chunk.position = SpatialWCS()


def _build_blueprints(uri):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uri The artifact URI for the file to be processed."""
    module = importlib.import_module(__name__)
    blueprint = ObsBlueprint(module=module)
    accumulate_bp(blueprint, uri)
    blueprints = {uri: blueprint}
    return blueprints


def _get_uri(args):
    if args.lineage:
        ignore, uri = mc.decompose_lineage(args.lineage[0])
    elif args.local:
        obs_id = mc.StorageName.remove_extensions(
            os.path.basename(args.local[0]))
        uri = DAOName(obs_id=obs_id).file_uri
    elif args.observation:
        uri = DAOName(obs_id=args.observation[1]).file_uri
    else:
        raise mc.CadcException(f'Could not define uri from these args {args}')
    return uri


def to_caom2():
    args = get_gen_proc_arg_parser().parse_args()
    uri = _get_uri(args)
    blueprints = _build_blueprints(uri)
    return gen_proc(args, blueprints)


def dao_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = to_caom2()
        logging.debug(f'Done {APPLICATION} processing.')
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)
