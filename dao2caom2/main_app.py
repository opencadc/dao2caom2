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
from enum import Enum

from caom2 import Observation, TargetType, DataProductType, ProductType
from caom2 import ObservationIntentType, CalibrationLevel, TypedSet, PlaneURI
from caom2 import ObservationURI
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from dao2caom2 import dao_name as dn


__all__ = ['dao_main_app', 'update', 'APPLICATION', 'to_caom2']


APPLICATION = 'dao2caom2'


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
    if dn.DAOName.is_processed(uri):
        result = CalibrationLevel.CALIBRATED
    return result


def get_data_product_type(header):
    obs_mode = header.get('OBSMODE')
    if obs_mode is not None and obs_mode.strip() == 'Imaging':
        data_product_type = DataProductType.IMAGE
    else:
        data_product_type = DataProductType.SPECTRUM
    return data_product_type


def get_energy_axis_function_naxis(parameters):
    uri = parameters.get('uri')
    header = parameters.get('header')
    naxis = 1
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        if dn.DAOName.is_processed(uri):
            naxis = _get_naxis1(header)
        else:
            dispaxis = _get_dispaxis(header)
            if dispaxis == 1:
                naxis = _get_naxis1(header)
            elif dispaxis == 2:
                naxis = _get_naxis2(header)
            else:
                raise mc.CadcException(
                    f'Could not find dispaxis for \'TODO\'')
    return naxis


def get_energy_axis_function_delta(parameters):
    uri = parameters.get('uri')
    header = parameters.get('header')
    cdelt = 1.0
    execution_path = _get_execution_path(parameters)
    if execution_path is ExecutionPath.SPECT_CALIBRATED:
        cdelt = header.get('CDELT1')
    # elif execution_path is ExecutionPath.SPECT_RAW:
    #     cdelt = header.get('DELTA_WL')
    else:
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


def get_energy_axis_function_refcoord_pix(parameters):
    header = parameters.get('header')
    crpix = 1.0
    execution_path = _get_execution_path(parameters)
    if execution_path is ExecutionPath.SPECT_CALIBRATED:
        crpix = header.get('CRPIX1')
    else:
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
    logging.error(f'get_energy_axis_function_refcoord_pix {execution_path} {crpix}')
    return crpix


def get_energy_axis_function_refcoord_val(parameters):
    header = parameters.get('header')
    execution_path = _get_execution_path(parameters)
    result = None
    if execution_path in [ExecutionPath.IMAGING, ExecutionPath.SPECT_RAW]:
        result = header.get('WAVELENG')
    elif execution_path is ExecutionPath.SPECT_CALIBRATED:
        result = header.get('CRVAL1')
    logging.error(f'get_energy_axis_function_refcoord val {result}')
    return result


def get_energy_resolving_power(parameters):
    resolving_power = None
    execution_path = _get_execution_path(parameters)
    numerator = get_energy_axis_function_refcoord_val(parameters)
    denominator = get_energy_axis_function_delta(parameters)
    if numerator is not None and denominator is not None:
        if execution_path is ExecutionPath.IMAGING:
            resolving_power = numerator / denominator
        elif execution_path in [ExecutionPath.SPECT_RAW,
                                ExecutionPath.SPECT_CALIBRATED]:
            resolving_power = numerator / (2.5 * denominator)
    return resolving_power


class ExecutionPath(Enum):
    IMAGING = 1,
    SPECT_CALIBRATED = 2,
    SPECT_RAW = 3


def _get_execution_path(parameters):
    uri = parameters.get('uri')
    header = parameters.get('header')
    obs_mode = _get_obs_mode(header)
    obs_type = _get_obs_type(header)
    if obs_mode == 'imaging':
        result = ExecutionPath.IMAGING
    else:
        result = ExecutionPath.SPECT_RAW
        if (dn.DAOName.is_processed(uri) and obs_type in ['object',
                                                          'comparison']):
            result = ExecutionPath.SPECT_CALIBRATED
    return result


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


def get_position_function_coord1_pix(parameters):
    header = parameters.get('header')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            calibration_level = get_calibration_level(parameters)
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CRPIX1')
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


def get_position_function_coord2_pix(parameters):
    header = parameters.get('header')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            calibration_level = get_calibration_level(parameters)
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CRPIX2')
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


def get_position_function_cd11(parameters):
    header = parameters.get('header')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    calibration_level = get_calibration_level(parameters)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = -0.001388
        else:
            if calibration_level is CalibrationLevel.RAW_STANDARD:
                platescale = mc.to_float(header.get('PLTSCALE'))
                pixsize = mc.to_float(header.get('PIXSIZE'))
                xbin = mc.to_float(header.get('XBIN'))
                if (platescale is not None and
                        pixsize is not None and
                        xbin is not None):
                    logging.error('am i here?')
                    result = platescale * pixsize * xbin / 3600000.0
            else:
                result = header.get('CD1_1')
    else:
        obs_type = header.get('OBSTYPE')
        if data_product_type == DataProductType.IMAGE:
            if obs_type == 'dark':
                platescale = mc.to_float(header.get('PLTSCALE', 0.0))
                pixsize = mc.to_float(header.get('PIXSIZE', 0.0))
                xbin = mc.to_float(header.get('XBIN', 0.0))
                logging.error('must be here then')
                result = platescale * pixsize * xbin / 3600000.0
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CD1_1')
        else:
            if obs_type == 'dark':
                result = -0.001388
    return result


def get_position_function_cd12(parameters):
    header = parameters.get('header')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    calibration_level = get_calibration_level(parameters)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 0.0
        else:
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CD1_2')
            else:
                result = 0.0
    else:
        obs_type = header.get('OBSTYPE')
        if obs_type == 'dark':
            result = 0.0
    return result


def get_position_function_cd21(parameters):
    header = parameters.get('header')
    result = None
    calibration_level = get_calibration_level(parameters)
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 0.0
        else:
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CD2_1')
            else:
                result = 0.0
    else:
        obs_type = header.get('OBSTYPE')
        if obs_type == 'dark':
            result = 0.0
    return result


def get_position_function_cd22(parameters):
    header = parameters.get('header')
    result = None
    artifact_product_type = get_artifact_product_type(header)
    data_product_type = get_data_product_type(header)
    calibration_level = get_calibration_level(parameters)
    if artifact_product_type is ProductType.SCIENCE:
        if data_product_type == DataProductType.SPECTRUM:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = 0.001388
        else:
            if calibration_level is CalibrationLevel.CALIBRATED:
                result = header.get('CD2_2')
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
    result = f'{observatory} {telescope}'
    if telescope is None or observatory is None:
        result = None
    return result


def get_time_axis_delta(header):
    exptime = get_time_exposure(header)
    return exptime / (24.0 * 3600.0)


def get_time_axis_val(header):
    return ac.get_datetime(header.get('DATE-OBS'))


def get_time_exposure(header):
    exptime = header.get('EXPTIME')
    ncombine = mc.to_float(header.get('NCOMBINE'))
    if ncombine is not None:
        # DB - approximation of exposure time for products (assume identical
        # EXPTIME)
        exptime *= ncombine
    logging.error(f'exptime {exptime} ncombine {ncombine}')
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
        if dispaxis is None:
            telescope = get_telescope_name(header)
            if telescope == 'DAO 1.2-m':
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


def _get_obs_mode(header):
    """
    1. OBSMODE keyword => 'imaging' vs 'spectroscopy'
    2. File name pattern

    All energy is CTYPE = WAVE, CUNIT = A(ngstrom)
    File Name Pattern
    all imaging:                       CRVAL    'WAVELENG'
                                       CDELT    'BANDPASS'
                                       CRPIX    1
                                       RP       CRVAL / CDELT

    processed photographic plate spectrum
    dao_[p]\d{3}_\d{6}(u|v|y|r|i|):
                                       CRVAL   'CRVAL1'
                                       CDELT   'CDELT1'
                                       CRPIX   'CRPIX1'
                                       RP      CRVAL/(2.5*CDELT)

    unprocessed DAO RETICON spectrum
    dao_[r]\d{3}_\d{4}_\d{6}:
                                       CRVAL   ''
                                       CDELT   ''
                                       CRPIX   ''
                                       RP      ''

    unprocessed DAO CCD spectrum
    dao_[c]\d{3}_\d{4}_\d{6}:
                                       CRVAL   'WAVELENG'
                                       CDELT   'DELTA_WL'
                                       CRPIX   'REFPIXEL'
                                       RP      CRVAL/(2.5*CDELT)

                                       If no 'WAVELENG' USE 'DATASEC'

                                       CDELT   'DISPERSI' * 15.0 * xbin/1000.0
                                       CRPIX   'DATASEC' + math

    processed DAO spectrum
    dao_[cr]\d{3}_\d{4}_\d{6}_[evBF], obs.type in ['object', 'comparison']:
                                       CRVAL   'CRVAL1'
                                       CDELT   'CDELT1'
                                       CRPIX   'CRPIX1'
                                       RP      CRVAL/(2.5*CDELT)

    :param header:
    :return:
    """
    obs_mode = header.get('OBSMODE')
    result = 'imaging'
    if '-slit' in obs_mode:
        result = 'spectroscopy'
    return result


def _get_obs_type(header):
    return header.get('OBSTYPE')


def _get_position(header):
    obs_type = _get_obs_type(header)
    ra_deg = None
    dec_deg = None
    if obs_type in ['comparison', 'dark', 'object']:
        if header.get('EQUINOX') is not None:
            # DB - 11-09-19 - precession with astropy
            ra = header.get('RA', 0)
            dec = header.get('DEC', 0)
            equinox = f'J{header.get("EQUINOX")}'
            fk5 = FK5(equinox=equinox)
            coord = SkyCoord(f'{ra} {dec}', unit=(u.hourangle, u.deg), frame=fk5)
            j2000 = FK5(equinox='J2000')
            result = coord.transform_to(j2000)
            ra_deg = result.ra.degree
            dec_deg = result.dec.degree
    return ra_deg, dec_deg


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
    # from dao2caom2.config
    bp.add_fits_attribute('Observation.metaRelease',  'DATE-OBS')

    bp.clear('Observation.algorithm.name')

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
    # from dao2caom2.config
    bp.add_fits_attribute('Plane.metaRelease',  'DATE-OBS')

    bp.clear('Plane.provenance.lastExecuted')
    bp.add_fits_attribute('Plane.provenance.lastExecuted', 'IRAF-TLM')
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
           'get_energy_axis_function_delta(parameters)')
    bp.set('Chunk.energy.axis.function.naxis',
           'get_energy_axis_function_naxis(parameters)')
    bp.set('Chunk.energy.axis.function.refCoord.pix',
           'get_energy_axis_function_refcoord_pix(parameters)')
    bp.set('Chunk.energy.axis.function.refCoord.val',
           'get_energy_axis_function_refcoord_val(parameters)')
    bp.set('Chunk.energy.specsys', 'TOPOCENT')
    bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
    bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
    bp.set('Chunk.energy.resolvingPower',
           'get_energy_resolving_power(parameters)')
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
           'get_position_function_coord1_pix(parameters)')
    bp.set('Chunk.position.axis.function.refCoord.coord1.val',
           'get_position_function_coord1_val(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix',
           'get_position_function_coord2_pix(parameters)')
    bp.set('Chunk.position.axis.function.refCoord.coord2.val',
           'get_position_function_coord2_val(header)')
    bp.set('Chunk.position.axis.function.cd11',
           'get_position_function_cd11(parameters)')
    bp.set('Chunk.position.axis.function.cd22',
           'get_position_function_cd22(parameters)')
    bp.set('Chunk.position.axis.function.cd12',
           'get_position_function_cd12(parameters)')
    bp.set('Chunk.position.axis.function.cd21',
           'get_position_function_cd21(parameters)')

    # derived observations
    if dn.DAOName.is_derived(uri):
        bp.set('CompositeObservation.members', 'get_members(header)')
        bp.add_fits_attribute('Observation.algorithm.name', 'PROCNAME')
    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes (an n:n
    relationship between TDM attributes and CAOM attributes). Must have this
    signature for import_module loading and execution.

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = kwargs.get('headers')
    fqn = kwargs.get('fqn')
    uri = kwargs.get('uri')
    dao_name = None
    if uri is not None:
        dao_name = dn.DAOName(artifact_uri=uri)
    if fqn is not None:
        dao_name = dn.DAOName(file_name=os.path.basename(fqn))

    if dao_name is None:
        raise mc.CadcException(f'Need one of fqn or uri defined for '
                               f'{observation.observation_id}')
    # correct the *_axis values
    for plane in observation.planes.values():
        for artifact in plane.artifacts.values():
            if artifact.uri.replace('.gz', '') != dao_name.file_uri.replace('.gz', ''):
                continue
            logging.error(f'update artifact {artifact.uri} file {dao_name.file_uri}')

            for part in artifact.parts.values():
                for chunk in part.chunks:
                    time_delta = get_time_axis_delta(headers[0])
                    cc.undo_astropy_cdfix_call(chunk, time_delta)

                    if plane.data_product_type == DataProductType.SPECTRUM:
                        chunk.observable_axis = 2
                        chunk.time_axis = 5
                        chunk.energy_axis = 1
                        if (dn.DAOName.is_unprocessed_reticon(artifact.uri) or
                            dn.DAOName.is_derived(artifact.uri) and
                                observation.type == 'flat'):
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
                            if (observation.type not in
                                    ['flat', 'comparison', 'dark']):
                                cc.reset_energy(chunk)
                    else:  # DataProductType.IMAGE
                        if dn.DAOName.is_processed_image(artifact.uri):
                            plane.provenance.producer = 'Spaceguard_C'
                        # no observable axis when image
                        cc.reset_observable(chunk)
                        if artifact.product_type == ProductType.CALIBRATION:
                            if observation.type != 'dark':
                                cc.reset_position(chunk)
                            if observation.type not in ['flat', 'dark']:
                                cc.reset_energy(chunk)

        if plane.product_id != dao_name.product_id:
            continue

        if observation.type == 'flat' and cc.is_composite(headers, 'FLAT_'):
            cc.update_plane_provenance(plane, headers, 'FLAT_', dn.COLLECTION,
                                       _repair_provenance_value,
                                       observation.observation_id)
        elif observation.type == 'bias' and cc.is_composite(headers, 'ZERO_'):
            cc.update_plane_provenance(plane, headers, 'ZERO_', dn.COLLECTION,
                                       _repair_provenance_value,
                                       observation.observation_id)

        if dn.DAOName.is_processed(dao_name.file_uri):
            _update_plane_provenance(observation, plane, headers)

        if (cc.is_composite(headers, 'FLAT_') or
                cc.is_composite(headers, 'ZERO_')):
            cc.update_observation_members(observation)

    logging.debug('Done update.')
    return observation


def _repair_provenance_value(value, obs_id):
    logging.debug(f'Begin _repair_provenance_value for {obs_id}')
    # values look like:
    # FLAT_1  = 'dao_c122_2007_000916.fits'
    # FLAT_2  = 'dao_c122_2007_000917.fits'
    # FLAT_3  = 'dao_c122_2007_000918.fits'
    # FLAT_4  = 'dao_c122_2007_000919.fits'
    #
    # OR
    #
    # ZERO_18 = 'dao_c122_2016_012728.fits'
    # ZERO_19 = 'dao_c122_2016_012729.fits'
    # ZERO_20 = 'dao_c122_2016_012730.fits'
    dao_name = dn.DAOName(file_name=value)
    prov_prod_id = dao_name.product_id
    prov_obs_id = dao_name.obs_id
    logging.debug(f'End _repair_provenance_value')
    return prov_obs_id, prov_prod_id


def _update_plane_provenance(observation, plane, headers):
    logging.debug(f'Begin _update_plane_provenance for {plane.product_id} with'
                  f'observation type: {observation.type}.')
    if observation.type in ['object', 'flat', 'comparison']:
        f_name = headers[0].get('BIAS')
        if f_name is not None:
            bias_name = dn.DAOName(file_name=f_name)
            plane_uri = _make_uris(bias_name.obs_id, bias_name.product_id)
            plane.provenance.inputs.add(plane_uri)
    if observation.type in ['object', 'comparison']:
        f_name = headers[0].get('FLAT')
        if f_name is not None:
            flat_name = dn.DAOName(file_name=f_name)
            plane_uri = _make_uris(flat_name.obs_id, flat_name.product_id)
            plane.provenance.inputs.add(plane_uri)
        # referral to raw plane
        plane_uri = _make_uris(observation.observation_id,
                               observation.observation_id)
        plane.provenance.inputs.add(plane_uri)
    if observation.type == 'object':
        f_name = headers[0].get('DCLOG1')
        if f_name is not None:
            ref_spec1_name = dn.DAOName(file_name=f_name.split()[2])
            plane_uri = _make_uris(ref_spec1_name.obs_id,
                                   ref_spec1_name.product_id)
            plane.provenance.inputs.add(plane_uri)
        if headers[0].get('DCLOG2') is not None:
            ref_spec1_name = dn.DAOName(
                file_name=headers[0].get('DCLOG2').split()[2])
            plane_uri = _make_uris(ref_spec1_name.obs_id,
                                   ref_spec1_name.product_id)
            plane.provenance.inputs.add(plane_uri)
    logging.debug(f'End _update_plane_provenance.')


def _make_uris(obs_id, product_id):
    obs_member_uri = ObservationURI(mc.CaomName.make_obs_uri_from_obs_id(
        dn.COLLECTION, obs_id))
    plane_uri = PlaneURI.get_plane_uri(obs_member_uri, product_id)
    return plane_uri


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


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DAO blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The artifact URI for the file to be processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.lineage:
        for ii in args.lineage:
            ignore, uri = mc.decompose_lineage(ii)
            result.append(uri)
    elif args.local:
        for ii in args.local:
            obs_id = mc.StorageName.remove_extensions(os.path.basename(ii))
            uri = dn.DAOName(obs_id=obs_id).file_uri
            result.append(uri)
    elif args.observation:
        uri = dn.DAOName(obs_id=args.observation[1]).file_uri
        result.append(uri)
    else:
        raise mc.CadcException(f'Could not define uri from these args {args}')
    return result


def to_caom2():
    args = get_gen_proc_arg_parser().parse_args()
    uris = _get_uris(args)
    blueprints = _build_blueprints(uris)
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
