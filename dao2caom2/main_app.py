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
from caom2 import ObservationIntentType
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc


__all__ = ['main_app', 'update', 'DAOName', 'COLLECTION', 'APPLICATION']


APPLICATION = 'dao2caom2'
COLLECTION = 'DAO'


class DAOName(ec.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support compressed raw files and uncompressed product files in storage
    - uncompressed product files have an extension added to the input
      names - lower case if there's a single input, upper case if there's
      multiple inputs, as far as I can tell right now.
    """

    DAO_NAME_PATTERN = '*'

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None):
        self.fname_in_ad = file_name
        super(DAOName, self).__init__(
            obs_id, COLLECTION, DAOName.DAO_NAME_PATTERN, fname_on_disk)

    def is_valid(self):
        return True


def get_artifact_product_type(header):
    obs_type = _get_obs_type(header)
    if obs_type == 'object':
        product_type = ProductType.SCIENCE
    else:
        product_type = ProductType.CALIBRATION
    return product_type


def get_calibration_level(uri):
    from caom2 import CalibrationLevel
    result = CalibrationLevel.RAW_STANDARD
    return result


def get_data_product_type(header):
    obs_mode = header.get('OBSMODE')
    if obs_mode is not None and obs_mode.strip() == 'Imaging':
        data_product_type = DataProductType.IMAGE
    else:
        data_product_type = DataProductType.SPECTRUM
    return data_product_type


def get_energy_function_naxis(header):
    naxis = 1
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        dispaxis = _get_dispaxis(header)
        if dispaxis == 1:
            naxis = header.get('NAXIS1')
        elif dispaxis == 2:
            naxis = header.get('NAXIS2')
        else:
            raise mc.CadcException('Could not find dispaxis for {}'.format('TODO'))
    return naxis


def get_energy_function_delta(header):
    cdelt = 1.0
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        wavelength = _get_wavelength(header)
        if wavelength is None:
            pass
        else:
            logging.error('hows the wavelength?')
            dispersion = header.get('DISPERSI')
            dispaxis = _get_dispaxis(header)
            if dispaxis == 1:
                xbin = mc.to_float(header.get('XBIN'))
            else:
                xbin = mc.to_float(header.get('YBIN'))
            cdelt = dispersion * 15.0 * xbin / 1000.0
            logging.error('cdelt is {} dispersion is {} xbin is {}'.format(cdelt, dispersion, xbin))
    else:
        cdelt = header.get('BANDPASS')
    return cdelt


def get_energy_function_pix(header):
    crpix = 1.0
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.SPECTRUM:
        wavelength = _get_wavelength(header)
        if wavelength is not None:
            datasec = re.sub(r'(\[)(\d+:\d+,\d+:\d+)(\])', r'\g<2>',
                             header.get('DATASEC'))
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
    band_pass = mc.to_float(header.get('BANDPASS'))
    wavelength = _get_wavelength(header)
    data_product_type = get_data_product_type(header)
    if data_product_type == DataProductType.IMAGE:
        resolving_power = wavelength / band_pass
    else:
        # assume 2.5 pixel wide resolution element
        if band_pass is None:
            cdelt = get_energy_function_delta(header)
            resolving_power = wavelength / (2.5 * cdelt)
        else:
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
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            result = _get_naxis1(header) / 2.0
    return result


def get_position_function_coord2_pix(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            result = 1.0
        else:
            result = _get_naxis2(header) / 2.0
    return result


def get_position_function_coord1_val(header):
    ra = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        ra, ignore_dec = _get_position(header)
    return ra


def get_position_function_coord2_val(header):
    dec = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
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
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            # I've set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = -0.001388
        else:
            platescale = mc.to_float(header.get('PLTSCALE'))
            pixsize = mc.to_float(header.get('PIXSIZE'))
            xbin = mc.to_float(header.get('XBIN'))
            result = platescale * pixsize * xbin / 3600000.0
    return result


def get_position_function_cd12(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            naxis1 = _get_naxis1(header)
            result = mc.to_float(naxis1) / 2.0
        else:
            result = 0.0
    return result


def get_position_function_cd21(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            naxis2 = _get_naxis2(header)
            result = mc.to_float(naxis2) / 2.0
        else:
            result = 0.0
    return result


def get_position_function_cd22(header):
    result = None
    artifact_product_type = get_artifact_product_type(header)
    if artifact_product_type is ProductType.SCIENCE:
        data_product_type = get_data_product_type(header)
        if data_product_type == DataProductType.SPECTRUM:
            # I've set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = 0.001388
        else:
            platescale = mc.to_float(header.get('PLTSCALE'))
            pixsize = mc.to_float(header.get('PIXSIZE'))
            xbin = mc.to_float(header.get('XBIN'))
            result = platescale * pixsize * xbin / 3600000.0
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
    return '{} {}'.format(observatory, telescope)


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
        logging.error('dispaxis is {}'.format(dispaxis))
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
        raise mc.CadcException(
            'Unexpected telescope value of {} for {}'.format(
                telescope, header.get('DAOPRGID')))


def _get_naxis1(header):
    return header.get('NAXIS1')


def _get_naxis2(header):
    return header.get('NAXIS2')


def _get_obs_type(header):
    return header.get('OBSTYPE')


def _get_position(header):
    ra = header.get('RA')
    dec = header.get('DEC')
    equinox = 'J{}'.format(header.get('EQUINOX'))
    fk5 = FK5(equinox=equinox)
    coord = SkyCoord(
        '{} {}'.format(ra, dec), unit=(u.hourangle, u.deg), frame=fk5)
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
    bp.set_default('Observation.proposal.id', 'No PRG ID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PINAME')

    bp.clear('Observation.environment.humidity')
    bp.add_fits_attribute('Observation.environment.humidity', 'REL_HUMI')
    bp.clear('Observation.environment.photometric')
    bp.set_default('Observation.environment.photometric', 'false')

    bp.set('Plane.dataProductType', 'get_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')

    bp.set('Plane.provenance.project', 'DAO Science Archive')
    bp.set('Plane.provenance.name', 'DAO unprocessed data')
    bp.set('Plane.provenance.producer', 'NRC Herzberg')
    bp.set('Plane.provenance.reference',
           'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/dao/')

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
           'get_energy_function_naxis(header)')
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
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

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
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        chunk.observable_axis = 2
                        chunk.time_axis = 5
                        chunk.energy_axis = 1
                        if artifact.product_type == ProductType.SCIENCE:
                            chunk.position_axis_1 = 3
                            chunk.position_axis_2 = 4
                        else:
                            chunk.position_axis_1 = None
                            chunk.position_axis_2 = None
                            chunk.position = None
                            # no energy for calibration?
                            if observation.type not in ['flat', 'comparison']:
                                chunk.energy_axis = None
                                chunk.energy = None
        else:  # DataProductType.IMAGE
            for artifact in plane.artifacts.values():
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        # no observable axis when image
                        chunk.observable_axis = None
                        chunk.observable = None
                        if artifact.product_type == ProductType.CALIBRATION:
                            chunk.position = None
                            chunk.position_axis_1 = None
                            chunk.position_axis_2 = None
                            if observation.type != 'flat':
                                chunk.energy_axis = None
                                chunk.energy = None

    logging.debug('Done update.')
    return observation


def update_chunk_position(chunk, headers):
    naxis1 = headers[0].get('NAXIS1')
    naxis2 = headers[0].get('NAXIS2')
    ra = headers[0].get('RA')
    dec = headers[0].get('DEC')
    equinox = headers[0].get('EQUINOX')
    platescale = mc.to_float(headers[0].get('PLTSCALE'))
    pixsize = mc.to_float(headers[0].get('PIXSIZE'))
    xbin = mc.to_float(headers[0].get('XBIN'))
    cd11 = platescale * pixsize * xbin / 3600000.0
    cd22 = cd11
    cd12 = mc.to_float(naxis1) / 2.0
    cd21 = cd12
    chunk.position_axis_1 = 1
    chunk.position_axis_2 = 2
    from caom2 import SpatialWCS
    chunk.position = SpatialWCS()


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
    result = None
    if args.observation:
        result = DAOName(obs_id=args.observation[1]).file_uri
    elif args.local:
        obs_id = ec.StorageName.remove_extensions(
            os.path.basename(args.local[0]))
        result = DAOName(obs_id=obs_id).file_uri
    elif args.lineage:
        result = args.lineage[0].split('/', 1)[1]
    else:
        raise mc.CadcException(
            'Could not define uri from these args {}'.format(args))
    return result


def main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        uri = _get_uri(args)
        blueprints = _build_blueprints(uri)
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
