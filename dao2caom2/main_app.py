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

"""
DB - github comments - 30-04-20

In the original dao2caom2 the FITS header OBSMODE value is checked first to
see if the file is an 'Imaging' or 'Spectroscopy' observation.

Use the INSTRUME value and if "Imager" is in the string it's an image. This is
better than depending on file naming patterns, since if I ever do a special
preview for spectropolarimeter observations it will capture that info as well.
(OBSMODE = single-slit spectroscopy for regular spectroscopy and
spectropolarimetry for historical reasons).

"""
import importlib
import logging
import os
import sys
import traceback

import astropy.io.fits
from astropy.coordinates import SkyCoord, FK5
import astropy.units as u

from dataclasses import dataclass

from cadcdata import FileInfo
from caom2 import Observation, DataProductType, ProductType
from caom2 import ObservationIntentType, TypedSet, PlaneURI
from caom2 import ObservationURI
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
# from caom2utils import gen_proc_with_headers
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from dao2caom2 import dao_name as dn
from dao2caom2 import telescopes


# __all__ = ['dao_main_app', 'update', 'APPLICATION', 'to_caom2']
__all__ = ['dao_main_app', 'APPLICATION']


APPLICATION = 'dao2caom2'


# def get_artifact_product_type(parameters):
#     header = parameters.get('header')
#     uri = parameters.get('uri')
#     return telescopes.get_current(uri).get_artifact_product_type(header)
#
#
# def get_energy_axis_function_naxis(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_energy_axis_function_naxis(header)
#
#
# def get_energy_axis_function_delta(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_energy_axis_function_delta(
#         header
#     )
#
#
# def get_energy_axis_function_refcoord_pix(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_energy_axis_function_refcoord_pix(
#         header
#     )
#
#
# def get_energy_axis_function_refcoord_val(parameters):
#     header = parameters.get('header')
#     uri = parameters.get('uri')
#     return telescopes.get_current(uri).get_energy_axis_function_refcoord_val(
#         header
#     )
#
#
# def get_energy_resolving_power(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_energy_resolving_power(header)
#
#
# def get_geo_x(uri):
#     x, ignore_y, ignore_z = _get_geo(uri)
#     return x
#
#
# def get_geo_y(uri):
#     ignore_x, y, ignore_z = _get_geo(uri)
#     return y
#
#
# def get_geo_z(uri):
#     ignore_x, ignore_y, z = _get_geo(uri)
#     return z
#
#
# def get_members(header):
#     # this function exists so fits2caom2 creates the correct Observation type
#     pass
#
#
# def get_obs_intent(header):
#     obs_type = _get_obs_type(header)
#     if obs_type == 'object':
#         intent = ObservationIntentType.SCIENCE
#     else:
#         intent = ObservationIntentType.CALIBRATION
#     return intent
#
#
# def get_position_function_coord1_pix(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_coord1_pix(
#         header
#     )
#
#
# def get_position_function_coord2_pix(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_coord2_pix(
#         header
#     )
#
#
# def get_position_function_coord1_val(header):
#     ra, ignore_dec = _get_position(header)
#     return ra
#
#
# def get_position_function_coord2_val(header):
#     ignore_ra, dec = _get_position(header)
#     return dec
#
#
# def get_position_function_cd11(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_cd11(header)
#
#
# def get_position_function_cd12(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_cd12(header)
#
#
# def get_position_function_cd21(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_cd21(header)
#
#
# def get_position_function_cd22(parameters):
#     uri = parameters.get('uri')
#     header = parameters.get('header')
#     return telescopes.get_current(uri).get_position_function_cd22(header)
#
#
# def get_skycam_release_date(header):
#     return ac.get_datetime(header.get('CLOCKVAL'))
#
#
# def get_telescope_name(uri):
#     return telescopes.get_current(uri).get_telescope_name()
#
#
# def get_time_axis_delta(header):
#     exptime = get_time_exposure(header)
#     return exptime / (24.0 * 3600.0)
#
#
# def get_time_axis_val(params):
#     header = params.get('header')
#     uri = params.get('uri')
#     return telescopes.get_current(uri).get_time_axis_val(header)
#
#
# def get_time_exposure(header):
#     exptime = mc.to_float(header.get('EXPTIME'))
#     ncombine = mc.to_float(header.get('NCOMBINE'))
#     if ncombine is not None:
#         # DB - approximation of exposure time for products (assume identical
#         # EXPTIME)
#         exptime *= ncombine
#     return exptime
#
#
# def get_time_resolution(header):
#     exptime = mc.to_float(header.get('EXPTIME'))
#     ncombine = mc.to_float(header.get('NCOMBINE'))
#     if ncombine is None:
#         ncombine = 1
#     else:
#         exptime = exptime * ncombine
#     return exptime / ncombine
#
#
# def _get_geo(uri):
#     return telescopes.get_current(uri).get_geo()
#
#
# def _get_obs_type(header):
#     return header.get('OBSTYPE')
#
#
# def _get_position(header):
#     obs_type = _get_obs_type(header)
#     ra_deg = None
#     dec_deg = None
#     if obs_type in ['comparison', 'dark', 'object']:
#         if header.get('EQUINOX') is not None:
#             # DB - 11-09-19 - precession with astropy
#             ra = header.get('RA', 0)
#             dec = header.get('DEC', 0)
#             equinox = f'J{header.get("EQUINOX")}'
#             fk5 = FK5(equinox=equinox)
#             coord = SkyCoord(
#                 f'{ra} {dec}', unit=(u.hourangle, u.deg), frame=fk5
#             )
#             j2000 = FK5(equinox='J2000')
#             result = coord.transform_to(j2000)
#             ra_deg = result.ra.degree
#             dec_deg = result.dec.degree
#     return ra_deg, dec_deg


def accumulate_bp(bp, uri, telescope_data):
    """Configure the telescope-specific ObsBlueprint at the CAOM model
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    # telescopes.factory(uri)
    # # for multi-planed/multi-artifact cases - ensure point to the
    # # correct instance
    # telescopes.get_current(uri).configure_axes(bp)
    # telescopes.get_current(uri).accumulate_bp(bp)
    logging.error(dir(telescope_data))
    telescope_data.configure_axes(bp)
    telescope_data.accumulate_bp(bp)

    meta_producer = mc.get_version(APPLICATION)
    bp.set('Observation.metaProducer', meta_producer)
    bp.set('Plane.metaProducer', meta_producer)
    bp.set('Artifact.metaProducer', meta_producer)
    bp.set('Chunk.metaProducer', meta_producer)

    bp.clear('Observation.algorithm.name')

    bp.set('Observation.telescope.name', 'get_telescope_name(uri)')
    bp.set('Observation.telescope.geoLocationX', 'get_geo_x(uri)')
    bp.set('Observation.telescope.geoLocationY', 'get_geo_y(uri)')
    bp.set('Observation.telescope.geoLocationZ', 'get_geo_z(uri)')

    bp.set('Chunk.time.axis.axis.ctype', 'TIME')
    bp.set('Chunk.time.axis.axis.cunit', 'd')
    bp.set('Chunk.time.axis.function.naxis', '1')
    bp.set('Chunk.time.axis.function.delta', 'get_time_axis_delta(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
    bp.set(
        'Chunk.time.axis.function.refCoord.val', 'get_time_axis_val(params)'
    )

    bp.set('Chunk.observable.axis.axis.ctype', 'FLUX')
    bp.set('Chunk.observable.axis.axis.cunit', 'COUNTS')
    bp.set('Chunk.observable.axis.function.refCoord.pix', 1)

    bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    bp.set('Chunk.energy.axis.axis.cunit', 'Angstrom')
    bp.set('Chunk.energy.specsys', 'TOPOCENT')
    bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
    bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')

    # derived observations
    if dn.DAOName.is_derived(uri):
        bp.set('DerivedObservation.members', 'get_members(header)')
        bp.add_fits_attribute('Observation.algorithm.name', 'PROCNAME')
    logging.debug('Done accumulate_bp.')


# def update(observation, **kwargs):
#     """Called to fill multiple CAOM model elements and/or attributes (an n:n
#     relationship between TDM attributes and CAOM attributes). Must have this
#     signature for import_module loading and execution.
#
#     :param observation A CAOM Observation model instance.
#     :param **kwargs Everything else."""
#     logging.debug('Begin update.')
#     mc.check_param(observation, Observation)
#
#     headers = kwargs.get('headers')
#     fqn = kwargs.get('fqn')
#     uri = kwargs.get('uri')
#     dao_name = None
#     if uri is not None:
#         dao_name = dn.DAOName(uri)
#     if fqn is not None and dao_name is None:
#         dao_name = dn.DAOName(os.path.basename(fqn))
#
#     if dao_name is None:
#         raise mc.CadcException(
#             f'Need one of fqn or uri defined for {observation.observation_id}'
#         )
#     # correct the *_axis values
#     for plane in observation.planes.values():
#         for artifact in plane.artifacts.values():
#             if artifact.uri.replace('.gz', '') != dao_name.file_uri.replace(
#                 '.gz', ''
#             ):
#                 continue
#
#             for part in artifact.parts.values():
#                 for chunk in part.chunks:
#                     time_delta = get_time_axis_delta(headers[0])
#                     cc.undo_astropy_cdfix_call(chunk, time_delta)
#
#                     if dao_name.file_name.startswith('d'):
#                         if plane.data_product_type == DataProductType.SPECTRUM:
#                             if (
#                                 dn.DAOName.is_unprocessed_reticon(artifact.uri)
#                                 or dn.DAOName.is_derived(artifact.uri)
#                                 and observation.type == 'flat'
#                             ):
#                                 cc.reset_energy(chunk)
#                             if (
#                                 artifact.product_type != ProductType.SCIENCE
#                             ):
#                                 if observation.type == 'dark':
#                                     chunk.position_axis_1 = 3
#                                     chunk.position_axis_2 = 4
#                                 else:
#                                     cc.reset_position(chunk)
#                                 # no energy for calibration?
#                                 if observation.type not in [
#                                     'flat',
#                                     'comparison',
#                                     'dark',
#                                 ]:
#                                     cc.reset_energy(chunk)
#                         else:  # DataProductType.IMAGE
#                             if dn.DAOName.override_provenance(artifact.uri):
#                                 plane.provenance.producer = 'Spaceguard_C'
#                             # no observable axis when image
#                             cc.reset_observable(chunk)
#                             if (
#                                     artifact.product_type
#                                     == ProductType.CALIBRATION
#                             ):
#                                 if observation.type != 'dark':
#                                     cc.reset_position(chunk)
#                                 if observation.type not in ['flat', 'dark']:
#                                     cc.reset_energy(chunk)
#                         if (
#                             chunk.energy is not None
#                             and not dn.DAOName.is_processed(artifact.uri)
#                             and headers[0].get('WAVELENG') is None
#                         ):
#                             # DB 16-02-21/04-03-21
#                             #  If WAVELENG isn’t present then all energy
#                             #  metadata should be ignored (spectra and images)
#                             cc.reset_energy(chunk)
#
#                         # WCS axis wisdom from Pat:
#                         #
#                         # In general, assigning axis indices above the value of
#                         # naxis is allowed but more or less pointless. The
#                         # only use case that would justify it is that in a FITS
#                         # file there could be a header with NAXIS=2 and
#                         # WCSAXES=4 which would tell the fits reader to look
#                         # for CTYPE1 through 4 and axes 3 and 4 are metadata.
#                         # Assign those values to Chunk only if you care about
#                         # capturing that the extra wcs metadata was really in
#                         # the fits header and so the order could be preserved;
#                         # in general do not assign the 3 and 4.
#
#                         naxis = headers[0].get('NAXIS')
#                         naxis1 = headers[0].get('NAXIS1')
#                         naxis2 = headers[0].get('NAXIS2')
#                         chunk.naxis = None
#                         chunk.position_axis_1 = None
#                         chunk.position_axis_2 = None
#                         chunk.energy_axis = None
#                         chunk.observable_axis = None
#                         chunk.time_axis = None
#                         if naxis is not None:
#                             if (
#                                 naxis1 is not None
#                                 and naxis2 is not None
#                                 and naxis == 2
#                                 and chunk.position is not None
#                                 and plane.data_product_type
#                                 is DataProductType.IMAGE
#                             ):
#                                 chunk.naxis = 2
#                                 chunk.position_axis_1 = 1
#                                 chunk.position_axis_2 = 2
#                             if (
#                                 naxis1 is not None
#                                 and naxis == 1
#                                 and chunk.energy is not None
#                             ):
#                                 chunk.naxis = 1
#                                 chunk.energy_axis = 1
#                     else:
#                         chunk.energy_axis = None
#
#         if plane.product_id != dao_name.product_id:
#             continue
#
#         # provenance: inputs vs members
#         #
#         # DB - 29-04-20
#         # The inconsistencies are consistent for both telescope for the
#         # derived observations: the processed, co-added flats (files with F
#         # suffix) and co-added biases (files with B suffix).  These should
#         # have the ‘members’ list added to the inputs.  From the definition of
#         # provenance:inputs I’m assuming for the science observations,
#         # processed comparison arcs, and processed flats having a composite
#         # flat and/or bias observation as an ‘input’ is okay rather than
#         # breaking these down into their individual members (since those
#         # derived observations will all be available in the archive with
#         # proper provenance provided).
#
#         if observation.type == 'flat' and cc.is_composite(headers, 'FLAT_'):
#             cc.update_plane_provenance(
#                 plane,
#                 headers,
#                 'FLAT_',
#                 dn.COLLECTION,
#                 _repair_provenance_value,
#                 observation.observation_id,
#             )
#         elif observation.type == 'bias' and cc.is_composite(headers, 'ZERO_'):
#             cc.update_plane_provenance(
#                 plane,
#                 headers,
#                 'ZERO_',
#                 dn.COLLECTION,
#                 _repair_provenance_value,
#                 observation.observation_id,
#             )
#
#         if dn.DAOName.is_processed(dao_name.file_uri):
#             _update_plane_provenance(observation, plane, headers)
#
#         if cc.is_composite(headers, 'FLAT_') or cc.is_composite(
#             headers, 'ZERO_'
#         ):
#             _update_observation_members(observation)
#
#     logging.debug('Done update.')
#     return observation
#
#
# def _repair_provenance_value(value, obs_id):
#     logging.debug(f'Begin _repair_provenance_value for {obs_id}')
#     # values look like:
#     # FLAT_1  = 'dao_c122_2007_000916.fits'
#     # FLAT_2  = 'dao_c122_2007_000917.fits'
#     # FLAT_3  = 'dao_c122_2007_000918.fits'
#     # FLAT_4  = 'dao_c122_2007_000919.fits'
#     #
#     # OR
#     #
#     # ZERO_18 = 'dao_c122_2016_012728.fits'
#     # ZERO_19 = 'dao_c122_2016_012729.fits'
#     # ZERO_20 = 'dao_c122_2016_012730.fits'
#     dao_name = dn.DAOName(value)
#     prov_prod_id = dao_name.product_id
#     prov_obs_id = dao_name.obs_id
#     logging.debug(f'End _repair_provenance_value')
#     return prov_obs_id, prov_prod_id
#
#
# def _update_observation_members(observation):
#     """
#     Must filter results because:
#     DB - 11-06-20
#     For the spectra there is a minor issue with members for master flat,
#     *_F, observations.  The master bias used in the processing, the *_B.fits
#     file, shouldn’t be a member for the master flats.
#
#     The master bias is in the list of inputs though:  Inputs for master flat
#     are the unprocessed flats and the master bias.  The master bias is
#     subtracted pixel-by-pixel from each unprocessed flat as part of the
#     processing before the flats are then co-added.
#
#     The composite/derived master flats (F) and master biases (B) should
#     never be members.  At least for any processing that is currently being
#     done.  For now the only members should those given by the
#     NCOMBINE x ZERO_# or FLAT_# keyword values.
#     """
#
#     def filter_fun(x):
#         result = True
#         if dn.DAOName.is_master_flat(observation.observation_id):
#             if dn.DAOName.is_master_bias(x.get_observation_uri().uri):
#                 result = False
#         return result
#
#     inputs = []
#     members_inputs = TypedSet(
#         ObservationURI,
#     )
#     for plane in observation.planes.values():
#         if (
#             plane.provenance is not None
#             and plane.provenance.inputs is not None
#         ):
#             inputs = filter(filter_fun, plane.provenance.inputs)
#
#     for entry in inputs:
#         members_inputs.add(entry.get_observation_uri())
#         logging.debug(f'Adding Observation URI {entry.get_observation_uri()}')
#     mc.update_typed_set(observation.members, members_inputs)
#
#
# def _update_plane_provenance(observation, plane, headers):
#     logging.debug(
#         f'Begin _update_plane_provenance for {plane.product_id} with'
#         f'observation type: {observation.type}.'
#     )
#     if observation.type in ['object', 'flat', 'comparison']:
#         f_name = headers[0].get('BIAS')
#         if f_name is not None:
#             bias_name = dn.DAOName(f_name)
#             plane_uri = _make_uris(bias_name.obs_id, bias_name.product_id)
#             plane.provenance.inputs.add(plane_uri)
#     if observation.type in ['object', 'comparison']:
#         f_name = headers[0].get('FLAT')
#         if f_name is not None:
#             flat_name = dn.DAOName(f_name)
#             plane_uri = _make_uris(flat_name.obs_id, flat_name.product_id)
#             plane.provenance.inputs.add(plane_uri)
#         # referral to raw plane
#         plane_uri = _make_uris(
#             observation.observation_id, observation.observation_id
#         )
#         plane.provenance.inputs.add(plane_uri)
#     if observation.type == 'object':
#         f_name = headers[0].get('DCLOG1')
#         if f_name is not None:
#             ref_spec1_name = dn.DAOName(f_name.split()[2])
#             plane_uri = _make_uris(
#                 ref_spec1_name.obs_id, ref_spec1_name.product_id
#             )
#             plane.provenance.inputs.add(plane_uri)
#         if headers[0].get('DCLOG2') is not None:
#             ref_spec1_name = dn.DAOName(headers[0].get('DCLOG2').split()[2])
#             plane_uri = _make_uris(
#                 ref_spec1_name.obs_id, ref_spec1_name.product_id
#             )
#             plane.provenance.inputs.add(plane_uri)
#     logging.debug(f'End _update_plane_provenance.')
#
#
# def _make_uris(obs_id, product_id):
#     obs_member_uri = ObservationURI(
#         mc.CaomName.make_obs_uri_from_obs_id(dn.COLLECTION, obs_id)
#     )
#     plane_uri = PlaneURI.get_plane_uri(obs_member_uri, product_id)
#     return plane_uri
#
#
# def _build_blueprints(uris, current_telescopes):
#     """This application relies on the caom2utils fits2caom2 ObsBlueprint
#     definition for mapping FITS file values to CAOM model element
#     attributes. This method builds the DAO blueprint for a single
#     artifact.
#
#     The blueprint handles the mapping of values with cardinality of 1:1
#     between the blueprint entries and the model attributes.
#
#     :param uris The artifact URI for the file to be processed."""
#     module = importlib.import_module(__name__)
#     blueprints = {}
#     for uri in uris:
#         blueprint = ObsBlueprint(module=module)
#         telescope_data = telescopes.factory(uri)
#         current_telescopes[uri] = telescope_data
#         # telescopes.factory(uri)
#         # # for multi-planed/multi-artifact cases - ensure point to the
#         # # correct instance
#         # telescopes.get_current(uri).configure_axes(bp)
#         # telescopes.get_current(uri).accumulate_bp(bp)
#         accumulate_bp(blueprint, uri, telescope_data)
#         blueprints[uri] = blueprint
#     return blueprints
#
#
# def _build_blueprints_with_client(xs, current_telescopes):
#     module = importlib.import_module(__name__)
#     for entry in xs:
#         blueprint = ObsBlueprint(module=module)
#         telescope_data = telescopes.factory_client(entry)
#         current_telescopes[entry.uri] = telescope_data
#         telescope_data.configure_axes(blueprint)
#         accumulate_bp(blueprint, entry.uri, telescope_data)
#         entry.blueprint = blueprint
#
#
# @dataclass
# class X:
#     uri: str
#     product_id: str
#     headers: [astropy.io.fits.Header]
#     blueprint: ObsBlueprint
#     file_info: FileInfo
#
#
# def _get_uris(args):
#     result = []
#     if args.lineage:
#         for ii in args.lineage:
#             ignore, uri = mc.decompose_lineage(ii)
#             result.append(uri)
#     elif args.local:
#         for ii in args.local:
#             dao_name = dn.DAOName(ii)
#             result.append(dao_name.file_uri)
#     else:
#         raise mc.CadcException(f'Could not define uri from these args {args}')
#     return result
#
#
# def _get_uris_with_client(args, header_client):
#     result = []
#     if args.local:
#         for ii in args.local:
#             dao_name = dn.DAOName(ii)
#             headers = header_client.get_head(ii)
#             y = X(dao_name.file_uri, dao_name.product_id, headers, None, None)
#             result.append(y)
#     elif args.lineage:
#         for ii in args.lineage:
#             product_id, uri = mc.decompose_lineage(ii)
#             headers = header_client.get_head(uri)
#             y = X(uri, product_id, headers, None, None)
#             result.append(y)
#     else:
#         raise mc.CadcException(f'Could not define uri from these args {args}')
#     return result
#
#
# def to_caom2():
#     args = get_gen_proc_arg_parser().parse_args()
#     # only need the telescopes instances for this invocation
#     telescopes.current = {}
#     uris = _get_uris(args)
#     blueprints = _build_blueprints(uris)
#     return gen_proc(args, blueprints)
#
#
# def to_caom2_with_client(header_client):
#     args = get_gen_proc_arg_parser().parse_args()
#     # only need the telescopes instances for this invocation
#     # telescopes.current = {}
#     current_telescopes = {}
#     xs = _get_uris_with_client(args, header_client)
#     _build_blueprints_with_client(xs, current_telescopes)
#     kwargs = {'blueprints': xs}
#     result = gen_proc_with_headers(args, **kwargs)
#     return result


def dao_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        # result = to_caom2()
        result = None
        logging.debug(f'Done {APPLICATION} processing.')
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)
