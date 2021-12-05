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

from caom2 import DataProductType, ProductType, TypedSet, ObservationURI
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from caom2utils import ObsBlueprint, GenericParser, FitsParser
from caom2utils import update_artifact_meta_no_client
from dao2caom2 import telescopes, dao_name


class Fits2caom2Visitor:
    def __init__(self, observation, **kwargs):
        self._observation = observation
        self._storage_name = kwargs.get('storage_name')
        self._headers_dict = kwargs.get('headers_dict')
        self._dump_config = False
        self._logger = logging.getLogger(self.__class__.__name__)

    def visit(self):
        for uri, values in self._headers_dict.items():
            headers = values[0]
            file_info = values[1]
            telescope_data = telescopes.factory_client(uri, headers)
            blueprint = ObsBlueprint(instantiated_class=telescope_data)
            telescope_data.configure_axes(blueprint)
            telescope_data.accumulate_bp(blueprint)
            self._logger.error(type(telescope_data))

            if self._dump_config:
                print(f'Blueprint for {uri}: {blueprint}')

            if len(headers) == 0:
                parser = GenericParser(blueprint, uri)
            else:
                parser = FitsParser(headers, blueprint, uri)

            parser.augment_observation(
                observation=self._observation,
                artifact_uri=uri,
                product_id=self._storage_name.product_id,
            )

            self._update(headers, telescope_data, file_info)
        return {'observation': 1}

    def _update(self, headers, telescope, file_info):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n
        relationship between TDM attributes and CAOM attributes). Must have this
        signature for import_module loading and execution.
        """
        logging.debug('Begin update.')
        # correct the *_axis values
        for plane in self._observation.planes.values():
            for artifact in plane.artifacts.values():
                if artifact.uri.replace('.gz', '') != self._storage_name.file_uri.replace(
                        '.gz', ''
                ):
                    logging.error(f'artifact uri {artifact.uri} file_uri {self._storage_name.file_uri}')
                    continue

                update_artifact_meta_no_client(artifact, file_info)
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        time_delta = telescope.get_time_axis_delta(ext=0)
                        cc.undo_astropy_cdfix_call(chunk, time_delta)

                        if self._storage_name.file_name.startswith('d'):
                            if plane.data_product_type == DataProductType.SPECTRUM:
                                if (
                                        dao_name.DAOName.is_unprocessed_reticon(artifact.uri)
                                        or dao_name.DAOName.is_derived(artifact.uri)
                                        and self._observation.type == 'flat'
                                ):
                                    cc.reset_energy(chunk)
                                if (
                                        artifact.product_type != ProductType.SCIENCE
                                ):
                                    if self._observation.type == 'dark':
                                        chunk.position_axis_1 = 3
                                        chunk.position_axis_2 = 4
                                    else:
                                        cc.reset_position(chunk)
                                    # no energy for calibration?
                                    if self._observation.type not in [
                                        'flat',
                                        'comparison',
                                        'dark',
                                    ]:
                                        cc.reset_energy(chunk)
                            else:  # DataProductType.IMAGE
                                if dao_name.DAOName.override_provenance(artifact.uri):
                                    plane.provenance.producer = 'Spaceguard_C'
                                # no observable axis when image
                                cc.reset_observable(chunk)
                                if (
                                        artifact.product_type
                                        == ProductType.CALIBRATION
                                ):
                                    if self._observation.type != 'dark':
                                        cc.reset_position(chunk)
                                    if self._observation.type not in ['flat', 'dark']:
                                        cc.reset_energy(chunk)
                            if (
                                    chunk.energy is not None
                                    and not dao_name.DAOName.is_processed(artifact.uri)
                                    and headers[0].get('WAVELENG') is None
                            ):
                                # DB 16-02-21/04-03-21
                                #  If WAVELENG isn’t present then all energy
                                #  metadata should be ignored (spectra and images)
                                cc.reset_energy(chunk)

                            # WCS axis wisdom from Pat:
                            #
                            # In general, assigning axis indices above the value of
                            # naxis is allowed but more or less pointless. The
                            # only use case that would justify it is that in a FITS
                            # file there could be a header with NAXIS=2 and
                            # WCSAXES=4 which would tell the fits reader to look
                            # for CTYPE1 through 4 and axes 3 and 4 are metadata.
                            # Assign those values to Chunk only if you care about
                            # capturing that the extra wcs metadata was really in
                            # the fits header and so the order could be preserved;
                            # in general do not assign the 3 and 4.

                            naxis = headers[0].get('NAXIS')
                            naxis1 = headers[0].get('NAXIS1')
                            naxis2 = headers[0].get('NAXIS2')
                            chunk.naxis = None
                            chunk.position_axis_1 = None
                            chunk.position_axis_2 = None
                            chunk.energy_axis = None
                            chunk.observable_axis = None
                            chunk.time_axis = None
                            if naxis is not None:
                                if (
                                        naxis1 is not None
                                        and naxis2 is not None
                                        and naxis == 2
                                        and chunk.position is not None
                                        and plane.data_product_type
                                        is DataProductType.IMAGE
                                ):
                                    chunk.naxis = 2
                                    chunk.position_axis_1 = 1
                                    chunk.position_axis_2 = 2
                                if (
                                        naxis1 is not None
                                        and naxis == 1
                                        and chunk.energy is not None
                                ):
                                    chunk.naxis = 1
                                    chunk.energy_axis = 1
                        else:
                            chunk.energy_axis = None

            if plane.product_id != self._storage_name.product_id:
                continue

            # provenance: inputs vs members
            #
            # DB - 29-04-20
            # The inconsistencies are consistent for both telescope for the
            # derived observations: the processed, co-added flats (files with F
            # suffix) and co-added biases (files with B suffix).  These should
            # have the ‘members’ list added to the inputs.  From the definition of
            # provenance:inputs I’m assuming for the science observations,
            # processed comparison arcs, and processed flats having a composite
            # flat and/or bias observation as an ‘input’ is okay rather than
            # breaking these down into their individual members (since those
            # derived observations will all be available in the archive with
            # proper provenance provided).

            if self._observation.type == 'flat' and cc.is_composite(headers, 'FLAT_'):
                cc.update_plane_provenance(
                    plane,
                    headers,
                    'FLAT_',
                    'DAO',  # TODO
                    _repair_provenance_value,
                    self._observation.observation_id,
                )
            elif self._observation.type == 'bias' and cc.is_composite(headers, 'ZERO_'):
                cc.update_plane_provenance(
                    plane,
                    headers,
                    'ZERO_',
                    'DAO',  # TODO
                    _repair_provenance_value,
                    self._observation.observation_id,
                )

            if dao_name.DAOName.is_processed(self._storage_name.file_uri):
                self._update_plane_provenance(plane, headers)

            if cc.is_composite(headers, 'FLAT_') or cc.is_composite(
                    headers, 'ZERO_'
            ):
                self._update_observation_members()

        logging.debug('Done update.')

    def _update_observation_members(self):
        """
        Must filter results because:
        DB - 11-06-20
        For the spectra there is a minor issue with members for master flat,
        *_F, observations.  The master bias used in the processing, the *_B.fits
        file, shouldn’t be a member for the master flats.

        The master bias is in the list of inputs though:  Inputs for master flat
        are the unprocessed flats and the master bias.  The master bias is
        subtracted pixel-by-pixel from each unprocessed flat as part of the
        processing before the flats are then co-added.

        The composite/derived master flats (F) and master biases (B) should
        never be members.  At least for any processing that is currently being
        done.  For now the only members should those given by the
        NCOMBINE x ZERO_# or FLAT_# keyword values.
        """

        def filter_fun(x):
            result = True
            if DAOName.is_master_flat(self._observation.observation_id):
                if DAOName.is_master_bias(x.get_observation_uri().uri):
                    result = False
            return result

        inputs = []
        members_inputs = TypedSet(ObservationURI,)
        for plane in self._observation.planes.values():
            if (
                    plane.provenance is not None
                    and plane.provenance.inputs is not None
            ):
                inputs = filter(filter_fun, plane.provenance.inputs)

        for entry in inputs:
            members_inputs.add(entry.get_observation_uri())
            logging.debug(f'Adding Observation URI {entry.get_observation_uri()}')
        mc.update_typed_set(self._observation.members, members_inputs)

    def _update_plane_provenance(self, plane, headers):
        logging.debug(
            f'Begin _update_plane_provenance for {plane.product_id} with'
            f'observation type: {self._observation.type}.'
        )
        if self._observation.type in ['object', 'flat', 'comparison']:
            f_name = headers[0].get('BIAS')
            if f_name is not None:
                bias_name = dao_name.DAOName(f_name)
                ignore, plane_uri = cc.make_plane_uri(bias_name.obs_id, bias_name.product_id,
                    'DAO')  # TODO
                plane.provenance.inputs.add(plane_uri)
        if self._observation.type in ['object', 'comparison']:
            f_name = headers[0].get('FLAT')
            if f_name is not None:
                flat_name = dao_name.DAOName(f_name)
                ignore, plane_uri = cc.make_plane_uri(flat_name.obs_id, flat_name.product_id,
                    'DAO')  # TODO
                plane.provenance.inputs.add(plane_uri)
            # referral to raw plane
            ignore, plane_uri = cc.make_plane_uri(
                self._observation.observation_id, self._observation.observation_id,
                'DAO',  # TODO
            )
            plane.provenance.inputs.add(plane_uri)
        if self._observation.type == 'object':
            f_name = headers[0].get('DCLOG1')
            if f_name is not None:
                ref_spec1_name = dao_name.DAOName(f_name.split()[2])
                ignore, plane_uri = cc.make_plane_uri(
                    ref_spec1_name.obs_id, ref_spec1_name.product_id,
                    'DAO',  # TODO
                )
                plane.provenance.inputs.add(plane_uri)
            if headers[0].get('DCLOG2') is not None:
                ref_spec1_name = dao_name.DAOName(headers[0].get('DCLOG2').split()[2])
                ignore, plane_uri = cc.make_plane_uri(
                    ref_spec1_name.obs_id, ref_spec1_name.product_id
                )
                plane.provenance.inputs.add(plane_uri)
        logging.debug(f'End _update_plane_provenance.')


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
    dn = dao_name.DAOName(value)
    prov_prod_id = dn.product_id
    prov_obs_id = dn.obs_id
    logging.debug(f'End _repair_provenance_value')
    return prov_obs_id, prov_prod_id


def visit(observation, **kwargs):
    s = Fits2caom2Visitor(observation, **kwargs)
    return s.visit()
