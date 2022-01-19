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

from astropy.coordinates import SkyCoord, FK5
import astropy.units as u

from caom2utils import update_artifact_meta
from caom2 import CalibrationLevel, DataProductType, TargetType, ProductType
from caom2 import ObservationIntentType, TypedSet, ObservationURI
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from dao2caom2.dao_name import DAOName, COLLECTION


__all__ = ['DAOTelescopeMapping']


APPLICATION = 'dao2caom2'


class DAOTelescopeMapping(cc.TelescopeMapping):
    """
    This class is the Spectrum implementation.

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

                                       # DB 04-02-21
                                       If no 'WAVELENG' no energy

                                       CRVAL   'WAVELENG'
                                       no `DELTA_WL', 'REFPIXEL'
                                       CDELT   'DISPERSI' * 15.0 * xbin/1000.0
                                       CRPIX   'DATASEC' + math

    processed DAO spectrum
    dao_[cr]\d{3}_\d{4}_\d{6}_[evBF], obs.type in ['object', 'comparison']:
                                       CRVAL   'CRVAL1'
                                       CDELT   'CDELT1'
                                       CRPIX   'CRPIX1'
                                       RP      CRVAL/(2.5*CDELT)

    """
    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def update(self, observation, file_info, caom_repo_client):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n
        relationship between TDM attributes and CAOM attributes). Must have this
        signature for import_module loading and execution.
        """
        logging.debug('Begin update.')
        # correct the *_axis values
        for plane in observation.planes.values():
            for artifact in plane.artifacts.values():
                if (
                    artifact.uri.replace('.gz', '')
                    != self._storage_name.file_uri.replace('.gz', '')
                ):
                    self._logger.debug(
                        f'Skipping artifact {self._storage_name.file_uri}'
                    )
                    continue

                update_artifact_meta(artifact, file_info)
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        time_delta = self.get_time_axis_delta(ext=0)
                        cc.undo_astropy_cdfix_call(chunk, time_delta)

                        if self._storage_name.file_name.startswith('d'):
                            if plane.data_product_type == DataProductType.SPECTRUM:
                                if (
                                    DAOName.is_unprocessed_reticon(artifact.uri)
                                    or DAOName.is_derived(artifact.uri)
                                    and observation.type == 'flat'
                                ):
                                    cc.reset_energy(chunk)
                                if (
                                    artifact.product_type != ProductType.SCIENCE
                                ):
                                    if observation.type == 'dark':
                                        chunk.position_axis_1 = 3
                                        chunk.position_axis_2 = 4
                                    else:
                                        cc.reset_position(chunk)
                                    # no energy for calibration?
                                    if observation.type not in [
                                        'flat',
                                        'comparison',
                                        'dark',
                                    ]:
                                        cc.reset_energy(chunk)
                            else:  # DataProductType.IMAGE
                                if DAOName.override_provenance(artifact.uri):
                                    plane.provenance.producer = 'Spaceguard_C'
                                # no observable axis when image
                                cc.reset_observable(chunk)
                                if (
                                        artifact.product_type
                                        == ProductType.CALIBRATION
                                ):
                                    if observation.type != 'dark':
                                        cc.reset_position(chunk)
                                    if observation.type not in ['flat', 'dark']:
                                        cc.reset_energy(chunk)
                            if (
                                chunk.energy is not None
                                and not DAOName.is_processed(artifact.uri)
                                and self._headers[0].get('WAVELENG') is None
                            ):
                                # DB 16-02-21/04-03-21
                                #  If WAVELENG isn’t present then all energy
                                #  metadata should be ignored (spectra and
                                #  images)
                                cc.reset_energy(chunk)

                            # WCS axis wisdom from Pat:
                            #
                            # In general, assigning axis indices above the
                            # value of naxis is allowed but more or less
                            # pointless. The only use case that would justify
                            # it is that in a FITS file there could be a
                            # header with NAXIS=2 and WCSAXES=4 which would
                            # tell the fits reader to look for CTYPE1 through
                            # 4 and axes 3 and 4 are metadata. Assign those
                            # values to Chunk only if you care about capturing
                            # that the extra wcs metadata was really in the
                            # fits header and so the order could be preserved;
                            # in general do not assign the 3 and 4.

                            naxis = self._headers[0].get('NAXIS')
                            naxis1 = self._headers[0].get('NAXIS1')
                            naxis2 = self._headers[0].get('NAXIS2')
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
            # derived observations: the processed, co-added flats (files with
            # F suffix) and co-added biases (files with B suffix).  These
            # should have the ‘members’ list added to the inputs.  From the
            # definition of provenance:inputs I’m assuming for the science
            # observations, processed comparison arcs, and processed flats
            # having a composite flat and/or bias observation as an ‘input’ is
            # okay rather than breaking these down into their individual
            # members (since those derived observations will all be available
            # in the archive with proper provenance provided).

            if (
                observation.type == 'flat'
                and cc.is_composite(self._headers, 'FLAT_')
            ):
                cc.update_plane_provenance(
                    plane,
                    self._headers,
                    'FLAT_',
                    'DAO',  # TODO
                    _repair_provenance_value,
                    observation.observation_id,
                )
            elif (
                observation.type == 'bias'
                and cc.is_composite(self._headers, 'ZERO_')
            ):
                cc.update_plane_provenance(
                    plane,
                    self._headers,
                    'ZERO_',
                    'DAO',  # TODO
                    _repair_provenance_value,
                    observation.observation_id,
                )

            if DAOName.is_processed(self._storage_name.file_uri):
                self._update_plane_provenance(
                    observation, plane, self._headers
                )

            if cc.is_composite(self._headers, 'FLAT_') or cc.is_composite(
                    self._headers, 'ZERO_'
            ):
                self._update_observation_members(observation)

        logging.debug('Done update.')
        return observation

    def _update_observation_members(self, observation):
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
            if DAOName.is_master_flat(observation.observation_id):
                if DAOName.is_master_bias(x.get_observation_uri().uri):
                    result = False
            return result

        inputs = []
        members_inputs = TypedSet(ObservationURI,)
        for plane in observation.planes.values():
            if (
                    plane.provenance is not None
                    and plane.provenance.inputs is not None
            ):
                inputs = filter(filter_fun, plane.provenance.inputs)

        for entry in inputs:
            members_inputs.add(entry.get_observation_uri())
            logging.debug(f'Adding Observation URI {entry.get_observation_uri()}')
        mc.update_typed_set(observation.members, members_inputs)

    def _update_plane_provenance(self, observation, plane, headers):
        self._logger.debug(
            f'Begin _update_plane_provenance for {plane.product_id} with'
            f'observation type: {observation.type}.'
        )
        if observation.type in ['object', 'flat', 'comparison']:
            f_name = headers[0].get('BIAS')
            if f_name is not None:
                bias_name = DAOName(f_name)
                ignore, plane_uri = cc.make_plane_uri(
                    bias_name.obs_id, bias_name.product_id, 'DAO'
                )  # TODO
                plane.provenance.inputs.add(plane_uri)
        if observation.type in ['object', 'comparison']:
            f_name = headers[0].get('FLAT')
            if f_name is not None:
                flat_name = DAOName(f_name)
                ignore, plane_uri = cc.make_plane_uri(
                    flat_name.obs_id, flat_name.product_id, 'DAO'
                )  # TODO
                plane.provenance.inputs.add(plane_uri)
            # referral to raw plane
            ignore, plane_uri = cc.make_plane_uri(
                observation.observation_id, observation.observation_id,
                'DAO',  # TODO
            )
            plane.provenance.inputs.add(plane_uri)
        if observation.type == 'object':
            f_name = headers[0].get('DCLOG1')
            if f_name is not None:
                ref_spec1_name = DAOName(f_name.split()[2])
                ignore, plane_uri = cc.make_plane_uri(
                    ref_spec1_name.obs_id, ref_spec1_name.product_id,
                    'DAO',  # TODO
                )
                plane.provenance.inputs.add(plane_uri)
            if headers[0].get('DCLOG2') is not None:
                ref_spec1_name = DAOName(headers[0].get('DCLOG2').split()[2])
                ignore, plane_uri = cc.make_plane_uri(
                    ref_spec1_name.obs_id,
                    ref_spec1_name.product_id,
                    COLLECTION,
                )
                plane.provenance.inputs.add(plane_uri)
        self._logger.debug(f'End _update_plane_provenance.')

    def _get_wavelength(self, ext):
        return mc.to_float(self._headers[ext].get('WAVELENG'))

    def configure_axes(self, bp):
        bp.configure_position_axes((1, 2))
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)
        bp.configure_observable_axis(5)

    def accumulate_blueprint(self, bp, application=None):
        self.configure_axes(bp)
        super().accumulate_blueprint(bp, APPLICATION)
        bp.set(
            'Chunk.energy.axis.function.delta',
            'get_energy_axis_function_delta()',
        )
        bp.set(
            'Chunk.energy.axis.function.naxis',
            'get_energy_axis_function_naxis()',
        )
        bp.set(
            'Chunk.energy.axis.function.refCoord.pix',
            'get_energy_axis_function_refcoord_pix()',
        )
        bp.clear('Chunk.energy.axis.function.refCoord.val')
        bp.add_fits_attribute(
            'Chunk.energy.axis.function.refCoord.val', 'WAVELENG'
        )

        bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
        bp.set(
            'Chunk.position.axis.function.cd11',
            'get_position_function_cd11()',
        )
        bp.set(
            'Chunk.position.axis.function.cd22',
            'get_position_function_cd22()',
        )
        bp.set(
            'Chunk.position.axis.function.cd12',
            'get_position_function_cd12()',
        )
        bp.set(
            'Chunk.position.axis.function.cd21',
            'get_position_function_cd21()',
        )
        bp.set('Observation.intent', 'get_obs_intent()')
        bp.clear('Observation.metaRelease')
        # from dao2caom2.config
        bp.add_fits_attribute('Observation.metaRelease', 'DATE-OBS')

        bp.set('Observation.target.type', TargetType.OBJECT)
        bp.clear('Observation.target.moving')
        bp.set_default('Observation.target.moving', 'false')
        bp.clear('Observation.target.standard')
        bp.set_default('Observation.target.standard', 'false')

        bp.clear('Observation.proposal.id')
        bp.add_fits_attribute('Observation.proposal.id', 'DAOPRGID')
        bp.clear('Observation.proposal.pi')
        bp.add_fits_attribute('Observation.proposal.pi', 'PINAME')

        bp.clear('Observation.environment.humidity')
        bp.add_fits_attribute('Observation.environment.humidity', 'REL_HUMI')
        bp.clear('Observation.environment.photometric')
        bp.set_default('Observation.environment.photometric', 'false')

        bp.set('Plane.dataProductType', DataProductType.SPECTRUM)
        bp.set('Plane.calibrationLevel', CalibrationLevel.RAW_STANDARD)
        bp.clear('Plane.metaRelease')
        # from dao2caom2.config
        bp.add_fits_attribute('Plane.metaRelease', 'DATE-OBS')

        bp.clear('Plane.provenance.lastExecuted')
        bp.add_fits_attribute('Plane.provenance.lastExecuted', 'IRAF-TLM')
        bp.set('Plane.provenance.project', 'DAO Science Archive')
        bp.clear('Plane.provenance.name')
        bp.add_fits_attribute('Plane.provenance.name', 'PROCNAME')
        bp.set_default('Plane.provenance.name', 'DAO unprocessed data')
        bp.set('Plane.provenance.producer', 'NRC Herzberg')
        bp.set(
            'Plane.provenance.reference',
            'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/dao/',
        )
        bp.clear('Plane.provenance.version')
        bp.add_fits_attribute('Plane.provenance.version', 'PROCVERS')
        bp.set_default('Plane.provenance.version', None)

        bp.set('Artifact.productType', 'get_artifact_product_type()')
        bp.set('Chunk.time.exposure', 'get_time_exposure(header)')
        bp.set('Chunk.time.resolution', 'get_time_resolution(header)')

        bp.set('Chunk.observable.axis.axis.ctype', 'FLUX')
        bp.set('Chunk.observable.axis.axis.cunit', 'COUNTS')
        bp.set('Chunk.observable.axis.function.refCoord.pix', 1)

        bp.set(
            'Chunk.energy.resolvingPower',
            'get_energy_resolving_power()',
        )
        bp.clear('Chunk.energy.bandpassName')
        bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')

        bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        bp.set(
            'Chunk.position.axis.function.refCoord.coord1.pix',
            'get_position_function_coord1_pix()',
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord1.val',
            'get_position_function_coord1_val(header)',
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord2.pix',
            'get_position_function_coord2_pix()',
        )
        bp.set(
            'Chunk.position.axis.function.refCoord.coord2.val',
            'get_position_function_coord2_val(header)',
        )
        self._accumulate_common_bp(bp)

    def get_artifact_product_type(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        if obs_type == 'object':
            product_type = ProductType.SCIENCE
        else:
            product_type = ProductType.CALIBRATION
        return product_type

    def get_energy_axis_function_delta(self, ext):
        wavelength = self._get_wavelength(ext)
        cdelt = self._headers[ext].get('DELTA_WL')
        if wavelength is None:
            cdelt = None
        else:
            if cdelt is None:
                dispersion = self._headers[ext].get('DISPERSI')
                dispaxis = self._get_dispaxis(ext)
                if dispaxis == 1:
                    xbin = mc.to_float(self._headers[ext].get('XBIN'))
                else:
                    xbin = mc.to_float(self._headers[ext].get('YBIN'))
                cdelt = dispersion * 15.0 * xbin / 1000.0
        return cdelt

    def get_energy_axis_function_refcoord_pix(self, ext):
        wavelength = self._get_wavelength(ext)
        if wavelength is None:
            crpix = None
        else:
            crpix = self._headers[ext].get('REFPIXEL')
            if crpix is None:
                temp = self._headers[ext].get('DATASEC')
                if temp is not None:
                    datasec = re.sub(
                        r'(\[)(\d+:\d+,\d+:\d+)(\])', r'\g<2>', temp
                    )
                    (dx, dy) = datasec.split(',')
                    (xl, xh) = dx.split(':')
                    (yl, yh) = dy.split(':')
                    dispaxis = self._get_dispaxis(ext)
                    if dispaxis == 1:
                        crpix = (int(xh) - int(xl)) / 2.0 + int(xl)
                    else:
                        crpix = (int(yh) - int(yl)) / 2.0 + int(yl)
        return crpix

    def get_energy_axis_function_naxis(self, ext):
        return 1

    def get_energy_resolving_power(self, ext):
        numerator = self._headers[ext].get('WAVELENG')
        denominator = self.get_energy_axis_function_delta(ext)
        return numerator / (2.5 * denominator)

    def get_geo(self):
        # DB 10-09-20
        # Google Maps to give you latitude/longitude if desired. 48.519497
        # and -123.416502.  Not sure of the elevation.
        return ac.get_location(48.519497, -123.416502, 210.0)

    def get_position_function_cd11(self, ext):
        result = self.get_position_function_cd22(ext)
        if result is not None:
            result = (-1) * result
        return result

    def get_position_function_cd12(self, ext):
        return self.get_position_function_cd21(ext)

    def get_position_function_cd21(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            result = 0.0
        return result

    def get_position_function_cd22(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            # DB - set entrance aperture to a fixed 5" by 5" because of lack
            # of detailed information
            result = 0.001388
        return result

    def get_position_function_coord1_pix(self, ext):
        return self.get_position_function_coord2_pix(ext)

    def get_position_function_coord2_pix(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        result = None
        if obs_type in ['dark', 'object']:
            result = 1.0
        return result

    def get_telescope_name(self, ign):
        return None

    def get_time_axis_val(self, ext):
        return ac.get_datetime(self._headers[ext].get('DATE-OBS'))

    def get_geo_x(self, ext):
        x, ign1, ign2 = self.get_geo()
        return x

    def get_geo_y(self, ext):
        x, y, ign2 = self.get_geo()
        return y

    def get_geo_z(self, ext):
        x, ign1, z = self.get_geo()
        return z

    def get_obs_intent(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        if obs_type == 'object':
            intent = ObservationIntentType.SCIENCE
        else:
            intent = ObservationIntentType.CALIBRATION
        return intent

    def get_time_exposure(self, ext):
        exptime = mc.to_float(self._headers[ext].get('EXPTIME'))
        ncombine = mc.to_float(self._headers[ext].get('NCOMBINE'))
        if ncombine is not None:
            # DB - approximation of exposure time for products (assume identical
            # EXPTIME)
            exptime *= ncombine
        return exptime

    def get_time_resolution(self, ext):
        exptime = mc.to_float(self._headers[ext].get('EXPTIME'))
        ncombine = mc.to_float(self._headers[ext].get('NCOMBINE'))
        if ncombine is None:
            ncombine = 1
        else:
            exptime = exptime * ncombine
        return exptime / ncombine

    def get_position_function_coord1_val(self, ext):
        ra, ignore_dec = self._get_position(ext)
        return ra

    def get_position_function_coord2_val(self, ext):
        ignore_ra, dec = self._get_position(ext)
        return dec

    def _get_position(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        ra_deg = None
        dec_deg = None
        if obs_type in ['comparison', 'dark', 'object']:
            if self._headers[ext].get('EQUINOX') is not None:
                # DB - 11-09-19 - precession with astropy
                ra = self._headers[ext].get('RA', 0)
                dec = self._headers[ext].get('DEC', 0)
                equinox = f'J{self._headers[ext].get("EQUINOX")}'
                fk5 = FK5(equinox=equinox)
                coord = SkyCoord(
                    f'{ra} {dec}', unit=(u.hourangle, u.deg), frame=fk5
                )
                j2000 = FK5(equinox='J2000')
                result = coord.transform_to(j2000)
                ra_deg = result.ra.degree
                dec_deg = result.dec.degree
        return ra_deg, dec_deg

    def get_time_axis_delta(self, ext):
        exptime = self.get_time_exposure(ext)
        return exptime / (24.0 * 3600.0)

    def _accumulate_common_bp(self, bp):
        meta_producer = mc.get_version('dao2caom2')
        bp.set('Observation.metaProducer', meta_producer)
        bp.set('Plane.metaProducer', meta_producer)
        bp.set('Artifact.metaProducer', meta_producer)
        bp.set('Chunk.metaProducer', meta_producer)

        bp.clear('Observation.algorithm.name')

        bp.set('Observation.telescope.name', 'get_telescope_name()')
        bp.set('Observation.telescope.geoLocationX', 'get_geo_x()')
        bp.set('Observation.telescope.geoLocationY', 'get_geo_y()')
        bp.set('Observation.telescope.geoLocationZ', 'get_geo_z()')

        bp.set('Artifact.releaseType', 'data')

        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.function.naxis', '1')
        bp.set(
            'Chunk.time.axis.function.delta', 'get_time_axis_delta()'
        )
        bp.set('Chunk.time.axis.function.refCoord.pix', '0.5')
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_axis_val(params)'
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
        if DAOName.is_derived(self._storage_name.file_uri):
            bp.set('DerivedObservation.members', [])
            bp.add_fits_attribute('Observation.algorithm.name', 'PROCNAME')


class Dao12Metre(DAOTelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def get_geo(self):
        return ac.get_location(48.52092, -123.42006, 225.0)

    def get_telescope_name(self, ext):
        return 'DAO 1.2-m'


class Dao18Metre(DAOTelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def get_geo(self):
        return ac.get_location(48.51967, -123.41833, 232.0)

    def get_telescope_name(self, ext):
        return 'DAO 1.8-m'


class SkyCam(DAOTelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def configure_axes(self, bp):
        # DB - 10-07-20
        # https://github.com/opencadc-metadata-curation/dao2caom2/issues/10
        bp.configure_time_axis(1)
        bp.configure_observable_axis(2)
        bp.configure_energy_axis(3)

    def accumulate_blueprint(self, bp, application=None):
        self.configure_axes(bp)
        meta_producer = mc.get_version(APPLICATION)
        bp.set('Observation.metaProducer', meta_producer)
        bp.set('Plane.metaProducer', meta_producer)
        bp.set('Artifact.metaProducer', meta_producer)
        bp.set('Chunk.metaProducer', meta_producer)

        bp.set('Observation.metaRelease', 'get_release_date()')
        bp.set('Observation.intent', ObservationIntentType.CALIBRATION)
        bp.set('Observation.instrument.name', 'Sky Camera')
        bp.set('Plane.calibrationLevel', 1)
        bp.set('Plane.dataProductType', DataProductType.IMAGE)
        bp.set('Plane.dataRelease', 'get_release_date()')
        bp.set('Plane.metaRelease', 'get_release_date()')
        bp.set('Plane.provenance.project', 'DAO Science Archive')
        bp.set(
            'Plane.provenance.producer',
            'NRC Herzberg Astronomy and Astrophysics Research Centre',
        )
        bp.set('Plane.provenance.name', 'DAO Sky Camera image')
        bp.set(
            'Plane.provenance.reference',
            'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/dao/',
        )
        bp.set('Artifact.productType', ProductType.CALIBRATION)
        bp.set('Artifact.releaseType', 'data')

        bp.set('Chunk.energy.axis.function.delta', 3000.0)
        bp.set('Chunk.energy.axis.function.naxis', 1)
        bp.set('Chunk.energy.axis.function.refCoord.pix', 0.5)
        bp.set('Chunk.energy.axis.function.refCoord.val', 4000.0)
        bp.set('Chunk.energy.resolvingPower', 5500.0 / 3000.0)

        bp.add_fits_attribute('Chunk.time.exposure', 'EXPTIME')
        bp.add_fits_attribute('Chunk.time.resolution', 'EXPTIME')
        self._accumulate_common_bp(bp)

    def get_release_date(self, ext):
        return ac.get_datetime(self._headers[ext].get('CLOCKVAL'))

    def get_telescope_name(self, ext):
        return 'DAO Skycam'

    def get_time_axis_val(self, ext):
        return ac.get_datetime(self._headers[ext].get('CLOCKVAL'))


class Imaging(DAOTelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def configure_axes(self, bp):
        bp.configure_position_axes((1, 2))
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, application)
        bp.set('Observation.target.type', TargetType.FIELD)
        bp.set('Plane.dataProductType', DataProductType.IMAGE)
        bp.clear('Chunk.energy.axis.function.delta')
        bp.add_fits_attribute('Chunk.energy.axis.function.delta', 'BANDPASS')
        bp.set('Chunk.energy.axis.function.refCoord.pix', 1.0)
        bp.clear('Chunk.position.axis.function.dimension.naxis1')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.dimension.naxis1', 'NAXIS1'
        )
        bp.clear('Chunk.position.axis.function.dimension.naxis2')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.dimension.naxis2', 'NAXI2'
        )

    def _get_position_by_scale_size_bin(self, ext):
        result = None
        platescale = mc.to_float(self._headers[ext].get('PLTSCALE'))
        pixsize = mc.to_float(self._headers[ext].get('PIXSIZE'))
        xbin = mc.to_float(self._headers[ext].get('XBIN'))
        if platescale is not None and pixsize is not None and xbin is not None:
            result = platescale * pixsize * xbin / 3600000.0
        return result

    def get_energy_resolving_power(self, ext):
        wavelength = self._headers[ext].get('WAVELENG')
        bandpass = self._headers[ext].get('BANDPASS')
        return wavelength / bandpass

    def get_position_function_cd11(self, ext):
        return self._get_position_by_scale_size_bin(ext)

    def get_position_function_cd22(self, ext):
        return self._get_position_by_scale_size_bin(ext)

    def get_position_function_coord1_pix(self, ext):
        return self._headers[ext].get('NAXIS1') / 2.0

    def get_position_function_coord2_pix(self, ext):
        return self._headers[ext].get('NAXIS2') / 2.0


class ProcessedImage(Imaging):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, application)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        bp.clear('Chunk.position.axis.function.cd11')
        bp.add_fits_attribute('Chunk.position.axis.function.cd11', 'CD1_1')
        bp.clear('Chunk.position.axis.function.cd22')
        bp.add_fits_attribute('Chunk.position.axis.function.cd22', 'CD2_2')
        bp.clear('Chunk.position.axis.function.cd12')
        bp.add_fits_attribute('Chunk.position.axis.function.cd12', 'CD1_2')
        bp.clear('Chunk.position.axis.function.cd21')
        bp.add_fits_attribute('Chunk.position.axis.function.cd21', 'CD2_1')
        bp.clear('Chunk.position.axis.function.refCoord.coord1.pix')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord1.pix', 'CRPIX1'
        )
        bp.clear('Chunk.position.axis.function.refCoord.coord2.pix')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord2.pix', 'CRPIX2'
        )


class ProcessedSpectrum(DAOTelescopeMapping):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def configure_axes(self, bp):
        bp.configure_position_axes((2, 3))
        bp.configure_time_axis(4)
        bp.configure_energy_axis(1)
        bp.configure_observable_axis(5)

    def accumulate_blueprint(self, bp, application=None):
        super().accumulate_blueprint(bp, application)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        if DAOName.is_derived(self._storage_name.file_uri):
            # original dao2caom2.py, l392, l400
            bp.set('Observation.target.type', None)
        bp.clear('Chunk.energy.axis.function.delta')
        bp.add_fits_attribute('Chunk.energy.axis.function.delta', 'CDELT1')
        bp.clear('Chunk.energy.axis.function.naxis')
        bp.add_fits_attribute('Chunk.energy.axis.function.naxis', 'NAXIS1')
        bp.clear('Chunk.energy.axis.function.refCoord.pix')
        bp.add_fits_attribute(
            'Chunk.energy.axis.function.refCoord.pix', 'CRPIX1'
        )
        bp.clear('Chunk.energy.axis.function.refCoord.val')
        bp.add_fits_attribute(
            'Chunk.energy.axis.function.refCoord.val', 'CRVAL1'
        )
        bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 1)


class Dao12MetreImage(Dao12Metre, Imaging):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)


class Dao12MetreProcessedImage(Dao12MetreImage, ProcessedImage):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)


class Dao12MetreSpectrum(Dao12Metre):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def _get_dispaxis(self, ext):
        return self._headers[ext].get('DISPAXIS', 2)

    def get_energy_axis_function_naxis(self, ext):
        dispaxis = self._get_dispaxis(ext)
        return self._headers[ext].get(f'NAXIS{dispaxis}')


class Dao12MetreProcessedSpectrum(Dao12MetreSpectrum, ProcessedSpectrum):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def get_energy_resolving_power(self, ext):
        obs_type = self._headers[ext].get('OBSTYPE')
        if obs_type in ['comparison', 'object']:
            numerator = self._headers[ext].get('CRVAL1')
            denominator = self._headers[ext].get('CDELT1')
        else:
            numerator = self._headers[ext].get('WAVELENG')
            denominator = self.get_energy_axis_function_delta(ext)
        return numerator / (2.5 * denominator)


class Dao18MetreImage(Dao18Metre, Imaging):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)


class Dao18MetreProcessedImage(Dao18MetreImage, ProcessedImage):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)


class Dao18MetreSpectrum(Dao18Metre):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def _get_dispaxis(self, ext):
        return self._headers[ext].get('DISPAXIS', 1)

    def get_energy_axis_function_naxis(self, ext):
        dispaxis = self._get_dispaxis(ext)
        return self._headers[ext].get(f'NAXIS{dispaxis}')


class Dao18MetreProcessedSpectrum(Dao18MetreSpectrum, ProcessedSpectrum):

    def __init__(self, storage_name, headers):
        super().__init__(storage_name, headers)

    def get_energy_resolving_power(self, ext):
        obs_type = self._headers[ext].get('OBS_TYPE')
        if obs_type in ['comparison', 'object']:
            numerator = self._headers[ext].get('CRVAL1')
            denominator = self._headers[ext].get('CDELT1')
        else:
            numerator = self._headers[ext].get('WAVELENG')
            denominator = self.get_energy_axis_function_delta(ext)
        return numerator / (2.5 * denominator)


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
    dn = DAOName(value)
    prov_prod_id = dn.product_id
    prov_obs_id = dn.obs_id
    logging.debug(f'End _repair_provenance_value')
    return prov_obs_id, prov_prod_id
