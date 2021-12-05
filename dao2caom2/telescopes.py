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
from caom2 import ObservationIntentType
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from dao2caom2.dao_name import DAOName
from dao2caom2.metadata import DefiningMetadataFinder


__all__ = ['factory', 'factory_client']


class Telescope:
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
    def __init__(self, headers):
        self._uri = None
        self._headers = headers
        self._logger = logging.getLogger(self.__class__.__name__)

    def _get_wavelength(self, ext):
        return mc.to_float(self._headers[ext].get('WAVELENG'))

    def configure_axes(self, bp):
        bp.configure_position_axes((1, 2))
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)
        bp.configure_observable_axis(5)

    def accumulate_bp(self, bp):
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
        if DAOName.is_derived(self._uri):
            bp.set('DerivedObservation.members', 'get_members()')
            bp.add_fits_attribute('Observation.algorithm.name', 'PROCNAME')

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
        from astropy.coordinates import SkyCoord, FK5
        import astropy.units as u
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

    @property
    def uri(self):
        return self._uri

    @uri.setter
    def uri(self, value):
        self._uri = value


class Dao12Metre(Telescope):

    def __init__(self, headers):
        super().__init__(headers)

    def get_geo(self):
        return ac.get_location(48.52092, -123.42006, 225.0)

    def get_telescope_name(self, ext):
        return 'DAO 1.2-m'


class Dao18Metre(Telescope):

    def __init__(self, headers):
        super().__init__(headers)

    def get_geo(self):
        return ac.get_location(48.51967, -123.41833, 232.0)

    def get_telescope_name(self, ext):
        return 'DAO 1.8-m'


class SkyCam(Telescope):

    def __init__(self, headers):
        super().__init__(headers)

    def configure_axes(self, bp):
        # DB - 10-07-20
        # https://github.com/opencadc-metadata-curation/dao2caom2/issues/10
        bp.configure_time_axis(1)
        bp.configure_observable_axis(2)
        bp.configure_energy_axis(3)

    def accumulate_bp(self, bp):
        bp.set('Observation.metaRelease', 'get_skycam_release_date()')
        bp.set('Observation.intent', ObservationIntentType.CALIBRATION)
        bp.set('Observation.instrument.name', 'Sky Camera')
        bp.set('Plane.calibrationLevel', 1)
        bp.set('Plane.dataProductType', DataProductType.IMAGE)
        bp.set('Plane.dataRelease', 'get_skycam_release_date()')
        bp.set('Plane.metaRelease', 'get_skycam_release_date()')
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

        bp.set('Chunk.energy.axis.function.delta', 3000.0)
        bp.set('Chunk.energy.axis.function.naxis', 1)
        bp.set('Chunk.energy.axis.function.refCoord.pix', 0.5)
        bp.set('Chunk.energy.axis.function.refCoord.val', 4000.0)
        bp.set('Chunk.energy.resolvingPower', 5500.0 / 3000.0)

        bp.add_fits_attribute('Chunk.time.exposure', 'EXPTIME')
        bp.add_fits_attribute('Chunk.time.resolution', 'EXPTIME')

    def get_telescope_name(self, ext):
        return 'DAO Skycam'

    def get_time_axis_val(self, ext):
        return ac.get_datetime(self._headers[ext].get('CLOCKVAL'))


class Imaging(Telescope):

    def __init__(self, headers):
        super().__init__(headers)

    def configure_axes(self, bp):
        bp.configure_position_axes((1, 2))
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)

    def accumulate_bp(self, bp):
        super().accumulate_bp(bp)
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

    def __init__(self, headers):
        super().__init__(headers)

    def accumulate_bp(self, bp):
        super().accumulate_bp(bp)
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


class ProcessedSpectrum(Telescope):

    def __init__(self, headers):
        super().__init__(headers)

    def configure_axes(self, bp):
        bp.configure_position_axes((2, 3))
        bp.configure_time_axis(4)
        bp.configure_energy_axis(1)
        bp.configure_observable_axis(5)

    def accumulate_bp(self, bp):
        super().accumulate_bp(bp)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        if DAOName.is_derived(self.uri):
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

    def __init__(self, headers):
        super().__init__(headers)


class Dao12MetreProcessedImage(Dao12MetreImage, ProcessedImage):

    def __init__(self, headers):
        super().__init__(headers)


class Dao12MetreSpectrum(Dao12Metre):

    def __init__(self, headers):
        super().__init__(headers)

    def _get_dispaxis(self, ext):
        return self._headers[ext].get('DISPAXIS', 2)

    def get_energy_axis_function_naxis(self, ext):
        dispaxis = self._get_dispaxis(ext)
        return self._headers[ext].get(f'NAXIS{dispaxis}')


class Dao12MetreProcessedSpectrum(Dao12MetreSpectrum, ProcessedSpectrum):

    def __init__(self, headers):
        super().__init__(headers)

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

    def __init__(self, headers):
        super().__init__(headers)


class Dao18MetreProcessedImage(Dao18MetreImage, ProcessedImage):

    def __init__(self, headers):
        super().__init__(headers)


class Dao18MetreSpectrum(Dao18Metre):

    def __init__(self, headers):
        super().__init__(headers)

    def _get_dispaxis(self, ext):
        return self._headers[ext].get('DISPAXIS', 1)

    def get_energy_axis_function_naxis(self, ext):
        dispaxis = self._get_dispaxis()
        return self._headers[ext].get(f'NAXIS{dispaxis}')


class Dao18MetreProcessedSpectrum(Dao18MetreSpectrum, ProcessedSpectrum):

    def __init__(self, headers):
        super().__init__(headers)

    def get_energy_resolving_power(self, ext):
        obs_type = self._headers[ext].get('OBS_TYPE')
        if obs_type in ['comparison', 'object']:
            numerator = self._headers[ext].get('CRVAL1')
            denominator = self._headers[ext].get('CDELT1')
        else:
            numerator = self._headers[ext].get('WAVELENG')
            denominator = self.get_energy_axis_function_delta(ext)
        return numerator / (2.5 * denominator)


# globally accessible pointer, placed here so they survive the
# importlib.import_module call from fits2caom2

# latest instances of Telescope
current = {}
# single instance of DefiningMetadataFinder
# defining_metadata_finder = None


def factory(uri, defining_metadata):
    # at this point, decompose the classes based on what can be determined
    # from the file name only
    ignore_scheme, ignore_path, f_name = mc.decompose_uri(uri)
    if f_name.startswith('a'):
        result = SkyCam()
    else:
        # defining_metadata = defining_metadata_finder.get(uri)
        f_id = DAOName.remove_extensions(f_name)
        if defining_metadata.data_product_type == DataProductType.IMAGE:
            if (
                re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122') or
                    f_name.startswith('dao_p122')
                ):
                    result = Dao12MetreProcessedImage()
                else:
                    result = Dao18MetreProcessedImage()
            elif (
                f_name.startswith('dao_c122') or
                f_name.startswith('dao_r122') or
                f_name.startswith('dao_p122')
            ):
                result = Dao12MetreImage()
            else:
                result = Dao18MetreImage()
        else:
            if (
                re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122') or
                    f_name.startswith('dao_p122')
                ):
                    result = Dao12MetreProcessedSpectrum()
                else:
                    result = Dao18MetreProcessedSpectrum()
            elif (
                f_name.startswith('dao_c122') or
                f_name.startswith('dao_r122') or
                f_name.startswith('dao_p122')
            ):
                result = Dao12MetreSpectrum()
            else:
                result = Dao18MetreSpectrum()

    result.uri = uri
#     global current
#     current[uri] = result
#
#
# def get_current(uri):
#     return current.get(uri)


def factory_client(uri, headers):
    # at this point, decompose the classes based on what can be determined
    # from the file name only
    ignore_scheme, ignore_path, f_name = mc.decompose_uri(uri)
    if f_name.startswith('a'):
        result = SkyCam(heaaders)
    else:
        # defining_metadata = defining_metadata_finder.get(uri)
        f_id = DAOName.remove_extensions(f_name)
        data_product_type = DefiningMetadataFinder.get_data_product_type(
            headers
        )
        if data_product_type == DataProductType.IMAGE:
            if (
                    re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                    re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                        f_name.startswith('dao_c122') or
                        f_name.startswith('dao_r122') or
                        f_name.startswith('dao_p122')
                ):
                    result = Dao12MetreProcessedImage(headers)
                else:
                    result = Dao18MetreProcessedImage(headers)
            elif (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122') or
                    f_name.startswith('dao_p122')
            ):
                result = Dao12MetreImage(headers)
            else:
                result = Dao18MetreImage(headers)
        else:
            if (
                    re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', f_id) or
                    re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', f_id)
            ):
                if (
                        f_name.startswith('dao_c122') or
                        f_name.startswith('dao_r122') or
                        f_name.startswith('dao_p122')
                ):
                    result = Dao12MetreProcessedSpectrum(headers)
                else:
                    result = Dao18MetreProcessedSpectrum(headers)
            elif (
                    f_name.startswith('dao_c122') or
                    f_name.startswith('dao_r122') or
                    f_name.startswith('dao_p122')
            ):
                result = Dao12MetreSpectrum(headers)
            else:
                result = Dao18MetreSpectrum(headers)

    result.uri = uri
    return result


# def get_current(uri):
#     return current.get(uri)
