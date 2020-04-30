# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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

import logging
import os

import numpy as np

from astropy.io import fits
from matplotlib import pylab

from caom2 import Observation, ReleaseType, ProductType
from caom2pipe import manage_composable as mc
from dao2caom2 import dao_name as dn


MIME_TYPE = 'image/jpeg'


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    cadc_client = kwargs.get('cadc_client')
    if cadc_client is None:
        logging.warning(
            'Visitor needs a cadc_client parameter to store previews.')
    stream = kwargs.get('stream')
    if stream is None:
        raise mc.CadcException('Visitor needs a stream parameter.')
    observable = kwargs.get('observable')
    if observable is None:
        raise mc.CadcException('Visitor needs a observable parameter.')
    science_file = kwargs.get('science_file')
    if science_file is None:
        raise mc.CadcException('Visitor needs a science_file parameter.')

    count = 0
    storage_name = dn.DAOName(file_name=science_file)
    for plane in observation.planes.values():
        delete_list = []
        for artifact in plane.artifacts.values():
            if artifact.uri.endswith(science_file):
                count += _do_prev(plane, working_dir, cadc_client,
                                  storage_name, stream, observable)
            if artifact.uri.endswith('.jpg'):
                delete_list.append(artifact.uri)

        for uri in delete_list:
            plane.artifacts.pop(uri)

    logging.info('Completed preview augmentation for {}.'.format(
        observation.observation_id))
    return {'artifacts': count}


def _do_prev(plane, working_dir, cadc_client, storage_name, stream,
             observable):
    count = 0
    science_fqn = os.path.join(working_dir, storage_name.file_name)
    preview = storage_name.prev
    preview_fqn = os.path.join(working_dir, preview)
    thumb = storage_name.thumb
    thumb_fqn = os.path.join(working_dir, thumb)
    hdu_list = fits.open(science_fqn)
    header = hdu_list[0].header

    if 'v' in storage_name.file_name:
        count += _do_cal_processed(
            hdu_list, header, science_fqn, storage_name, preview_fqn,
            thumb_fqn)
    else:
        count += _do_sci(hdu_list, header, storage_name, science_fqn,
                         preview_fqn, thumb_fqn)

    _augment(plane, storage_name.prev_uri, preview_fqn, ProductType.PREVIEW)
    _augment(plane, storage_name.thumb_uri, thumb_fqn, ProductType.THUMBNAIL)
    _store_smalls(cadc_client, working_dir, stream, preview, thumb,
                  observable.metrics)

    return count


def _do_cal_processed(hdu_list, header, science_fqn, storage_name,
                      preview_fqn, thumb_fqn):
    logging.debug(f'Do calibration preview augmentation with {science_fqn}')
    count = 0
    object_type = header.get('OBJECT')
    if object_type is not None:
        crval1 = header.get('CRVAL1')
        crpix1 = header.get('CRPIX1')
        cd1_1 = header.get('CD1_1')
        naxis1 = header.get('NAXIS1')
        logging.info(f'Object: {object_type}')

        # if daoplate, daoPlate = False if 'v' in science_fqn:
        if 'v' in science_fqn:
            flux = hdu_list[0].data
        else:
            flux = hdu_list[0].data[0]

        hdu_list.close()
        wl = []
        for i in range(0, naxis1):
            wl.append(crval1 + cd1_1 * (float(i) - crpix1 - 1.0))

        wln = np.array(wl)
        _do_write_files(wln, flux, 'Wavelength ($\AA$)',
                        f'{storage_name.file_id}: {object_type}', thumb_fqn,
                        preview_fqn)
        count = 2
    return count


def _do_write_files(wln, flux, x_label, title, thumb_fqn, preview_fqn):
    pylab.clf()
    pylab.grid(True)
    pylab.plot(wln, flux, color='k' )
    pylab.xlabel(x_label, color='k')
    pylab.ylabel('Intensity', color='k')
    pylab.xlim(wln.min(), wln.max())
    pylab.ylim(0, flux.max())
    pylab.title(title, color='k', fontweight='bold')
    temp_fn = 'temp.jpg'
    pylab.savefig(temp_fn, format='jpg')
    mc.exec_cmd(f'convert -resize 256x256 {temp_fn} {thumb_fqn}')
    mc.exec_cmd(f'convert -resize 1024x1024 {temp_fn} {preview_fqn}')


def _do_sci(hdu_list, header, storage_name, science_fqn, preview_fqn, thumb_fqn):
    logging.debug(f'Do science preview augmentation with {science_fqn}')
    count = 0
    detector = header.get('DETECTOR')
    instrument = header.get('INSTRUME')
    if detector in ['SITe-4', 'UBC-1', 'SITe-2', 'SITe-5', 'E2V-1', 'E2V-4']:
        # unprocessed CCD data
        if detector == 'SITe-4':
            axis = 'NAXIS2'
            naxis1 = mc.to_int(header.get(axis))
            xc = naxis1/2
            xs = 512
            xoffset = xc - xs/2
            rotate = '90.0'
            geometry = '256x' + str(xs) + '+1+' + str(xoffset)
            resize1 = 'x1024'
            resize2 = 'x256'
        else:
            axis = 'NAXIS1'
            naxis1 =mc.to_int(header.get(axis))
            xc = naxis1/2
            xs = 512
            xoffset = xc - xs/2
            rotate = '0.0'
            geometry = str(xs) + 'x256+' + str(xoffset) + '+1'
            resize1 = '1024x1024'
            resize2 = '256'

        if 'Imager' in instrument:
            mc.exec_cmd(f'convert -resize 1024x1024 '
                        f'-normalize -negate {science_fqn} {preview_fqn}')
            mc.exec_cmd(f'convert -resize 256x256 -normalize -negate '
                        f'{science_fqn} {thumb_fqn}')
        else:
            mc.exec_cmd(f'convert -resize {resize1} -rotate {rotate} '
                        f'-normalize -negate {science_fqn} {preview_fqn}')
            mc.exec_cmd(f'convert -crop {geometry} -resize {resize2} '
                        f'-rotate {rotate} -normalize '
                        f'-negate {science_fqn} {thumb_fqn}')
        count = 2
    else:
        # unprocessed RETICON spectrum
        object_type = header.get('OBJECT')
        if object_type is not None:
            naxis1 = header.get('NAXIS1')
            logging.info(f'Object: {object_type}')

            signal = hdu_list[0].data[0]
            baseline = hdu_list[0].data[1]
            flux = np.subtract(signal, baseline)
            wl = []
            for i in range(0, naxis1):
                wl.append(i + 1)
            wln = np.array(wl)
            _do_write_files(wln, flux, 'Pixel',
                            f'{storage_name.file_id}: {object_type}',
                            thumb_fqn, preview_fqn)
            count = 2
    hdu_list.close()
    return count


def _augment(plane, uri, fqn, product_type):
    temp = None
    if uri in plane.artifacts:
        temp = plane.artifacts[uri]
    plane.artifacts[uri] = mc.get_artifact_metadata(
        fqn, product_type, ReleaseType.DATA, uri, temp)


def _store_smalls(cadc_client, working_directory, stream, preview_fname,
                  thumb_fname, metrics):
    if cadc_client is not None:
        mc.data_put(cadc_client, working_directory, preview_fname,
                    dn.COLLECTION, stream, mime_type=MIME_TYPE,
                    metrics=metrics)
        mc.data_put(cadc_client, working_directory, thumb_fname, dn.COLLECTION,
                    stream, mime_type=MIME_TYPE, metrics=metrics)

        # if there's no client, leave the files behind - this is consistent
        # with TaskType.SCRAPE expectations
        mc.delete_list_of_files([f'{working_directory}/{preview_fname}',
                                 f'{working_directory}/{thumb_fname}',
                                 f'{working_directory}/temp.jpg'])
