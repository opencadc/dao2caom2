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

"""
DB - 06-05-20
Thumbnails and Previews are proprietary for science datasets.

"""
import logging
import os

import numpy as np

from astropy.io import fits
from matplotlib import pylab

from caom2 import ReleaseType, ProductType
from caom2pipe import manage_composable as mc
from dao2caom2 import dao_name as dn


class DAOPreview(mc.PreviewVisitor):

    def __init__(self, **kwargs):
        super(DAOPreview, self).__init__(
            dn.COLLECTION, ReleaseType.DATA, **kwargs)

    def generate_plots(self, obs_id):
        count = 0
        storage_name = dn.DAOName(file_name=self._science_file)
        science_fqn = os.path.join(self._working_dir, storage_name.file_name)
        preview = storage_name.prev
        preview_fqn = os.path.join(self._working_dir, preview)
        thumb = storage_name.thumb
        thumb_fqn = os.path.join(self._working_dir, thumb)
        hdu_list = fits.open(science_fqn)
        header = hdu_list[0].header

        if ('e' in storage_name.file_name or 'p' in storage_name.file_name or
                'v' in storage_name.file_name):
            count += self._do_cal_processed(
                hdu_list, header, science_fqn, storage_name, preview_fqn,
                thumb_fqn, obs_id)
        else:
            count += self._do_sci(hdu_list, header, storage_name, science_fqn,
                                  preview_fqn, thumb_fqn)

        hdu_list.close()

        if count == 2:
            self.add_preview(storage_name.thumb_uri, storage_name.thumb,
                             ProductType.THUMBNAIL)
            self.add_preview(storage_name.prev_uri, storage_name.prev,
                             ProductType.PREVIEW)
            self.add_to_delete(thumb_fqn)
            self.add_to_delete(preview_fqn)
        return count

    def _do_cal_processed(self, hdu_list, header, science_fqn, storage_name,
                          preview_fqn, thumb_fqn, obs_id):
        logging.debug(f'Do calibration preview augmentation with {science_fqn}')
        count = 0
        object_type = header.get('OBJECT')
        if object_type is None:
            self._logger.warning(f'Stopping preview generation. No object '
                                 f'type for {obs_id}.')
        else:
            crval1 = header.get('CRVAL1')
            crpix1 = header.get('CRPIX1')
            cd1_1 = header.get('CD1_1')
            naxis1 = header.get('NAXIS1')
            logging.info(f'Object: {object_type}')

            # if daoplate, daoPlate = False if 'v' in science_fqn:
            if 'v' in science_fqn or 'e' in science_fqn:
                flux = hdu_list[0].data
            else:
                flux = hdu_list[0].data[0]

            hdu_list.close()
            wl = []
            for i in range(0, naxis1):
                wl.append(crval1 + cd1_1 * (float(i) - crpix1 - 1.0))

            wln = np.array(wl)
            self._write_files_to_disk(wln, flux, 'Wavelength ($\AA$)',
                                      f'{storage_name.file_id}: {object_type}',
                                      thumb_fqn, preview_fqn)
            count = 2
        return count

    def _do_sci(self, hdu_list, header, storage_name, science_fqn, preview_fqn,
                thumb_fqn):
        logging.debug(f'Do science preview augmentation with {science_fqn}')
        count = 0
        detector = header.get('DETECTOR')
        instrument = header.get('INSTRUME')
        if detector in ['SITe-4', 'UBC-1', 'SITe-2', 'SITe-5', 'E2V-1',
                        'E2V-4']:
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
                naxis1 = mc.to_int(header.get(axis))
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
                self._write_files_to_disk(
                    wln, flux, 'Pixel',
                    f'{storage_name.file_id}: {object_type}', thumb_fqn,
                    preview_fqn)
                count = 2
        return count

    def _write_files_to_disk(self, wln, flux, x_label, title, thumb_fqn,
                             preview_fqn):
        pylab.clf()
        pylab.grid(True)
        pylab.plot(wln, flux, color='k')
        pylab.xlabel(x_label, color='k')
        pylab.ylabel('Intensity', color='k')
        pylab.xlim(wln.min(), wln.max())
        pylab.ylim(0, flux.max())
        pylab.title(title, color='k', fontweight='bold')
        temp_fn = 'temp.jpg'
        pylab.savefig(temp_fn, format='jpg')
        mc.exec_cmd(f'convert -resize 256x256 {temp_fn} {thumb_fqn}')
        mc.exec_cmd(f'convert -resize 1024x1024 {temp_fn} {preview_fqn}')
        self.add_to_delete(f'{self._working_dir}/{temp_fn}')


def visit(observation, **kwargs):
    previewer = DAOPreview(**kwargs)
    dao_name = dn.DAOName(file_name=previewer.science_file)
    return previewer.visit(observation, dao_name)
