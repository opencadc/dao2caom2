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

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.visualization import ZScaleInterval, SqrtStretch, ImageNormalize
from matplotlib import pylab
from PIL import Image

from caom2 import ReleaseType, ProductType
from caom2pipe import manage_composable as mc


class DAOPreview(mc.PreviewVisitor):
    def __init__(self, mime_type, **kwargs):
        super(DAOPreview, self).__init__(ReleaseType.DATA, mime_type, **kwargs)
        self._ext = 0

    def generate_plots(self, obs_id):
        count = 0
        self._logger.info(f'Building preview and thumbnail with {self._science_fqn}')
        self._hdu_list = fits.open(self._science_fqn)
        header = self._hdu_list[self._ext].header

        if (
            'e' in self._storage_name.file_name
            or 'p' in self._storage_name.file_name
            or 'v' in self._storage_name.file_name
        ):
            count += self._do_cal_processed(header, obs_id)
        elif self._storage_name.file_name.startswith('a'):
            count += self._do_skycam()
        else:
            count += self._do_sci(header)

        if count == 2:
            self.add_preview(
                self._storage_name.thumb_uri,
                self._storage_name.thumb,
                ProductType.THUMBNAIL,
            )
            self.add_preview(
                self._storage_name.prev_uri,
                self._storage_name.prev,
                ProductType.PREVIEW,
            )
            self.add_to_delete(self._thumb_fqn)
            self.add_to_delete(self._preview_fqn)
        return count

    def _do_cal_processed(self, header, obs_id):
        self._logger.debug(f'Do calibration preview augmentation with {self._science_fqn}')
        count = 0
        object_type = header.get('OBJECT')
        if object_type is None:
            self._logger.warning(f'Stopping preview generation. No object type for {obs_id}.')
        else:
            crval1 = header.get('CRVAL1')
            crpix1 = header.get('CRPIX1')
            cd1_1 = header.get('CD1_1')
            naxis1 = header.get('NAXIS1')
            self._logger.info(f'Object: {object_type}')

            if 'v' in self._science_fqn or 'e' in self._science_fqn:
                flux = self._hdu_list[self._ext].data
            else:
                flux = self._hdu_list[self._ext].data[0]

            wl = []
            for i in range(0, naxis1):
                wl.append(crval1 + cd1_1 * (float(i) - crpix1 - 1.0))

            wln = np.array(wl)
            self._write_files_to_disk(
                wln,
                flux,
                'Wavelength ($\AA$)',
                f'{self._storage_name.file_id}: {object_type}',
            )
            count = 2
        return count

    def _do_sci(self, header):
        self._logger.debug(f'Do science preview augmentation with {self._science_fqn}')
        count = 0
        detector = header.get('DETECTOR')
        if detector.upper() == 'RETICON':
            # unprocessed RETICON spectrum
            object_type = header.get('OBJECT')
            if object_type is not None:
                self._logger.info(f'Object: {object_type}')
                naxis1 = header.get('NAXIS1')
                signal = self._hdu_list[self._ext].data[0]
                baseline = self._hdu_list[self._ext].data[1]
                flux = np.subtract(signal, baseline)
                wl = []
                for i in range(0, naxis1):
                    wl.append(i + 1)
                wln = np.array(wl)
                self._write_files_to_disk(wln, flux, 'Pixel', f'{self._storage_name.file_id}: {object_type}')
                count = 2
        else:
            instrument = header.get('INSTRUME')
            if 'Imager' in instrument:
                preview_cmd = f'convert -resize 1024x1024 -normalize -negate {self._science_fqn} {self._preview_fqn}'
                thumbnail_cmd = f'convert -resize 256x256 -normalize -negate {self._science_fqn} {self._thumb_fqn}'
            else:
                # unprocessed CCD data
                dispaxis = mc.to_int(header.get('DISPAXIS'))
                if dispaxis == 2:
                    axis = 'NAXIS2'
                    naxis1 = mc.to_int(header.get(axis))
                    xc = naxis1 / 2
                    xs = 512
                    xoffset = xc - xs / 2
                    rotate = '90.0'
                    geometry = '256x' + str(xs) + '+1+' + str(xoffset)
                    resize1 = 'x1024'
                    resize2 = 'x256'
                else:
                    axis = 'NAXIS1'
                    naxis1 = mc.to_int(header.get(axis))
                    xc = naxis1 / 2
                    xs = 512
                    xoffset = xc - xs / 2
                    rotate = '0.0'
                    geometry = str(xs) + 'x256+' + str(xoffset) + '+1'
                    resize1 = '1024x1024'
                    resize2 = '256'

                preview_cmd = (
                    f'convert -resize {resize1} -rotate {rotate} -normalize -negate {self._science_fqn} '
                    f'{self._preview_fqn}'
                )
                thumbnail_cmd = (
                    f'convert -crop {geometry} -resize {resize2} -rotate {rotate} -normalize -negate '
                    f'{self._science_fqn} {self._thumb_fqn}'
                )
            self._logger.debug(f'Preview Generation: {preview_cmd}')
            mc.exec_cmd(preview_cmd)
            self._logger.debug(f'Thumbnail Generation: {thumbnail_cmd}')
            mc.exec_cmd(thumbnail_cmd)
            count = 2
        return count

    def _do_skycam(self):
        self._hdulist = fits.open(self._science_fqn)
        image_data = self._hdulist[self._ext].data
        norm = ImageNormalize(
            image_data[330:950, 215:750],
            interval=ZScaleInterval(),
            stretch=SqrtStretch(),
        )
        plt.imshow(image_data, cmap='gray', norm=norm)
        plt.gca().invert_yaxis()
        plt.axis('off')
        with np.errstate(invalid='ignore'):
            plt.savefig(
                self._preview_fqn,
                dpi=100,
                bbox_inches='tight',
                pad_inches=0,
                format='png',
            )
        count = 1
        count += self._gen_thumbnail()
        return count

    def _write_files_to_disk(self, wln, flux, x_label, title):
        pylab.clf()
        pylab.grid(True)
        pylab.plot(wln, flux, color='k')
        pylab.xlabel(x_label, color='k')
        pylab.ylabel('Intensity', color='k')
        pylab.xlim(wln.min(), wln.max())
        pylab.ylim(0, flux.max())
        pylab.title(title, color='k', fontweight='bold')
        temp_fn = 'temp.png'
        pylab.savefig(temp_fn, format='png')
        img = Image.open(temp_fn)
        img.thumbnail((1024, 1024))
        img.save(self._preview_fqn)
        img.thumbnail((256, 256))
        img.save(self._thumb_fqn)
        self.add_to_delete(f'{self._working_dir}/{temp_fn}')


def visit(observation, **kwargs):
    previewer = DAOPreview(mime_type='image/png', **kwargs)
    return previewer.visit(observation)
