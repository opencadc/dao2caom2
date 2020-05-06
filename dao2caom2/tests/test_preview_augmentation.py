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

import glob
import os

from mock import patch, Mock

from caom2pipe import manage_composable as mc
from dao2caom2 import preview_augmentation, dao_name
import test_main_app

TEST_FILES_DIR = '/test_files'
REJECTED_FILE = os.path.join(test_main_app.TEST_DATA_DIR, 'rejected.yml')


@patch('caom2pipe.manage_composable.data_put')
def test_visit(ad_put_mock):

    # this should result in three new artifacts being added to every plane:
    # one for a thumbnail and two for previews (one zoom)

    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    cadc_client_mock = Mock()

    test_files = {
        # processed spectrum
        'visit_obs_start_v.xml':
            ['dao_c122_2007_000882_v.fits', 'dao_c122_2007_000882.fits.gz'],
        # processed image
        'visit_obs_start_a.xml':
            ['dao_c182_2016_004034_a.fits', 'dao_c182_2016_004034.fits.gz'],
        'dao_r182_1989_000369.xml': ['dao_r182_1989_000369.fits.gz'],
        'dao_r122_1989_003111.xml': ['dao_r122_1989_003111.fits.gz'],
        'dao_c122_2017_011124.xml': ['dao_c122_2017_011124.fits.gz'],
        'dao_c182_2017_010870.xml': ['dao_c182_2017_010870.fits.gz'],
        'dao_c182_2017_019322.xml': ['dao_c182_2017_019322.fits.gz'],
        'dao_c182_2017_016292.xml': ['dao_c182_2017_016292.fits.gz'],
        'visit_obs_start_e.xml':
            ['dao_c122_2007_000881.fits.gz', 'dao_c122_2007_000881_e.fits']
    }

    test_checksums = {
        'ad:DAO/dao_c122_2007_000882_v_prev.jpg':
            'md5:72ed80fafb213084299d2c037f75de9e',
        'ad:DAO/dao_c122_2007_000882_v_prev_256.jpg':
            'md5:b33f81f6b182d6172e3a5054afbb2c6f',
        'ad:DAO/dao_c122_2007_000882_prev.jpg':
            'md5:2112cc01d7a786226b402baa965b5d4a',
        'ad:DAO/dao_c122_2007_000882_prev_256.jpg':
            'md5:e40445f4b615b79433579971c1d76cbd',
        'ad:DAO/dao_c182_2016_004034_a_prev.jpg':
            'md5:f8c3cbad294d59fa1b971039c3c201f0',
        'ad:DAO/dao_c182_2016_004034_a_prev_256.jpg':
            'md5:23e0f270b575ac85f1749f12080788f9',
        'ad:DAO/dao_c182_2016_004034_prev.jpg':
            'md5:cef242b052395a8407928e4afbf9f542',
        'ad:DAO/dao_c182_2016_004034_prev_256.jpg':
            'md5:0e69e3359137b14578f2c7c742c47456',
        'ad:DAO/dao_r182_1989_000369_prev.jpg':
            'md5:6aecde907aeb40e377e95c62d82730e3',
        'ad:DAO/dao_r182_1989_000369_prev_256.jpg':
            'md5:704c5b8ed3e5e7551f8761c396269044',
        'ad:DAO/dao_r122_1989_003111_prev.jpg':
            'md5:f1ddc6905cfbe56c2f1adff64a53114e',
        'ad:DAO/dao_r122_1989_003111_prev_256.jpg':
            'md5:58c66377dae7a612b1ef3a4c56ecef1b',
        'ad:DAO/dao_c122_2017_011124_prev.jpg':
            'md5:a398bc712f003b085e03d270a7ba31eb',
        'ad:DAO/dao_c122_2017_011124_prev_256.jpg':
            'md5:9da055e9754300da6e5c310880d29892',
        'ad:DAO/dao_c182_2017_010870_prev.jpg':
            'md5:16ea25f015601180cdc0bff30a69ad9b',
        'ad:DAO/dao_c182_2017_010870_prev_256.jpg':
            'md5:8b4c0009485a8e282657f98b2a573d98',
        'ad:DAO/dao_c182_2017_019322_prev.jpg':
            'md5:38c44f703df7a548b40e4a93bfacd905',
        'ad:DAO/dao_c182_2017_019322_prev_256.jpg':
            'md5:3069a4cff264f4b7f1157b1c0c2a22da',
        'ad:DAO/dao_c182_2017_016292_prev.jpg':
            'md5:88266cc671db81ad3804ed9f00662f51',
        'ad:DAO/dao_c182_2017_016292_prev_256.jpg':
            'md5:41c5e368839076d5efbcdad09be094af',
        'ad:DAO/dao_c122_2007_000881_prev.jpg':
            'md5:16ad2bc4747fb2360813a4cc34f0ad84',
        'ad:DAO/dao_c122_2007_000881_prev_256.jpg':
            'md5:33ac01e8211e535552495b496dc286fb',
        'ad:DAO/dao_c122_2007_000881_e_prev.jpg':
            'md5:854bf08b08ec703fcca4fd57d1890360',
        'ad:DAO/dao_c122_2007_000881_e_prev_256.jpg':
            'md5:43764c2acc3dd93d3c4b7086a4eee58c',
    }

    kwargs = {'working_directory': TEST_FILES_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'observable': test_observable}

    for entry in glob.glob(f'{TEST_FILES_DIR}/*.jpg'):
        os.unlink(entry)

    for key, value in test_files.items():
        obs = mc.read_obs_from_file(
            f'{test_main_app.TEST_DATA_DIR}/previews/{key}')
        for f_name in value:
            kwargs['science_file'] = f_name

            try:
                ignore = preview_augmentation.visit(obs, **kwargs)
            except Exception as e:
                assert False, f'{str(e)}'

            test_name = dao_name.DAOName(file_name=f_name)
            f_name_list = [test_name.prev_uri, test_name.thumb_uri]
            for p in f_name_list:
                artifact = obs.planes[test_name.product_id].artifacts[p]
                assert artifact.content_checksum.uri == test_checksums[p], \
                    f'wrong checksum {p} {artifact.content_checksum} ' \
                    f'{test_checksums[p]}'
    # assert False
