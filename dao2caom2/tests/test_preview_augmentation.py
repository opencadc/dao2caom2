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

from caom2pipe import manage_composable as mc
from dao2caom2 import preview_augmentation, dao_name
import test_fits2caom2_augmentation

TEST_FILES_DIR = '/test_files'
REJECTED_FILE = os.path.join(
    test_fits2caom2_augmentation.TEST_DATA_DIR, 'rejected.yml'
)


def test_visit(test_config, tmp_path):

    # this should result in three new artifacts being added to every plane:
    # one for a thumbnail and two for previews (one zoom)

    test_config.rejected_fqn = f'{tmp_path}/{test_config.rejected_file_name}'
    test_observable = mc.Observable(test_config)

    test_files = {
        # processed spectrum
        'visit_obs_start_v.xml': [
            'dao_c122_2007_000882_v.fits',
            'dao_c122_2007_000882.fits.gz',
        ],
        # processed image
        'visit_obs_start_a.xml': [
            'dao_c182_2016_004034_a.fits',
            'dao_c182_2016_004034.fits.gz',
        ],
        'dao_r182_1989_000369.xml': ['dao_r182_1989_000369.fits.gz'],
        'dao_r122_1989_003111.xml': ['dao_r122_1989_003111.fits.gz'],
        'dao_c122_2017_011124.xml': ['dao_c122_2017_011124.fits.gz'],
        'dao_c182_2017_010870.xml': ['dao_c182_2017_010870.fits.gz'],
        'dao_c182_2017_019322.xml': ['dao_c182_2017_019322.fits.gz'],
        'dao_c182_2017_016292.xml': ['dao_c182_2017_016292.fits.gz'],
        'visit_obs_start_e.xml': [
            'dao_c122_2007_000881.fits.gz',
            'dao_c122_2007_000881_e.fits',
        ],
        'sky_cam_start.xml': ['a2020_06_17_07_00_01.fits'],
    }

    test_checksums = {
        'cadc:DAO/dao_c122_2007_000881_e_256.png': 'md5:e68a3c389bbdf53cde0cca073c7a37c1',
        'cadc:DAO/dao_c122_2007_000881_e_1024.png': 'md5:48065b2feb0fda3d7561c02b2890efaa',
        'cadc:DAO/dao_c122_2007_000882_256.png': 'md5:1268e5b79e463ca75580d7e538d88a2a',
        'cadc:DAO/dao_c122_2007_000882_1024.png': 'md5:5694ee03486899d3b14b404caab22505',
        'cadc:DAO/dao_c122_2017_011124_256.png': 'md5:9793efcdeec043be8f0c77c7bf875cca',
        'cadc:DAO/dao_c122_2017_011124_1024.png': 'md5:4d77005756fbed7fc3ae525591cc07db',
        'cadc:DAO/dao_c182_2016_004034_256.png': 'md5:5c6b6816f8c0a54a0ff041b29247bdd5',
        'cadc:DAO/dao_c182_2016_004034_1024.png': 'md5:455bddac53488acca46758fbc1d2b096',
        'cadc:DAO/dao_c182_2017_010870_256.png': 'md5:106a6e0837fd00e971a37462266c2f66',
        'cadc:DAO/dao_c182_2017_010870_1024.png': 'md5:bc92d495aa0e56eb4ea35f409932391e',
        'cadc:DAO/dao_c182_2017_016292_256.png': 'md5:ef6bcc6501db0f64dfc7741f5ab37870',
        'cadc:DAO/dao_c182_2017_016292_1024.png': 'md5:6c425b3606dd598be33b61dfef51d132',
        'cadc:DAO/dao_c182_2017_019322_256.png': 'md5:dfca85e675f29a33ecdbabc91ed274bf',
        'cadc:DAO/dao_c182_2017_019322_1024.png': 'md5:2818ddf6acb0f36655fbdfc9e19098e7',
        'cadc:DAO/dao_r122_1989_003111_256.png': 'md5:7cb2344a84d195228d7fb866a83d031a',
        'cadc:DAO/dao_r122_1989_003111_1024.png': 'md5:a864bc14955bc1f8524e8829c3ee2ee9',
        'cadc:DAO/dao_r182_1989_000369_256.png': 'md5:5ef8f53784376b31a12dccb8671aca5e',
        'cadc:DAO/dao_r182_1989_000369_1024.png': 'md5:2e14aefda2a74b5575f6a4d8d1b1b3cb',
        'cadc:DAO/dao_c122_2007_000882_v_1024.png': 'md5:37cfb96bce27da93ede75eadabad7b8f',
        'cadc:DAO/a2020_06_17_07_00_01_1024.png': 'md5:37cfb96bce27da93ede75eadabad7b8f',
    }

    kwargs = {
        'working_directory': TEST_FILES_DIR,
        'cadc_client': None,
        'stream': 'stream',
        'observable': test_observable,
    }

    for entry in glob.glob(f'{TEST_FILES_DIR}/*.png'):
        os.unlink(entry)

    for key, value in test_files.items():
        obs = mc.read_obs_from_file(
            f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/previews/{key}'
        )
        for f_name in value:
            test_name = dao_name.DAOName(f_name)
            kwargs['storage_name'] = test_name

            try:
                ignore = preview_augmentation.visit(obs, **kwargs)
                f_name_list = [test_name.prev_uri, test_name.thumb_uri]
                for p in f_name_list:
                    artifact = obs.planes[test_name.product_id].artifacts[p]
                    # assert artifact.content_checksum.uri ==
                    # test_checksums[p], \
                    #     f'wrong checksum {p} {artifact.content_checksum} ' \
                    #     f'{test_checksums[p]}'
            except Exception as e:
                assert False, f'key {key} value {value} f_name {f_name} {str(e)}'
    # assert False
