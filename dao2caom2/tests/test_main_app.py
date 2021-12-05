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

import warnings
from astropy.utils.exceptions import AstropyUserWarning

from cadcdata import FileInfo
from caom2 import DataProductType
from dao2caom2 import main_app, APPLICATION, COLLECTION, DAOName
from dao2caom2 import metadata, telescopes, PRODUCT_COLLECTION
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

from mock import patch, Mock

import os
import sys

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')


def pytest_generate_tests(metafunc):
    files = []
    if os.path.exists(TEST_DATA_DIR):
        files = [
            os.path.join(TEST_DATA_DIR, name)
            for name in os.listdir(TEST_DATA_DIR)
            if name.endswith('header')
        ]
    metafunc.parametrize('test_name', files)


@patch('caom2utils.data_util.HeaderReader')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('dao2caom2.metadata.DefiningMetadataFinder.check_local')
@patch('caom2utils.data_util.StorageClientWrapper')
def test_main_app(
    data_client_mock,
    local_headers_mock,
    util_headers_mock,
    header_client_mock,
    test_name,
):
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    header_client_mock.get_head.side_effect = ac.make_headers_from_file
    local_headers_mock.side_effect = _local_headers
    util_headers_mock.side_effect = ac.make_headers_from_file
    config = mc.Config()
    config.use_local_files = True
    config.data_sources = [TEST_DATA_DIR]
    config.task_types = [mc.TaskType.SCRAPE]
    dmf = metadata.DefiningMetadataFinder(Mock(), config)
    telescopes.defining_metadata_finder = dmf

    data_client_mock.return_value.info.side_effect = _get_file_info
    basename = os.path.basename(test_name)
    dao_name = DAOName(basename.replace('.header', '.gz'))
    collection = (
        PRODUCT_COLLECTION if DAOName.is_processed(test_name) else COLLECTION
    )

    obs_path = f'{TEST_DATA_DIR}/{dao_name.file_id}.expected.xml'
    output_file = f'{TEST_DATA_DIR}/{dao_name.file_id}.actual.xml'

    sys.argv = (
        f'{APPLICATION} --no_validate --local {test_name} '
        f'--observation {collection} {dao_name.obs_id} -o {output_file} '
        f'--plugin {PLUGIN} --module {PLUGIN} --lineage '
        f'{dao_name.lineage}'
    ).split()
    print(sys.argv)
    try:
        main_app.to_caom2_with_client(header_client_mock)
    except Exception as e:
        assert False, e

    compare_result = mc.compare_observations(output_file, obs_path)
    if compare_result is not None:
        raise AssertionError(compare_result)
    # assert False  # cause I want to see logging messages


def _get_file_info(file_id):
    return FileInfo(id=file_id, file_type='application/fits')


def _local_headers(uri):
    ign1, ign2, f_name = mc.decompose_uri(uri)
    headers = ac.make_headers_from_file(
        f'{TEST_DATA_DIR}/{f_name.replace(".gz", ".header")}'
    )
    temp = headers[0].get('OBSMODE')
    dpt = (
        DataProductType.SPECTRUM if '-slit' in temp else DataProductType.IMAGE
    )
    return metadata.DefiningMetadata(dpt, uri)
