# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

import logging
import traceback

from dao2caom2 import main_app, APPLICATION, COLLECTION, DAOName
from dao2caom2 import metadata, telescopes
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

import os
import sys

from mock import patch, Mock
import test_main_app

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'processed'
LOOKUP = {
    'dao_c122_2007_000882': ['dao_c122_2007_000882_v.fits.gz'],
    'dao_c122_2007_000881': [
        'dao_c122_2007_000881.fits',
        'dao_c122_2007_000881_e.fits.gz',
    ],
    'dao_c182_2016_004034': [
        'dao_c182_2016_004034.fits',
        'dao_c182_2016_004034_a.fits.gz',
    ],
}


def pytest_generate_tests(metafunc):
    obs_id_list = LOOKUP.keys()
    metafunc.parametrize('test_name', obs_id_list)


@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2utils.data_util.StorageClientWrapper')
def test_multi_plane(data_client_mock, local_headers_mock, test_name):
    local_headers_mock.side_effect = ac.make_headers_from_file
    config = mc.Config()
    config.use_local_files = True
    config.data_sources = [f'{test_main_app.TEST_DATA_DIR}/processed']
    config.task_types = [mc.TaskType.SCRAPE]
    dmf = metadata.DefiningMetadataFinder(Mock(), config)
    telescopes.defining_metadata_finder = dmf

    dao_name = DAOName(file_name=f'{LOOKUP[test_name][0]}.fits')
    lineage = _get_lineage(dao_name.obs_id)
    expected_fqn = (
        f'{test_main_app.TEST_DATA_DIR}/{DIR_NAME}/'
        f'{dao_name.obs_id}.expected.xml'
    )
    actual_fqn = (
        f'{test_main_app.TEST_DATA_DIR}/{DIR_NAME}/'
        f'{dao_name.obs_id}.actual.xml'
    )

    local = _get_local(test_name)
    data_client_mock.return_value.info.side_effect = (
        test_main_app._get_file_info
    )

    if os.path.exists(actual_fqn):
        os.remove(actual_fqn)

    sys.argv = (
        f'{APPLICATION} --quiet --no_validate --local {local} '
        f'--observation {COLLECTION} {dao_name.obs_id} '
        f'--plugin {test_main_app.PLUGIN} --module {test_main_app.PLUGIN} '
        f'--out {actual_fqn} --lineage {lineage}'
    ).split()
    print(sys.argv)
    try:
        main_app.to_caom2()
    except Exception as e:
        logging.error(traceback.format_exc())

    compare_result = mc.compare_observations(actual_fqn, expected_fqn)
    if compare_result is not None:
        raise AssertionError(compare_result)
    # assert False  # cause I want to see logging messages


def _get_lineage(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        fits = mc.get_lineage(
            COLLECTION, mc.StorageName.remove_extensions(ii), ii
        )
        result = f'{result} {fits}'
    return result


def _get_local(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        result = f'{result} {test_main_app.TEST_DATA_DIR}/{DIR_NAME}/{ii}'
    return result
