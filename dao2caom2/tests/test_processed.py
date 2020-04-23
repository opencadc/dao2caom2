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
import pytest

from dao2caom2 import main_app, APPLICATION, COLLECTION, DAOName
from caom2.diff import get_differences
from caom2pipe import manage_composable as mc

import os
import sys

from mock import patch

pytest.main(args=['-s', os.path.abspath(__file__)])
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'processed'
LOOKUP = {'dao_c122_2007_000882': ['dao_c122_2007_000882_v']
          }


def pytest_generate_tests(metafunc):
    obs_id_list = LOOKUP.keys()
    metafunc.parametrize('test_name', obs_id_list)


def test_multi_plane(test_name):
    dao_name = DAOName(file_name='{}.fits.header'.format(LOOKUP[test_name][0]))
    lineage = _get_lineage(dao_name.obs_id)
    actual_fqn = '{}/{}/{}.actual.xml'.format(
        TEST_DATA_DIR, DIR_NAME, dao_name.obs_id)

    local = _get_local(test_name)
    plugin = PLUGIN

    expected_fqn = '{}/{}/{}.expected.xml'.format(
        TEST_DATA_DIR, DIR_NAME, dao_name.obs_id)
    expected = mc.read_obs_from_file(expected_fqn)

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock:
        def gfi(archive, file_id):
            fid = mc.StorageName.remove_extensions(file_id)
            temp = expected.planes[fid].artifacts[
                'ad:{}/{}'.format(archive, file_id)]
            return {'size': temp.content_length,
                    'md5sum': temp.content_checksum.uri,
                    'type': temp.content_type}

        data_client_mock.return_value.get_file_info.side_effect = gfi

        if os.path.exists(actual_fqn):
            os.remove(actual_fqn)

        sys.argv = \
            ('{} --quiet --no_validate --local {} --observation {} {} '
             '--plugin {} --module {} --out {} --lineage {}'.
             format(APPLICATION, local, COLLECTION, dao_name.obs_id, plugin,
                    plugin, actual_fqn, lineage)).split()
        print(sys.argv)
        main_app()
        actual = mc.read_obs_from_file(actual_fqn)
        result = get_differences(expected, actual, 'Observation')
        if result:
            msg = 'Differences found obs id {} \n{}'.format(
                expected.observation_id, '\n'.join([r for r in result]))
            raise AssertionError(msg)
        # assert False  # cause I want to see logging messages


def _get_lineage(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        fits = mc.get_lineage(COLLECTION, ii, '{}.fits.gz'.format(ii))
        result = '{} {}'.format(result, fits)
    return result


def _get_local(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        result = '{} {}/{}/{}.fits.header'.format(result, TEST_DATA_DIR,
                                                  DIR_NAME, ii)
    return result
