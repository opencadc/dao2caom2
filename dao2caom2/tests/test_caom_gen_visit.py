# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

from os import listdir, unlink
from os.path import basename, dirname, exists, join, realpath

from caom2.diff import get_differences
from cadcdata import FileInfo
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from dao2caom2 import DAOName
from dao2caom2 import fits2caom2_augmentation

THIS_DIR = dirname(realpath(__file__))
TEST_DATA_DIR = join(THIS_DIR, 'data')


def pytest_generate_tests(metafunc):
    files = []
    if exists(TEST_DATA_DIR):
        files = [join(TEST_DATA_DIR, name) for name in listdir(TEST_DATA_DIR) if name.endswith('header')]
    metafunc.parametrize('test_name', files)


def test_visitor(test_name, test_config, tmp_path):
    test_config.change_working_directory(tmp_path.as_posix())
    dao_name = DAOName([test_name.replace('.header', '.gz')])
    file_info = FileInfo(id=dao_name.file_uri, file_type='application/fits')
    headers = ac.make_headers_from_file(test_name)
    dao_name.metadata = {dao_name.file_uri: headers}
    dao_name.file_info = {dao_name.file_uri: file_info}
    test_reporter = mc.ExecutionReporter2(test_config)
    kwargs = {
        'storage_name': dao_name,
        'config': test_config,
        'reporter': test_reporter,
    }
    expected_fqn = f'{TEST_DATA_DIR}/{dao_name.file_id}.expected.xml'
    actual_fqn = expected_fqn.replace('expected', 'actual')
    if exists(actual_fqn):
        unlink(actual_fqn)
    observation = None
    observation = fits2caom2_augmentation.visit(observation, **kwargs)

    if observation is None:
        assert False, f'observation construction failed for {test_name}'
    else:
        if exists(expected_fqn):
            expected = mc.read_obs_from_file(expected_fqn)
            compare_result = get_differences(expected, observation)
            if compare_result is not None:
                mc.write_obs_to_file(observation, actual_fqn)
                compare_text = '\n'.join([r for r in compare_result])
                msg = f'Differences found in observation {expected.observation_id}\n{compare_text}'
                raise AssertionError(msg)
        else:
            mc.write_obs_to_file(observation, actual_fqn)
            raise AssertionError(f'{expected_fqn} does not exist.')
