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

from mock import Mock, patch

from cadcdata import FileInfo
from caom2pipe import manage_composable as mc
from dao2caom2 import data_source


@patch('caom2pipe.astro_composable.check_fits', autospec=True)
@patch('caom2pipe.client_composable.vault_info', autospec=True)
def test_dao_transfer_check_fits_verify(vault_info_mock, check_fits_mock):
    test_match_file_info = FileInfo(
        id='vos:abc/def.fits',
        md5sum='ghi',
    )
    test_different_file_info = FileInfo(
        id='vos:abc/def.fits',
        md5sum='ghi',
    )
    test_file_info = [test_match_file_info, test_different_file_info]

    test_data_client = Mock(autospec=True)
    test_vos_client = Mock(autospec=True)

    test_config = mc.Config()
    test_config.data_source_extensions = ['.fits.gz']
    test_config.data_sources = ['vos:DAO/Archive/Incoming']
    test_config.cleanup_failure_destination = 'vos:DAO/failure'
    test_config.cleanup_success_destination = 'vos:DAO/success'

    def _mock_listdir(entry):
        if entry.endswith('Incoming'):
            return ['dao123.fits.gz', 'dao456.fits', 'Yesterday']
        else:
            return []

    test_vos_client.listdir.side_effect = _mock_listdir
    test_vos_client.isdir.side_effect = [
        False, False, True, False, False, True
    ]

    test_result_1 = False
    test_result_2 = False

    vault_info_mock.return_value = test_file_info
    test_data_client.info.return_value = test_file_info

    for case in [True, False]:
        check_fits_mock.return_value = not case
        test_config.cleanup_files_when_storing = case

        test_subject = data_source.DAOVaultDataSource(
            test_config, test_vos_client, test_data_client
        )
        assert test_subject is not None, 'expect ctor to work'
        test_result = test_subject.get_work()

        assert test_result is not None, 'expect a work list'
        if not case:
            assert len(test_result) == 1, 'wrong work list entries'
            assert (
                test_result[0] == 'vos:DAO/Archive/Incoming/dao123.fits.gz'
            ), 'wrong work entry'
            test_result_1 = True
        else:
            assert len(test_result) == 0, 'check fits should fail'
            assert test_vos_client.move.called, 'move to failure'
            test_vos_client.move.assert_called_with(
                'vos:DAO/Archive/Incoming/dao123.fits.gz',
                'vos:DAO/failure/dao123.fits.gz',
            ), 'wrong move args'
            test_result_2 = True

        assert test_vos_client.isdir.call_count == 3, 'wrong is_dir count'
        test_vos_client.isdir.reset_mock()

    assert test_result_1, 'followed one result path at least once'
    assert test_result_2, 'followed no result path at least once'

    # test the case when the md5sums are the same, so the transfer does
    # not occur, but the file ends up in the success location
    test_vos_client.isdir.side_effect = [False, False, True]
    test_config.cleanup_files_when_storing = True
    test_config.store_modified_files_only = True
    check_fits_mock.return_value = True
    vault_info_mock.return_value = test_match_file_info
    test_data_client.info.return_value = test_different_file_info

    second_test_subject = data_source.DAOVaultDataSource(
        test_config, test_vos_client, test_data_client
    )
    assert second_test_subject is not None, 'second ctor fails'
    second_test_result = second_test_subject.get_work()
    assert second_test_result is not None, 'expect a second result'
    assert len(second_test_result) == 0, 'should be no successes'
    assert test_vos_client.move.called, 'expect a success move call'
    test_vos_client.move.assert_called_with(
        'vos:DAO/Archive/Incoming/dao123.fits.gz',
        'vos:DAO/success/dao123.fits.gz',
    ), 'wrong success move args'
