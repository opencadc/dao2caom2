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

import os
import test_fits2caom2_augmentation

from collections import deque
from mock import Mock, patch
from tempfile import TemporaryDirectory

from cadcdata import FileInfo
from caom2 import SimpleObservation
from caom2pipe import manage_composable as mc
from dao2caom2 import composable, dao_name

F_NAME_LIST = [
    'data_report.txt',
    'failure_log.txt',
    'rejected.yml',
    'retries.txt',
    'success_log.txt',
]


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run(run_mock, access_mock):
    access_mock.return_value = 'https://localhost'
    test_f_id = 'test_file_id'
    test_obs_id = test_f_id
    test_f_name = f'{test_f_id}.fits'
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_fits2caom2_augmentation.TEST_DATA_DIR)
    try:
        # execution
        composable._run()
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(
            test_storage, dao_name.DAOName
        ), type(test_storage)
        assert test_storage.obs_id == test_obs_id, 'wrong obs id'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert (
                test_storage.fname_on_disk == test_f_name
        ), 'wrong fname on disk'
        assert test_storage.url is None, 'wrong url'
    finally:
        os.getcwd = getcwd_orig
        # clean up the summary report text file
        # clean up the files created as a by-product of a run
        for f_name in F_NAME_LIST:
            fqn = os.path.join(
                test_fits2caom2_augmentation.TEST_DATA_DIR, f_name
            )
            if os.path.exists(fqn):
                os.unlink(fqn)


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('dao2caom2.composable.Client', autospec=True)
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_vo(run_mock, vo_client_mock, access_mock):
    access_mock.return_value = 'https://localhost'
    test_obs_id = 'sky_cam_image'
    test_f_name = f'{test_obs_id}.fits.gz'
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_fits2caom2_augmentation.TEST_DATA_DIR)
    vo_client_mock.return_value.listdir.return_value = [
        'sky_cam_image.fits.gz',
    ]
    vo_client_mock.return_value.isdir.return_value = False

    try:
        # execution
        composable._run_vo()
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(
            test_storage, dao_name.DAOName
        ), type(test_storage)
        assert test_storage.obs_id == test_obs_id, 'wrong obs id'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert (
                test_storage.fname_on_disk == test_f_name
        ), 'wrong fname on disk'
        assert test_storage.url is None, 'wrong url'
        assert (
                test_storage.source_names ==
                ['vos:goliaths/DAOTest/sky_cam_image.fits.gz']
        ), 'wrong source names'
        assert (
                test_storage.destination_uris ==
                ['cadc:DAO/sky_cam_image.fits.gz']
        ), 'wrong destination uris'
    finally:
        os.getcwd = getcwd_orig
        # clean up the summary report text file
        # clean up the files created as a by-product of a run
        for f_name in F_NAME_LIST:
            fqn = os.path.join(
                test_fits2caom2_augmentation.TEST_DATA_DIR, f_name
            )
            if os.path.exists(fqn):
                os.unlink(fqn)


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_headers')
@patch('caom2pipe.reader_composable.FileMetadataReader._retrieve_file_info')
@patch('caom2pipe.data_source_composable.LocalFilesDataSource.clean_up')
@patch('caom2pipe.data_source_composable.LocalFilesDataSource.get_work')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_store_ingest(
    data_client_mock,
    repo_client_mock,
    get_work_mock,
    cleanup_mock,
    reader_file_info_mock,
    reader_headers_mock,
    access_mock,
):
    access_mock.return_value = 'https://localhost'
    temp_deque = deque()
    temp_deque.append('/data/abc.fits')
    get_work_mock.return_value = temp_deque
    repo_client_mock.return_value.read.return_value = None
    reader_headers_mock.return_value = []
    reader_file_info_mock.return_value = FileInfo(
        id='cadc:DAOCADC/abc.fits',
        file_type='application/fits',
        md5sum='md5:def',
    )

    cwd = os.getcwd()
    with TemporaryDirectory() as tmp_dir_name:
        os.chdir(tmp_dir_name)
        test_config = mc.Config()
        test_config.working_directory = tmp_dir_name
        test_config.task_types = [mc.TaskType.STORE, mc.TaskType.INGEST]
        test_config.use_local_files = True
        test_config.cleanup_files_when_storing = True
        test_config.cleanup_failure_destination = '/data/failure'
        test_config.cleanup_success_destination = '/data/success'
        test_config.data_sources = ['/data']
        test_config.data_source_extensions = ['.fits']
        test_config.logging_level = 'INFO'
        test_config.proxy_file_name = 'cadcproxy.pem'
        test_config.proxy_fqn = f'{tmp_dir_name}/cadcproxy.pem'
        test_config.features.supports_latest_client = True
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=tmp_dir_name)
        try:
            test_result = composable._run()
            assert test_result is not None, 'expect result'
            assert test_result == 0, 'expect success'
            assert repo_client_mock.return_value.read.called, 'read called'
            # make sure data is not really being written to CADC storage :)
            assert (
                data_client_mock.return_value.put.called
            ), 'put should be called'
            assert (
                    data_client_mock.return_value.put.call_count == 1
            ), 'wrong number of puts'
            data_client_mock.return_value.put.assert_called_with(
                '/data', 'cadc:DAO/abc.fits', None
            )
            assert cleanup_mock.called, 'cleanup'
            cleanup_mock.assert_called_with(
                '/data/abc.fits', 0
            ), 'wrong cleanup args'
            assert repo_client_mock.return_value.create.called, 'create call'
            args, kwargs = repo_client_mock.return_value.create.call_args
            assert isinstance(args[0], SimpleObservation), 'wrong arg type'
        finally:
            os.getcwd = getcwd_orig
            os.chdir(cwd)
