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

from mock import Mock, patch
from dao2caom2 import composable, dao_name, COLLECTION

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
        assert (
                test_storage.lineage ==
                f'sky_camera_image/ad:{COLLECTION}/{test_f_name}'
        ), 'wrong lineage'
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
                test_storage.lineage ==
                f'sky_camera_image/ad:{COLLECTION}/{test_f_name}'
        ), 'wrong lineage'
        assert (
                test_storage.source_names ==
                ['vos:goliaths/DAOTest/sky_cam_image.fits.gz']
        ), 'wrong source names'
        assert (
                test_storage.destination_uris ==
                ['ad:DAO/sky_cam_image.fits.gz']
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
