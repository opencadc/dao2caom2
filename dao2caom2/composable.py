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

"""
Application to create CAOM2 observations from DAO FITS files. Based on
code and configuration from wcaom2archive/dao2caom2.
"""

import logging
import sys
import traceback

from vos import Client

from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe import name_builder_composable as nbc
from caom2pipe import manage_composable as mc
from caom2pipe import reader_composable as rdc
from caom2pipe import run_composable as rc
from dao2caom2 import dao_name, preview_augmentation
from dao2caom2 import cleanup_augmentation, data_source, transfer
from dao2caom2 import fits2caom2_augmentation

DAO_BOOKMARK = 'dao_timestamp'


META_VISITORS = [fits2caom2_augmentation, cleanup_augmentation]
DATA_VISITORS = [preview_augmentation]


def _common():
    config = mc.Config()
    config.get_executors()
    clients = clc.ClientCollection(config)
    metadata_reader = None
    if config.use_local_files:
        metadata_reader = rdc.FileMetadataReader()
    name_builder = nbc.EntryBuilder(dao_name.DAOName)
    return config, clients, name_builder, metadata_reader


def _run():
    """
    Uses a todo file to identify the work to be done.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config, clients, name_builder, metadata_reader = _common()
    files_source = None
    if config.use_local_files:
        if config.cleanup_files_when_storing:
            files_source = data_source.DAOLocalFilesDataSource(
                config, clients.data_client, metadata_reader
            )
    else:
        files_source = dsc.TodoFileDataSource(config)
    return rc.run_by_todo(
        name_builder=name_builder,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        clients=clients,
        config=config,
        sources=[files_source],
        metadata_reader=metadata_reader,
    )


def run():
    """Wraps _run in exception handling, with sys.exit calls."""
    try:
        result = _run()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_vo():
    """
    Uses a VOS listdir to identify the work to be done.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config, clients, name_builder, metadata_reader = _common()
    vos_client = Client(vospace_certfile=config.proxy_file_name)
    if metadata_reader is None:
        metadata_reader = rdc.VaultReader(vos_client)
    clients.vo_client = vos_client
    source = data_source.DAOVaultDataSource(
        config, clients.vo_client, clients.data_client, metadata_reader
    )
    store_transferrer = transfer.VoFitsCleanupTransfer(vos_client, config)
    return rc.run_by_todo(
        name_builder=name_builder,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        sources=[source],
        clients=clients,
        store_transfer=store_transferrer,
        metadata_reader=metadata_reader,
    )


def run_vo():
    """Wraps _run_vo in exception handling, with sys.exit calls."""
    try:
        result = _run_vo()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_state():
    """Uses a state file with a timestamp to control which entries will be
    processed.
    """
    config, clients, name_builder, metadata_reader = _common()
    files_source = None
    if config.use_local_files:
        if config.cleanup_files_when_storing:
            files_source = data_source.DAOLocalFilesDataSource(
                config, clients.data_client, metadata_reader
            )
    else:
        files_source = dsc.ListDirTimeBoxDataSource(config)
    return rc.run_by_state(
        name_builder=name_builder,
        bookmark_name=DAO_BOOKMARK,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        sources=[files_source],
        clients=clients,
        metadata_reader=metadata_reader,
    )


def run_state():
    """Wraps _run_state in exception handling."""
    try:
        _run_state()
        sys.exit(0)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
