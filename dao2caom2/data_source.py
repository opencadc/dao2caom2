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

import logging
import traceback

from os.path import basename, join
from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe import manage_composable as mc


class DAOVaultDataSource(dsc.VaultListDirDataSource):
    def __init__(self, config, vo_client, cadc_client, recursive=True):
        super(DAOVaultDataSource, self).__init__(vo_client, config)
        self._cleanup_when_storing = config.cleanup_files_when_storing
        self._cleanup_failure_directory = config.cleanup_failure_destination
        self._cleanup_success_directory = config.cleanup_success_destination
        self._store_modified_files_only = config.store_modified_files_only
        self._supports_latest_client = config.features.supports_latest_client
        self._archive = config.archive
        self._collection = config.collection
        self._recursive = recursive
        self._cadc_client = cadc_client
        self._work = []
        self._logger = logging.getLogger(self.__class__.__name__)

    def clean_up(self):
        if self._cleanup_when_storing:
            for fqn in self._work:
                check_result, vos_meta = self._check_md5sum(fqn)
                if check_result:
                    # if vos_meta is None, it's already been cleaned up,
                    # due to astropy fits verify failure cleanup
                    if vos_meta is not None:
                        # the transfer itself failed, so track as a failure
                        self._move_action(
                            fqn, self._cleanup_failure_directory
                        )
                else:
                    self._move_action(fqn, self._cleanup_success_directory)

    def default_filter(self, entry, entry_fqn):
        copy_file = False
        for extension in self._data_source_extensions:
            if entry.endswith(extension):
                if entry.startswith('.'):
                    # skip dot files
                    copy_file = False
                elif self._store_modified_files_only:
                    # only transfer files with a different MD5 checksum
                    copy_file, ignore_meta = self._check_md5sum(entry_fqn)
                    if not copy_file and self._cleanup_when_storing:
                        self._move_action(
                            entry_fqn, self._cleanup_success_directory
                        )
                else:
                    copy_file = True
                break
        return copy_file

    def get_work(self):
        self._logger.debug(f'Begin get_work.')
        for source in self._source_directories:
            self._logger.info(f'Look in {source} for work.')
            self._find_work(source)
        self._logger.debug('End get_work')
        return self._work

    def _check_md5sum(self, entry_fqn):
        # get the metadata from VOS
        result = True
        vos_meta = clc.vault_info(self._vo_client, entry_fqn)
        # get the metadata at CADC
        f_name = basename(entry_fqn)
        scheme = 'cadc' if self._supports_latest_client else 'ad'
        cadc_name = mc.build_uri(self._collection, f_name, scheme)
        cadc_meta = self._cadc_client.info(cadc_name)
        if cadc_meta is not None and vos_meta.md5sum == cadc_meta.md5sum:
            self._logger.warning(
                f'{entry_fqn} has the same md5sum at CADC. Not transferring.'
            )
            result = False
        return result, vos_meta

    def _find_work(self, entry):
        dir_listing = self._vo_client.listdir(entry)
        for dir_entry in dir_listing:
            dir_entry_fqn = f'{entry}/{dir_entry}'
            if self._vo_client.isdir(dir_entry_fqn) and self._recursive:
                self._find_work(dir_entry_fqn)
            else:
                if self.default_filter(dir_entry, dir_entry_fqn):
                    self._logger.info(f'Adding {dir_entry_fqn} to work list.')
                    self._work.append(dir_entry_fqn)

    def _move_action(self, fqn, destination):
        """
        :param fqn: VOS URI, includes a file name
        :param destination: VOS URI, points to a directory
        :return:
        """
        # if move when storing is enabled, move to an after-action location
        if self._cleanup_when_storing:
            try:
                # if the destination is a fully-qualified name, an
                # over-write will succeed
                f_name = basename(fqn)
                dest_fqn = join(destination, f_name)
                if self._vo_client.status(dest_fqn):
                    self._logger.warning(
                        f'Removing {dest_fqn} prior to over-write.'
                    )
                    self._vo_client.delete(dest_fqn)
                self._logger.warning(f'Moving {fqn} to {dest_fqn}')
                self._vo_client.move(fqn, dest_fqn)
            except Exception as e:
                self._logger.debug(traceback.format_exc())
                self._logger.error(
                    f'Failed to move {fqn} to {destination}'
                )
                raise mc.CadcException(e)
