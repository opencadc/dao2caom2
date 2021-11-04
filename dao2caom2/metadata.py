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

from dataclasses import dataclass
from os.path import exists, join

from caom2utils import data_util
from caom2 import DataProductType
from caom2pipe import client_composable as clc
from caom2pipe import manage_composable as mc

__all__ = ['DefiningMetadata', 'DefiningMetadataFinder']


@dataclass
class DefiningMetadata:
    """The metadata that is required to know 'what to do next' for
    ingesting a file.

    The data_product_type changes the rules for most metadata calculation.
    """
    data_product_type: DataProductType
    uri: str


class DefiningMetadataFinder:
    def __init__(self, clients, config):
        self._clients = clients
        self._data_sources = config.data_sources
        self._use_local_files = config.use_local_files
        self._connected = mc.TaskType.SCRAPE not in config.task_types
        self._logger = logging.getLogger(self.__class__.__name__)

    def _get_data_product_type(self, headers):
        obs_mode = headers[0].get('OBSMODE')
        # DB 30-04-20
        # I added an obsmodes hash to allow it to be imaging or Imaging but
        # all should be 'Imaging'.
        data_product_type = DataProductType.IMAGE
        if '-slit' in obs_mode:
            data_product_type = DataProductType.SPECTRUM
        self._logger.debug('End _get_data_product_type')
        return data_product_type

    def check_caom2(self, uri):
        self._logger.debug(f'Begin check_caom2 for {uri}')
        query_string = f"""
        SELECT P.dataProductType
        FROM caom2.Observation AS O
        JOIN caom2.Plane AS P on P.obsID = O.obsID
        JOIN caom2.Artifact AS A on A.planeID = P.planeID
        WHERE A.uri = '{uri}' 
        """
        table = clc.query_tap_client(query_string, self._clients.query_client)
        result = None
        if len(table) == 1:
            result = DefiningMetadata(table[0]['dataProductType'], uri)
        self._logger.debug('End check_caom2')
        return result

    def check_local(self, uri):
        self._logger.debug(f'Begin check_local for {uri}')
        ignore_scheme, ignore_path, f_name = mc.decompose_uri(uri)
        result = None
        for entry in self._data_sources:
            fqn = join(entry, f_name)
            if exists(fqn):
                self._logger.debug(f'Looking in {fqn} for headers.')
                headers = data_util.get_local_headers_from_fits(fqn)
                result = DefiningMetadata(
                    self._get_data_product_type(headers), uri
                )
                break
        self._logger.debug('End check_local')
        return result

    def check_remote(self, uri):
        self._logger.debug(f'Begin check_remote for {uri}')
        headers = self._clients.data_client.get_head(uri)
        data_product_type = self._get_data_product_type(headers)
        self._logger.debug('End check_remote')
        return DefiningMetadata(data_product_type, uri)

    def get(self, uri):
        """
        :param uri: str CADC storage reference
        :return: DefiningMetadata, if it can be found by a service, for the
           uri.
        """
        if self._connected:
            result = None
            if self._use_local_files:
                result = self.check_local(uri)
            if result is None:
                result = self.check_caom2(uri)
            if result is None:
                result = self.check_remote(uri)
        else:
            result = self.check_local(uri)
        return result
