# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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
import re

from os.path import basename
from urllib.parse import urlparse

from caom2pipe import manage_composable as mc

__all__ = ['COLLECTION', 'DAOName', 'PRODUCT_COLLECTION']


COLLECTION = 'DAO'
PRODUCT_COLLECTION = 'DAOCADC'


class DAOName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support compressed and uncompressed product files
    - uncompressed product files have an extension added to the input
      names - lower case if there's a single input, upper case if there's
      multiple inputs, as far as I can tell right now.

    DB 2-11-21
    Processed observation IDs are the same raw observation IDs, even though
    the collections are different. That way, if someone just searches for a
    specific observation ID but doesn’t specify the collection then they
    will see both the unprocessed and processed dataset (if it exists) listed
    in the results page.
    """

    DAO_NAME_PATTERN = '*'

    def __init__(
        self,
        entry=None,
    ):
        self._file_name = mc.CaomName.extract_file_name(entry).replace(
            '.header', ''
        )
        obs_id = DAOName.get_obs_id(self._file_name)
        file_id = mc.StorageName.remove_extensions(self._file_name)
        self._collection = (
            PRODUCT_COLLECTION if DAOName.is_processed(entry) else COLLECTION
        )
        super(DAOName, self).__init__(
            obs_id,
            self._collection,
            DAOName.DAO_NAME_PATTERN,
            self._file_name,
            entry=entry,
            compression='',
        )
        self._file_id = file_id
        self._source_names = [entry]
        self._destination_uris = [self.file_uri]
        self._logger = logging.getLogger(__name__)
        self._logger.debug(self)

    @property
    def file_id(self):
        return self._file_id

    @property
    def file_name(self):
        return self._file_name

    @property
    def file_uri(self):
        return f'{self.scheme}:{self.collection}/{self.file_name}'

    def is_valid(self):
        return True

    @property
    def prev(self):
        """The preview file name for the file."""
        return '{}_1024.png'.format(self.file_id)

    @property
    def product_id(self):
        if self._file_id.startswith('d'):
            return self._file_id
        else:
            return 'sky_camera_image'

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return '{}_256.png'.format(self.file_id)

    @staticmethod
    def get_obs_id(file_name):
        # observation ID differs from the file ID for processed data, except
        # for composite processed observations (master biases and flats)
        file_id = mc.StorageName.remove_extensions(file_name)
        obs_id = file_id
        if re.match('dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aev]', file_id):
            obs_id = file_id[0:-2]
        return obs_id

    @staticmethod
    def is_derived(entry):
        # entry is a uri
        result = False
        if re.match('ad:DAOCADC/dao_[c]\\d{3}_\\d{4}_\\d{6}_[BF].\\w', entry):
            result = True
        return result

    @staticmethod
    def is_master_bias(entry):
        return '_B' in entry

    @staticmethod
    def is_master_flat(entry):
        return '_F' in entry

    @staticmethod
    def is_processed(entry):
        file_name = mc.CaomName.extract_file_name(entry)
        file_id = DAOName.remove_extensions(file_name)
        result = False
        if re.match(
            'dao_[cr]\\d{3}_\\d{4}_\\d{6}_[aevBF]', file_id
        ) or re.match('dao_[p]\\d{3}_\\d{6}(u|v|y|r|i|)', file_id):
            result = True
        return result

    @staticmethod
    def is_unprocessed_reticon(entry):
        file_name = mc.CaomName.extract_file_name(entry)
        file_id = DAOName.remove_extensions(file_name)
        result = False
        if re.match('dao_[r]\\d{3}_\\d{4}_\\d{6}', file_id):
            result = True
        return result

    @staticmethod
    def override_provenance(entry):
        # entry is a uri
        result = False
        if re.match('ad:DAOCADC/dao_[c]\\d{3}_\\d{4}_\\d{6}_[aBF].\\w', entry):
            result = True
        return result
