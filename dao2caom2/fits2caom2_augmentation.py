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

from caom2pipe import caom_composable as cc
from caom2 import DataProductType
from dao2caom2 import telescopes, dao_name


class DAOFits2caom2Visitor(cc.Fits2caom2VisitorRunnerMeta):

    @staticmethod
    def get_data_product_type(headers):
        obs_mode = headers[0].get('OBSMODE')
        # DB 30-04-20
        # I added an obsmodes hash to allow it to be imaging or Imaging but
        # all should be 'Imaging'.
        data_product_type = DataProductType.IMAGE
        if '-slit' in obs_mode:
            data_product_type = DataProductType.SPECTRUM
        return data_product_type

    def _get_mappings(self, dest_uri):
        if self._storage_name.file_name.startswith('a'):
            result = telescopes.SkyCam(
                self._storage_name, self._clients, self._reporter, self._observation, self._config
            )
        else:
            data_product_type = DAOFits2caom2Visitor.get_data_product_type(self._storage_name.metadata.get(dest_uri))
            if data_product_type == DataProductType.IMAGE:
                if dao_name.DAOName.is_processed(self._storage_name.file_id):
                    if self._storage_name.is_12_metre:
                        result = telescopes.Dao12MetreProcessedImage(
                            self._storage_name,
                            self._clients,
                            self._reporter,
                            self._observation,
                            self._config,
                        )
                    else:
                        result = telescopes.Dao18MetreProcessedImage(
                            self._storage_name,
                            self._clients,
                            self._reporter,
                            self._observation,
                            self._config,
                        )
                elif self._storage_name.is_12_metre:
                    result = telescopes.Dao12MetreImage(
                        self._storage_name, self._clients, self._reporter, self._observation, self._config
                    )
                else:
                    result = telescopes.Dao18MetreImage(
                        self._storage_name, self._clients, self._reporter, self._observation, self._config
                    )
            else:
                if dao_name.DAOName.is_processed(self._storage_name.file_id):
                    if self._storage_name.is_12_metre:
                        result = telescopes.Dao12MetreProcessedSpectrum(
                            self._storage_name,
                            self._clients,
                            self._reporter,
                            self._observation,
                            self._config,
                        )
                    else:
                        result = telescopes.Dao18MetreProcessedSpectrum(
                            self._storage_name,
                            self._clients,
                            self._reporter,
                            self._observation,
                            self._config,
                        )
                elif self._storage_name.is_12_metre:
                    result = telescopes.Dao12MetreSpectrum(
                        self._storage_name, self._clients, self._reporter, self._observation, self._config
                    )
                else:
                    result = telescopes.Dao18MetreSpectrum(
                        self._storage_name, self._clients, self._reporter, self._observation, self._config
                    )

        self._logger.debug(f'Created {result.__class__.__name__} instance.')
        return [result]


def visit(observation, **kwargs):
    return DAOFits2caom2Visitor(observation, **kwargs).visit()
