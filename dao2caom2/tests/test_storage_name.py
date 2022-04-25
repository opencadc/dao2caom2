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

from caom2pipe import name_builder_composable as nbc
from dao2caom2 import DAOName, COLLECTION, PRODUCT_COLLECTION


def test_is_valid():
    assert DAOName('anything').is_valid()


def test_processed():
    test_subject = DAOName('dao_c122_2020_004100_v.fits.gz')
    assert test_subject is not None, 'expect a value'
    assert test_subject.obs_id == 'dao_c122_2020_004100', 'wrong obs id'
    assert (
        test_subject.file_name == 'dao_c122_2020_004100_v.fits.gz'
    ), 'wrong file name'
    assert (
        test_subject.product_id == 'dao_c122_2020_004100_v'
    ), 'wrong product id'
    assert test_subject.file_id == 'dao_c122_2020_004100_v', 'wrong file id'
    assert test_subject.source_names == [
        'dao_c122_2020_004100_v.fits.gz'
    ], 'wrong source names'
    assert test_subject.destination_uris == [
        f'cadc:{PRODUCT_COLLECTION}/dao_c122_2020_004100_v.fits'
    ], 'wrong destination uris'
    assert (
        test_subject.file_uri
        == f'cadc:{PRODUCT_COLLECTION}/dao_c122_2020_004100_v.fits'
    ), 'wrong file uri'
    assert test_subject.collection == PRODUCT_COLLECTION, 'wrong collection'


def test_raw():
    test_result = DAOName('cadc:DAO/dao_c182_2018_015013.fits')
    assert test_result is not None, 'expect a result'
    assert test_result.obs_id == 'dao_c182_2018_015013'
    assert test_result.file_name == 'dao_c182_2018_015013.fits'
    assert test_result.file_id == 'dao_c182_2018_015013'
    assert test_result.source_names == ['cadc:DAO/dao_c182_2018_015013.fits']
    assert test_result.destination_uris == [
        'cadc:DAO/dao_c182_2018_015013.fits'
    ], 'wrong destination uris'

    test_result_2 = DAOName(
        '/usr/src/app/dao_c182_2018_015013/dao_c182_2018_015013.fits.gz'
    )
    assert test_result_2.source_names == [
        '/usr/src/app/dao_c182_2018_015013/dao_c182_2018_015013.fits.gz'
    ]
    assert test_result_2.destination_uris == [
        f'cadc:{COLLECTION}/dao_c182_2018_015013.fits'
    ], 'wrong destination uris'
