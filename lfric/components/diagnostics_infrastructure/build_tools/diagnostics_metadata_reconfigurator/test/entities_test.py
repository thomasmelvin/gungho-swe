##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Unit tests for entities module"""
from entities import NonSpatialDimension


def test_non_spatial_dim_add_unit(caplog):
    """Testing that add_unit() works correctly"""
    nsd_1 = NonSpatialDimension(name='test_dim',
                                definition=['1', '2', '3'],
                                unit='test_unit')

    nsd_1.add_unit('test_unit')

    nsd_2 = NonSpatialDimension(name='test_dim',
                                definition=['1', '2', '3'],
                                unit=None)

    nsd_2.add_unit('test_unit')

    assert nsd_1.unit == 'test_unit'
    assert nsd_2.unit == 'test_unit'
    assert caplog.text == ''


def test_non_spatial_dim_add_unit_2(caplog):
    """Testing that add_unit() provides the correct errors"""
    nsd_1 = NonSpatialDimension(name='test_dim',
                                definition=['1', '2', '3'],
                                unit='correct_unit')

    nsd_1.add_unit('not_the_correct_unit')

    assert nsd_1.unit == 'correct_unit'
    assert ("There is a mismatch in non-spatial dimension units of measure."
            in caplog.text)
    assert "test_dim correct_unit not_the_correct_unit" in caplog.text
