##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Tests for the diagnostics_metadata_collection module"""
import pytest
from diagnostics_metadata_collection import Metadata


def test_get_non_spatial_dimension_1():
    """Tests get_non_spatial_dimension functions when used correctly"""

    dummy_metadata_object = Metadata()

    dummy_metadata_object._non_spatial_dimensions = {
        "test_1": "It's over",
        "test_2": "9000!"
    }

    result_1 = dummy_metadata_object.get_non_spatial_dimension("test_1")
    result_2 = dummy_metadata_object.get_non_spatial_dimension("test_2")

    assert result_1 == "It's over"
    assert result_2 == "9000!"


def test_get_non_spatial_dimension_2(caplog):
    """Tests that get_non_spatial_dimension returns correct error code when
     trying to a non_spatial_dimension that does not exist"""

    dummy_metadata_object = Metadata()

    dummy_metadata_object._non_spatial_dimensions = {
        "test_1": "It's over",
        "test_2": "9000!"
    }
    with pytest.raises(KeyError):
        dummy_metadata_object.get_non_spatial_dimension(
            "What's the scouter say about his power level?")

    assert ("Metadata mismatch: What's the scouter say about his power level?"
            in caplog.text)
