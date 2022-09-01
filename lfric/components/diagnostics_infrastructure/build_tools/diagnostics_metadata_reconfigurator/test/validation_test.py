##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Test metadata_validator"""

import logging
import pytest

from entities import Field, FieldGroup, VerticalDimension
from diagnostics_metadata_collection import Metadata
from metadata_validator import validate_metadata, InvalidMetadataError


def test_validate_metadata_valid(caplog):
    """
    Check that function validate_metadata does not return an error when
    passed valid diagnostic field metadata
    """

    metadata = Metadata()

    field_group_name = "TestFieldGroup"
    metadata.add_field_group(FieldGroup(name=field_group_name))

    # Fields with valid naming metadata
    valid_field1 = Field(unique_id="valid_field1",
                         field_group_id=field_group_name,
                         long_name='valid_long_name1')
    valid_field2 = Field(unique_id="valid_field2",
                         field_group_id=field_group_name,
                         standard_name="valid_standard_name2")
    valid_field3 = Field(unique_id="valid_field3",
                         field_group_id=field_group_name,
                         long_name='valid_long_name3',
                         standard_name="valid_standard_name")

    metadata.add_field(valid_field1)
    metadata.add_field(valid_field2)
    metadata.add_field(valid_field3)

    validate_metadata(metadata)
    assert "'valid_field1" not in caplog.text
    assert "'valid_field2" not in caplog.text
    assert "'valid_field3" not in caplog.text


def test_validate_metadata_invalid_names(caplog):
    """
    Check that metadata is invalid when long name and standard name are
    not specified, or set to either None or an empty string
    """

    metadata = Metadata()

    field_group_name = "TestFieldGroup"
    metadata.add_field_group(FieldGroup(name=field_group_name))

    invalid_field1 = Field(unique_id="invalid_field1",
                           field_group_id=field_group_name,)
    invalid_field2 = Field(unique_id="invalid_field2",
                           field_group_id=field_group_name,
                           long_name="",
                           standard_name="")
    invalid_field3 = Field(unique_id="invalid_field3",
                           field_group_id=field_group_name,
                           long_name=None,
                           standard_name=None)

    valid_field1 = Field(unique_id="valid_field",
                         field_group_id=field_group_name,
                         long_name='valid_long_name')
    valid_field2 = Field(unique_id="valid_field2",
                         field_group_id=field_group_name,
                         standard_name="valid_standard_name")

    metadata.add_field(invalid_field1)
    metadata.add_field(invalid_field2)
    metadata.add_field(invalid_field3)
    metadata.add_field(valid_field1)
    metadata.add_field(valid_field2)

    with pytest.raises(InvalidMetadataError):
        validate_metadata(metadata)

    # Invalid fields should be logged as errors, valid fields should not
    assert "'valid_field1" not in caplog.text
    assert "'valid_field2" not in caplog.text

    assert "'invalid_field1" in caplog.text
    assert "'invalid_field2" in caplog.text
    assert "'invalid_field3" in caplog.text


def test_validate_metadata_invalid_names_force(caplog):
    """
    Check that metadata is invalid when long name and standard name are
    not specified or set to either None or an empty string, but that no
    error when is raised when force is set to True
    """

    metadata = Metadata()

    field_group_name = "TestFieldGroup"
    metadata.add_field_group(FieldGroup(name=field_group_name))

    invalid_field1 = Field(unique_id="invalid_field1",
                           field_group_id=field_group_name,)
    invalid_field2 = Field(unique_id="invalid_field2",
                           field_group_id=field_group_name,
                           long_name="",
                           standard_name="")
    invalid_field3 = Field(unique_id="invalid_field3",
                           field_group_id=field_group_name,
                           long_name=None,
                           standard_name=None)

    valid_field1 = Field(unique_id="valid_field",
                         field_group_id=field_group_name,
                         long_name='valid_long_name')
    valid_field2 = Field(unique_id="valid_field2",
                         field_group_id=field_group_name,
                         standard_name="valid_standard_name")

    metadata.add_field(invalid_field1)
    metadata.add_field(invalid_field2)
    metadata.add_field(invalid_field3)
    metadata.add_field(valid_field1)
    metadata.add_field(valid_field2)

    # Suppress logging at warning and below
    with caplog.at_level(logging.ERROR):
        # With force there should should be no exception raised and invalid
        # fields should be warnings and so not appear in the log
        validate_metadata(metadata, True)
        assert "'invalid_field1" not in caplog.text
        assert "'invalid_field2" not in caplog.text
        assert "'invalid_field3" not in caplog.text

        # Without force the exception should be raised and invalid fields
        # logged as errors, so they should appear now
        with pytest.raises(InvalidMetadataError):
            validate_metadata(metadata)

        assert "'valid_field1" not in caplog.text
        assert "'valid_field2" not in caplog.text

        assert "'invalid_field1" in caplog.text
        assert "'invalid_field2" in caplog.text
        assert "'invalid_field3" in caplog.text


def test_validate_metadata_undefined_vertical_dimension(caplog):
    """
    Check that metadata fails validation when vertical dimension for a group
    is not defined in rose-app.conf
    """

    metadata = Metadata()

    field_group_name = "TestFieldGroup"
    metadata.add_field_group(FieldGroup(name=field_group_name))

    vertical_dimension = VerticalDimension(
            unique_id="model_vert_axis_1_half_levels",
            name="model_vert_axis_1_half_levels",
            positive_direction="up",
            primary_axis='false',
            level_definition=[1.5],
            number_of_layers=1,
            units="m")

    invalid_field1 = Field(
            unique_id="invalid_field1",
            field_group_id=field_group_name,
            long_name="invalid_field1",
            vertical_dimension_id="model_vert_axis_2_half_levels")
    valid_field1 = Field(
            unique_id="valid_field1",
            field_group_id=field_group_name,
            long_name="valid_field1",
            vertical_dimension_id="model_vert_axis_1_half_levels")

    metadata.add_vertical_dimension(vertical_dimension)
    metadata.add_field(invalid_field1)
    metadata.add_field(valid_field1)

    with pytest.raises(InvalidMetadataError):
        validate_metadata(metadata)

    # Invalid fields should be logged as errors, valid fields should not
    assert "'valid_field1" not in caplog.text
    assert "'invalid_field1" in caplog.text


def test_validate_metadata_undefined_vertical_dimension_force(caplog):
    """
    Check that the force options prevents an InvalidMetadataError when a
    field's vertical dimension is invalid
    """

    metadata = Metadata()

    field_group_name = "TestFieldGroup"
    metadata.add_field_group(FieldGroup(name=field_group_name))

    vertical_dimension = VerticalDimension(
            unique_id="model_vert_axis_1_half_levels",
            name="model_vert_axis_1_half_levels",
            positive_direction="up",
            primary_axis='false',
            level_definition=[1.5],
            number_of_layers=1,
            units="m")

    invalid_field1 = Field(
            unique_id="invalid_field1",
            field_group_id=field_group_name,
            long_name="invalid_field1",
            vertical_dimension_id="model_vert_axis_2_half_levels")
    valid_field1 = Field(
            unique_id="valid_field1",
            field_group_id=field_group_name,
            long_name="valid_field1",
            vertical_dimension_id="model_vert_axis_1_half_levels")

    metadata.add_vertical_dimension(vertical_dimension)
    metadata.add_field(invalid_field1)
    metadata.add_field(valid_field1)

    # Suppress logging at warning and below
    with caplog.at_level(logging.ERROR):

        # With force there should should be no exception raised and invalid
        # fields should be warnings and so not appear in the log
        validate_metadata(metadata, True)
        assert "'valid_field1" not in caplog.text
        assert "'invalid_field1" not in caplog.text

        # Without force the exception should be raised and invalid fields
        # logged as errors, so they should appear now
        with pytest.raises(InvalidMetadataError):
            validate_metadata(metadata)
        assert "'valid_field1" not in caplog.text
        assert "'invalid_field1" in caplog.text


def test_validate_vert_axes(caplog):
    """Test error raised when there are multiple primary axes defined"""

    axis_1_full = VerticalDimension("model_vert_axis_1_full_levels",
                                    name="model_vert_axis_1_full_levels",
                                    domain_top=42,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_1_half = VerticalDimension("model_vert_axis_1_half_levels",
                                    name="model_vert_axis_1_half_levels",
                                    domain_top=42,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_2_full = VerticalDimension("model_vert_axis_2_full_levels",
                                    name="model_vert_axis_2_full_levels",
                                    domain_top=43,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_2_half = VerticalDimension("model_vert_axis_2_half_levels",
                                    name="model_vert_axis_2_half_levels",
                                    domain_top=43,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')

    metadata = Metadata()
    metadata.add_vertical_dimension(axis_1_full)
    metadata.add_vertical_dimension(axis_1_half)
    metadata.add_vertical_dimension(axis_2_full)
    metadata.add_vertical_dimension(axis_2_half)

    with pytest.raises(InvalidMetadataError):
        validate_metadata(metadata, False)

    assert "Maximum 2 primary vertical axis (4 found)" in caplog.text


def test_validate_vert_axes_force(caplog):
    """Test error not raised when there are multiple primary axes defined and
    force is True"""

    axis_1_full = VerticalDimension("model_vert_axis_1_full_levels",
                                    name="model_vert_axis_1_full_levels",
                                    domain_top=42,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_1_half = VerticalDimension("model_vert_axis_1_half_levels",
                                    name="model_vert_axis_1_half_levels",
                                    domain_top=42,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_2_full = VerticalDimension("model_vert_axis_2_full_levels",
                                    name="model_vert_axis_2_full_levels",
                                    domain_top=43,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')
    axis_2_half = VerticalDimension("model_vert_axis_2_half_levels",
                                    name="model_vert_axis_2_half_levels",
                                    domain_top=43,
                                    extrusion_method='uniform',
                                    positive_direction='up',
                                    primary_axis='true',
                                    number_of_layers=3,
                                    units='m')

    metadata = Metadata()
    metadata.add_vertical_dimension(axis_1_full)
    metadata.add_vertical_dimension(axis_1_half)
    metadata.add_vertical_dimension(axis_2_full)
    metadata.add_vertical_dimension(axis_2_half)

    validate_metadata(metadata, True)

    assert "Maximum 2 primary vertical axis (4 found)" in caplog.text
