##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
"""
Validation step for diagnostic fields

Checks that all fields have either a long or a standard name and that their
vertical dimension is in the metadata collection of valid vertical dimensions
"""

from logging import getLogger, WARNING, ERROR

from diagnostics_metadata_collection import Metadata

LOGGER = getLogger("reconfigurator.validator")


class InvalidMetadataError(ValueError):
    """Error showing field(s) have invalid metadata"""


def validate_metadata(metadata: Metadata, force: bool = False) -> None:
    """Validator that reports any invalid metadata"""

    invalid_vert_axes = validate_vert_axes(metadata, force)
    invalid_fields = validate_fields(metadata, force)

    if (invalid_vert_axes or invalid_fields) and not force:
        raise InvalidMetadataError("Metadata invalid")

    LOGGER.info("Metadata valid")


def validate_fields(metadata: Metadata, force: bool = False) -> bool:
    """Validator that reports any diagnostic field(s) with invalid metadata"""
    invalid_fields = False
    valid_vert_dims = [vert_dim.unique_id
                       for vert_dim in metadata.get_vertical_dimensions()]
    LOGGER.info("Checking if fields have long or standard name")

    for field in metadata.get_fields():
        LOGGER.debug("Checking %s", field.unique_id)

        # Check that field has long or standard name
        if not (field.long_name or field.standard_name):
            LOGGER.log(WARNING if force else ERROR,
                       "Field '%s' missing long and standard names",
                       field.unique_id)
            invalid_fields = True

        # Check that vertical dimension is in vertical dimension collection
        if field.vertical_dimension_id:
            if field.vertical_dimension_id not in valid_vert_dims:
                LOGGER.log(WARNING if force else ERROR,
                           "Vertical dimension '%s' for field '%s' not "
                           "defined",
                           field.vertical_dimension_id, field.unique_id)
                invalid_fields = True
            # If dimension valid then set active to True so it will be output
            else:
                metadata.activate_vertical_dimension(
                    field.vertical_dimension_id)

    return invalid_fields


def validate_vert_axes(metadata: Metadata, force: bool = False) -> bool:
    """Validator that checks there aren't too many primary vertical axes"""
    LOGGER.info("Checking number of primary axes")

    n_primary_axes = sum(vert_dim.primary_axis == 'true' for vert_dim in
                         metadata.get_vertical_dimensions())

    # 2 due to there being one each for full and half levels
    invalid_vert_axes = n_primary_axes > 2

    if invalid_vert_axes:
        LOGGER.log(WARNING if force else ERROR,
                   "Maximum 2 primary vertical axis (%i found)",
                   n_primary_axes)

    return invalid_vert_axes
