##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Class for storing and operating on metadata fields and their output streams.
"""
import logging

from typing import Dict

from entities import Field, FieldGroup, OutputStream, OutputStreamField, \
        VerticalDimension, Grid, NonSpatialDimension

LOGGER = logging.getLogger("diag_metadata_collection")


class Metadata:
    """Store data for fields and output streams"""

    def __init__(self):
        self._field_groups: Dict[str, FieldGroup] = {}
        self._output_streams: Dict[int, OutputStream] = {}
        self._fields: Dict[str, Field] = {}
        self._vertical_dimensions: Dict[str, VerticalDimension] = {}
        self._non_spatial_dimensions: Dict[str, NonSpatialDimension] = {}
        self._grids: Dict[str, Grid] = {}

    def add_field(self, field: Field):
        """
        Add a Field object to the collection of fields and the ID to the
        given field group
        :param field: Field object to add to metadata collection
        """
        self._fields.update({field.unique_id: field})
        self._field_groups[field.field_group_id].add_field(field)

    def get_field(self, field_id: str) -> Field:
        """
        :param field_id: Identifier for field
        :return: A Field object with the ID given
        """
        return self._fields[field_id]

    def get_fields(self) -> [Field]:
        """:return: A sorted list of fields"""
        return sorted(self._fields.values(),
                      key=lambda field: field.unique_id)

    def add_field_group(self, field_group: FieldGroup):
        """Add a FieldGroup to the collection of field groups"""
        self._field_groups.update({field_group.name: field_group})

    def get_field_groups(self) -> [FieldGroup]:
        """:return: A sorted list of fields contained in the field group"""
        return sorted(self._field_groups.values(),
                      key=lambda field_group: field_group.name)

    def add_output_stream(self, stream: OutputStream):
        """Add an OutputStream object to the collection of output streams"""
        self._output_streams.update({stream.unique_id: stream})

    def get_output_streams(self) -> [OutputStream]:
        """:return: A sorted list of all output streams"""
        return sorted(self._output_streams.values(),
                      key=lambda streams: streams.unique_id)

    def add_output_stream_field(self, output_stream_field: OutputStreamField,
                                stream_unique_id: int):
        """
        Add an output stream field to the output stream corresponding to the
        given output stream ID
        :param output_stream_field: Output stream field object to add to stream
        :param stream_unique_id: Integer identifying stream to add field to
        """
        self._output_streams[stream_unique_id].add_field(output_stream_field)

    def add_vertical_dimension(self, vertical_dimension: VerticalDimension) \
            -> str:
        """Add a vertical dimension if it is not a duplicate of one already
        stored
        :param vertical_dimension: Vertical dimension object to add to metadata
        :return: Unique id of vertical dimension after de-duplication"""
        vertical_dimension.validate()
        vertical_dimension.format_level_definition()
        is_duplicate = False
        unique_id = vertical_dimension.unique_id
        for dimension in self._vertical_dimensions.values():
            if vertical_dimension == dimension:
                is_duplicate = True
                unique_id = dimension.unique_id
                break
        if not is_duplicate:
            self._vertical_dimensions.update(
                {vertical_dimension.unique_id: vertical_dimension}
            )
        return unique_id

    def get_vertical_dimensions(self) -> [VerticalDimension]:
        """:return: A sorted list of vertical dimensions"""
        return sorted(self._vertical_dimensions.values(),
                      key=lambda vertical_dim: vertical_dim.unique_id)

    def activate_vertical_dimension(self, vert_dim_id: str) -> None:
        """
        Set 'active' attribute of vertical dimension to True
        :param vert_dim_id: Identifier for vertical_dimension
        """
        if vert_dim_id in self._vertical_dimensions:
            vert_dim = self._vertical_dimensions[vert_dim_id]
            vert_dim.active = True
        else:
            raise KeyError(f"Vertical dimension '{vert_dim_id}' not found")

    def count_fixed_vert_axes(self) -> int:
        """:return: Number of fixed vertical axes in the metadata collection"""
        return sum(axis.startswith("fixed_vert_axis_")
                   for axis in self._vertical_dimensions)

    def add_grid(self, grid: Grid) -> str:
        """Add a grid if it is not a duplicate of one already stored
        :param grid: Grid object to add to metadata collection
        :return: Unique id of grid after de-duplication"""
        is_duplicate = False
        unique_id = grid.unique_id
        for stored_grid in self._grids.values():
            if grid == stored_grid:
                is_duplicate = True
                unique_id = stored_grid.unique_id
                break
        if not is_duplicate:
            self._grids.update({grid.unique_id: grid})
        return unique_id

    def get_grids(self) -> [Grid]:
        """:return: A sorted list of grids"""
        return sorted(self._grids.values(), key=lambda grid: grid.unique_id)

    def get_non_spatial_dimensions(self) -> [NonSpatialDimension]:
        """:return: A list of non-spatial dimensions"""
        return self._non_spatial_dimensions

    def get_non_spatial_dimension(self, key) -> NonSpatialDimension:
        """:return: A specific non-spatial dimension"""
        if key not in self._non_spatial_dimensions:
            LOGGER.error("""Non-spatial dimension not recognised.
                         Metadata mismatch: %s""", key)
        return self._non_spatial_dimensions[key]

    def add_non_spatial_dimension(self, non_spatial_dimension):
        """Add a non-spatial dimension if it is not a duplicate of one already
        stored.
        :param non_spatial_dimension: non_spatial dimension object to add to
        metadata"""
        is_duplicate = False
        for dimension in self._non_spatial_dimensions.values():
            if non_spatial_dimension == dimension:
                is_duplicate = True
                break
        if not is_duplicate:
            self._non_spatial_dimensions.update(
                {non_spatial_dimension.name: non_spatial_dimension}
            )
