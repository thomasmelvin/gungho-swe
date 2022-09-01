##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Populates a Metadata object from data in a rose suite and a JSON file
containing immutable metadata
"""

import copy
import hashlib
import json
import re
import logging
from pathlib import Path
from typing import TextIO, Union

from entities import Field, FieldGroup, OutputStream, OutputStreamField, \
    VerticalDimension, Grid, NonSpatialDimension
from diagnostics_metadata_collection import Metadata

KEY_VALUE_RE = re.compile(r'^([a-zA-Z_]+)=\'?([\w.]+)\'?$')
KEY_VALUE_LIST_RE = re.compile(
    r'^(\w+)=(\'?[\w.]+\'?(?:, ?\'?[\w.]+\'?)*)(?!,)')

FIELD_CONFIG_RE = re.compile(r'^\[field_config:([a-zA-Z_]+):([a-zA-Z_]+)\]$')
FIELD_RE = re.compile(
    r'^([a-zA-Z]+(?:_[a-zA-Z]+)*__[a-zA-Z]+(?:_[a-zA-Z]+)*)=(true|false)$')
ADDITIONAL_INPUT_RE = re.compile(
    r'^[a-zA-Z_]+__[a-zA-Z_]+__([a-zA-Z]+)=(true|false)$')
BASE_OUTPUT_RE = re.compile(r'^\[output_stream\(([0-9]+)\)\]$')
OUTPUT_FIELD_RE = re.compile(
    r'^\[output_stream\(([0-9]+)\):field\(([0-9]+)\)\]$')
VERTICAL_DIMENSION_RE = re.compile(r'^\[vertical_dimension\(([0-9]+)\)\]$')
NON_SPATIAL_DIMENSION_RE = re.compile(r'^\[non_spatial_dimensions]$')

REQUIRED_METADATA = ['unique_id', 'long_name', 'standard_name',
                     'units', 'grid_ref', 'function_space', 'mesh_id', 'order',
                     'io_driver', 'data_type']
LOGGER = logging.getLogger("reconfigurator.extractor")


def create_md5_checksum(obj) -> str:
    """
    :param obj: Object to hash
    :return: String containing checksum
    """
    obj_str = json.dumps(obj, sort_keys=True)
    checksum_hex = hashlib.md5(obj_str.encode('utf-8')).hexdigest()
    return f'md5: {checksum_hex}'


class MetadataExtractor:
    """
    Get fields and output streams which are configured in a Rose suite and
    extract their immutable metadata from a JSON file
        Usage:
        >>> extractor = MetadataExtractor('path/to/suite/', 'immutable.json')
        >>> extractor.extract_metadata()
    """

    def __init__(self, diagnostic_config_path: str,
                 immutable_data_file: Union[str, 'Path']):
        LOGGER.debug("Initialising metadata object")
        self._metadata = Metadata()
        immutable_data_path = Path(immutable_data_file).expanduser()
        self._immutable_metadata = self._get_immutable_data(
            immutable_data_path)

        self._rose_app_path = Path(diagnostic_config_path).expanduser()

        if not self._rose_app_path.exists():
            raise IOError("Could not find rose-app.conf")

    @staticmethod
    def _get_immutable_data(file_path) -> dict:
        """
        Returns immutable metadata from the JSON file given.
        Validates the MD5 checksum to ensure the file has not been modified by
        hand after it was generated.

        :param file_path: Path to the immutable data JSON file
        :return: Multidimensional dictionary containing immutable metadata
        """
        with open(file_path, 'r') as file_handle:
            LOGGER.debug("Opening immutable metadata")
            file_data = json.load(file_handle)
            if 'checksum' not in file_data:
                raise KeyError("Can't find checksum "
                               "to validate immutable data")

            # Generate a checksum from a copy of the data without the checksum
            # and check it against the original
            immutable_metadata = copy.deepcopy(file_data)
            del immutable_metadata['checksum']
            checksum = create_md5_checksum(immutable_metadata)
            if file_data['checksum'] != checksum:
                raise RuntimeError("Immutable data has been modified by hand\n"
                                   f"Expected checksum: "
                                   f"{file_data['checksum']}")
            LOGGER.info("Immutable metadata checksum valid")
            return immutable_metadata

    def _add_immutable_field_metadata(self, field: Field) -> Field:
        """
        Populates the given metadata field with required metadata attributes
        taken from the immutable metadata

        :param field: Object containing field's data
        """
        # Get the current field's immutable metadata
        section_name, group_name = field.field_group_id.split('__')
        metadata_section = self._immutable_metadata["meta_data"]["sections"][
            section_name]["groups"][group_name]["fields"][field.unique_id]

        # Add the immutable metadata to the field
        for attr in metadata_section:
            if attr in REQUIRED_METADATA:
                setattr(field, attr, metadata_section[attr])

        # Add fixed vertical dimension if field has one
        vert_dim_meta = metadata_section.get("vertical_dimension", None)
        if vert_dim_meta is not None:
            # Create immutable vertical dimension
            if vert_dim_meta.get("level_definition") is not None:
                positive = {"POSITIVE_UP": "up", "POSITIVE_DOWN": "down"}
                id_number = self._metadata.count_fixed_vert_axes() + 1
                vertical_dimension = VerticalDimension(
                    "fixed_vert_axis_" + str(id_number),
                    name="fixed_vert_axis_" + str(id_number),
                    positive_direction=positive[vert_dim_meta["positive"]],
                    primary_axis='false',
                    level_definition=vert_dim_meta["level_definition"],
                    number_of_layers=len(
                        vert_dim_meta["level_definition"]),
                    units=vert_dim_meta["units"]
                )

                vert_dim_id = self._metadata.add_vertical_dimension(
                    vertical_dimension)

                field.vertical_dimension_id = vert_dim_id

            else:
                # Add information about top_arg and bottom_arg in future ticket
                field.vertical_dimension_id = "mutable"

        immutable_non_spatial_dimension_meta = metadata_section.get(
            "non_spatial_dimension", None)
        for value in immutable_non_spatial_dimension_meta.values():
            unit = None
            definition = None
            if 'label_definition' in value and 'axis_definition' in value:
                LOGGER.error("Non-spatial dimension %s has both a "
                             "'label_definition' and an 'axis_definition'."
                             " This is invalid.", value["name"])
                raise ValueError("Non-spatial dimension %s has both a "
                                 "'label_definition' and an 'axis_definition.'"
                                 " This is invalid." % value["name"])

            if 'label_definition' in value:
                definition = value['label_definition']
            if 'axis_definition' in value:
                definition = value['axis_definition']

            if 'unit' in value.keys():
                unit = value['unit']

            nsd = NonSpatialDimension(name=value['name'],
                                      definition=definition,
                                      unit=unit)

            field.non_spatial_dimension.append(nsd)
            self._metadata.add_non_spatial_dimension(nsd)

        # Set domain reference based on function space
        if field.function_space in ["W0"]:
            field.domain_ref = "node"
        elif field.function_space in ["W2H", "W2"]:
            field.domain_ref = "edge"
        elif field.function_space.upper() in ["W3", "WTHETA"]:
            field.domain_ref = "face"
        else:
            raise ValueError(f"Invalid function space "
                             f"{field.function_space}")

        # Hard-code following attributes for now
        field.mesh_id = 1
        field.io_driver = 'WRITE_FIELD_FACE'
        return field

    def _parse_field_config(self, file_pointer: TextIO, section_name: str,
                            group_name: str) -> str:
        """
        Parse the [field_config:...] section of the rose-app the file pointer
        is pointing to and set the active status and checksum status of each
        metadata field according to the values found

        :param file_pointer: IO object for the rose-app.conf file
        :param section_name: String containing the name the metadata
                               section currently being processed
        :param group_name: String containing the name of the metadata group
                            currently being processed
        :return: The line the file pointer is currently pointing to so
                 calling method can continue parsing where this method left off
        """
        field_group_id = section_name + "__" + group_name
        field_group = FieldGroup(field_group_id)
        self._metadata.add_field_group(field_group)

        # Process fields within group
        line = file_pointer.readline()

        # Stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            field_match = FIELD_RE.search(line)
            key_value_match = KEY_VALUE_RE.search(line)

            # Matches field__unique_id=(true|false)
            # Enables field and looks for any additional configuration options
            if field_match is not None:
                field_id, active = field_match.groups()
                field = Field(field_id, field_group_id,
                              active=active.lower() == 'true')

                # Process additional configuration for field
                # currently only checksum is supported
                line = file_pointer.readline()
                match = ADDITIONAL_INPUT_RE.search(line)

                while line and line != '[' and match:
                    input_name, value = match.groups()

                    if input_name == 'checksum':
                        field.checksum = value.lower() == 'true'
                    else:
                        raise ValueError("Unrecognised additional field input "
                                         f"'{input_name}' found for field "
                                         f"'{field_id}'")
                    line = file_pointer.readline()
                    match = ADDITIONAL_INPUT_RE.search(line)

                field = self._add_immutable_field_metadata(field)
                self._metadata.add_field(field)

            # Matches key=value
            # Used to find vertical dimension id for fields in this group
            elif key_value_match is not None:
                key, value = key_value_match.groups()

                if key == "vertical_dimension_for_group":
                    setattr(field_group, "vertical_dimension_id",
                            "model_vert_axis_" + value)
                line = file_pointer.readline()

            else:
                line = file_pointer.readline()

        # Add group vertical dimension to fields on mutable vertical axis
        for field_id in field_group.get_fields():
            field = self._metadata.get_field(field_id)

            if field.vertical_dimension_id == "mutable":
                if not field_group.vertical_dimension_id:
                    raise Exception(f"Field group '{field_group.name}' "
                                    f"missing vertical dimension")

                if field.function_space in ['W0', 'W2H', 'Wtheta']:
                    field.vertical_dimension_id = \
                        field_group.vertical_dimension_id + '_full_levels'
                elif field.function_space in ['W3']:
                    field.vertical_dimension_id = \
                        field_group.vertical_dimension_id + '_half_levels'
                else:
                    raise ValueError(f"Invalid function space "
                                     f"{field.function_space}")

            axes = []
            if field.vertical_dimension_id is not None:
                axes.append(field.vertical_dimension_id)

            if field.non_spatial_dimension:
                for nsd in field.non_spatial_dimension:
                    axes.append(nsd.name)

            grid = Grid(field.domain_ref, axes)
            grid_id = self._metadata.add_grid(grid)
            field.grid_ref = grid_id

        return line

    def _parse_output_stream(self, file_pointer: TextIO, stream_id: int) \
            -> str:
        """
        Read through an [output_stream(x)] section of the rose-app and store
        information about the output_stream being declared

        :param file_pointer: IO object pointing to the output_stream(x) section
        :param stream_id: Identifier for current output_stream being processed
        :return: Current line of the file pointer to allow calling method to
                  continue parsing where this method left off
        """
        line = file_pointer.readline()
        stream = OutputStream(stream_id)

        # Stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            # Check if line has key-value pair
            match = KEY_VALUE_RE.search(line)

            if match is not None:
                key, value = match.groups()
                recognised_keys = {'name': 'name', 'timestep': 'timestep'}

                if key in recognised_keys:
                    setattr(stream, recognised_keys[key], value)
                else:
                    LOGGER.warning("Key '%s' unrecognised in "
                                   "[output_stream(%s)]", key, stream_id)
                    LOGGER.debug("Recognised keys are: %s",
                                 ', '.join(recognised_keys.keys()))
            line = file_pointer.readline()
        self._metadata.add_output_stream(stream)
        return line

    def _parse_output_stream_field(self, file_pointer: TextIO, stream_id: int,
                                   output_stream_field_id: int) -> str:
        """
        Read an output_stream(x):field(y) section (in the rose-app) and store
        info about the output stream field in the Metadata object

        :param file_pointer: IO object pointing to the start of an
                              output_stream(x):field(y) section
        :param stream_id: Identifier for current output_stream being processed
        :param output_stream_field_id: Identifier for output_stream field
                                being processed
        :return: Current line of the file pointer to allow calling method to
                  continue parsing where this method left off
        """
        output_stream_field = OutputStreamField(output_stream_field_id)
        line = file_pointer.readline()

        # Stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            # Check if line has key-value pair
            match = KEY_VALUE_RE.search(line)

            if match is not None:
                key, value = match.groups()
                recognised_keys = {'id': 'field_ref', 'temporal': 'temporal'}

                if key in recognised_keys:
                    setattr(output_stream_field, recognised_keys[key], value)
                else:
                    LOGGER.warning("Key '%s' unrecognised in "
                                   "[output_stream(%s):field(%s)]",
                                   key, stream_id, output_stream_field_id)
                    LOGGER.debug("Recognised keys are: %s",
                                 ', '.join(recognised_keys.keys()))
            line = file_pointer.readline()

        # Check whether field has been configured
        try:
            self._metadata.get_field(output_stream_field.field_ref)
        except KeyError:
            LOGGER.warning("Output stream field '%s' has no field config",
                           output_stream_field.field_ref)

        self._metadata.add_output_stream_field(output_stream_field, stream_id)
        return line

    def _parse_non_spatial_dimension(self, file_pointer: TextIO) -> str:

        """
        Read a non-spatial dimension section (in the rose-app) and store
        info about the non-spatial dimension in the Metadata object.
        :param file_pointer: IO object pointing to the start of an
                             output_stream(x):field(y) section
        :return: Current line of the file pointer to allow calling method to
                 continue parsing where this method left off
        """

        line = file_pointer.readline()

        while line and line[0] != '[':

            # Check if line has key-list pair
            match = KEY_VALUE_LIST_RE.search(line)

            if match:
                key, values = match.groups()
                values = list(values.split(","))

                # Add any immutable meta data from the JSON file
                # Add unit of measure
                if self._immutable_metadata['meta_data'][
                        'non_spatial_dimensions'][key].get('unit', False):
                    unit = self._immutable_metadata['meta_data'][
                        'non_spatial_dimensions'][key]['unit']
                    self._metadata.get_non_spatial_dimension(key).add_unit(
                        unit)
                self._metadata.get_non_spatial_dimension(key).add_definition(
                    values)

            line = file_pointer.readline()

        return line

    def _parse_vertical_dimension(self, file_pointer: TextIO,
                                  dimension_id: int) -> str:
        """
        Read a vertical_dimension(x) section (in the rose-app) and store
        info about the vertical dimension in the Metadata object. Creates one
        vertical axis object for fields using 'half_levels' and one for
        'full_levels'

        :param file_pointer: IO object pointing to the start of an
                              output_stream(x):field(y) section
        :param dimension_id: Identifier for current dimension being processed
        :return: Current line of the file pointer to allow calling method to
                 continue parsing where this method left off
        """
        vert_dim_half = VerticalDimension(
            f"model_vert_axis_{dimension_id}_half_levels")
        vert_dim_full = VerticalDimension(
            f"model_vert_axis_{dimension_id}_full_levels")
        line = file_pointer.readline()

        # Stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            # Check if line has key-value pair
            single_value_match = KEY_VALUE_RE.search(line)
            # Check if line has key-list pair
            list_match = KEY_VALUE_LIST_RE.search(line)

            if single_value_match is not None:
                key, value = single_value_match.groups()
                recognised_keys = {'domain_top': 'domain_top',
                                   'extrusion_method': 'extrusion_method',
                                   'number_of_layers': 'number_of_layers',
                                   'positive': 'positive_direction',
                                   'primary_axis': 'primary_axis',
                                   'units': 'units'}
                if key in recognised_keys:
                    setattr(vert_dim_half, recognised_keys[key], value)
                    setattr(vert_dim_full, recognised_keys[key], value)
                elif key == 'name':
                    setattr(vert_dim_half, key, value + '_half_levels')
                    setattr(vert_dim_full, key, value + '_full_levels')

            elif list_match is not None:
                key, values = list_match.groups()
                values = [float(i) for i in values.split(",")]

                if key == "level_definition":
                    setattr(vert_dim_half, "level_definition", values)
                    setattr(vert_dim_full, "level_definition", values)
                else:
                    LOGGER.warning("Key '%s' unrecognised as a list in "
                                   "[vertical_dimension(%s)]",
                                   key, dimension_id)

            line = file_pointer.readline()

        self._metadata.add_vertical_dimension(vert_dim_half)
        self._metadata.add_vertical_dimension(vert_dim_full)

        return line

    def _parse_rose_app(self):
        """
        Read through the rose-app.conf to find which fields are present in the
        current Rose suite configuration, store their data, and store
        information about which output streams they are destined for
        """
        with open(self._rose_app_path, 'r') as file_pointer:
            LOGGER.debug("Opening rose-app.conf")
            line = file_pointer.readline()

            # File pointer is advanced by each individual method as it
            # processes its own section (updated line is returned by each
            # method pointing to start of next section)
            while line:
                # If line is start of a new section
                if line[0] == '[':
                    LOGGER.debug("Parsing section: %s", line[:-1])

                    # Matches [field_config:(section):(field_group)]
                    match = FIELD_CONFIG_RE.search(line)
                    if match is not None:
                        line = self._parse_field_config(
                            file_pointer,
                            match.group(1),
                            match.group(2)
                        )
                        continue

                    # Matches [output_stream(x)]
                    match = BASE_OUTPUT_RE.search(line)
                    if match is not None:
                        line = self._parse_output_stream(
                            file_pointer,
                            int(match.group(1))
                        )
                        continue

                    # Matches [output_stream(x):field(y)]
                    match = OUTPUT_FIELD_RE.search(line)
                    if match is not None:
                        line = self._parse_output_stream_field(
                            file_pointer,
                            int(match.group(1)),
                            int(match.group(2))
                        )
                        continue

                    # Matches [vertical_dimension(x)]
                    match = VERTICAL_DIMENSION_RE.search(line)
                    if match is not None:
                        line = self._parse_vertical_dimension(
                            file_pointer,
                            int(match.group(1)))
                        continue

                    # Matches [non_spatial_dimension]
                    match = NON_SPATIAL_DIMENSION_RE.search(line)
                    if match:
                        line = self._parse_non_spatial_dimension(file_pointer)
                        continue

                line = file_pointer.readline()

            # Log number of fields in configs
            LOGGER.info("Found %s fields in %s field groups",
                        len(self._metadata.get_fields()),
                        len(self._metadata.get_field_groups()))

            # Log number of fields in output streams
            stream_fields = 0
            for stream in self._metadata.get_output_streams():
                stream_fields += len(stream.get_fields())
            LOGGER.info("Found %s fields in %s output streams", stream_fields,
                        len(self._metadata.get_output_streams()))

    def extract_metadata(self):
        """
        Populate a Metadata object with data extracted from the rose-app file
        and the immutable data

        :return: A newly populated Metadata object
        """
        self._parse_rose_app()

        return self._metadata
