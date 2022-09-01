##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Output metadata from a Metadata object to an XIOS iodef.xml file"""

from pathlib import Path
from logging import getLogger
from xml.etree import ElementTree
from xml.dom import minidom

from entities import Field, FieldGroup, OutputStream, VerticalDimension, Grid,\
    NonSpatialDimension
from diagnostics_metadata_collection import Metadata
from extrusion_namelist_updater import update_extrusion_namelist

FIELD_NODE_VARIABLE_TYPES = {
    'mesh_id': 'int', 'function_space': 'string',
    'element_order': 'int', 'io_driver': 'string',
    'field_type': 'string', 'checksum': 'bool'
}

# mappings for attribute name in iodef : name in Metadata
FIELD_NODE_VARIABLES = {
    'mesh_id': 'mesh_id',
    'function_space': 'function_space',
    'element_order': 'order',
    'io_driver': 'io_driver',
    'field_type': 'data_type',
    'checksum': 'checksum'
}
FIELD_NODE_MANDATORY_ATTRS = {
    'id': 'unique_id',
    'unit': 'units',
    'grid_ref': 'grid_ref'
}
FIELD_NODE_OPTIONAL_ATTRS = {
    'long_name': 'long_name',
    'standard_name': 'standard_name',
    'freq_op': 'timestep'
}
FIELD_GROUP_NODE_OPTIONAL_ATTRS = {
    'freq_op': 'timestep',
    'operation': 'temporal'
}
FILE_NODE_ATTRS = {
    'id': 'name',
    'name': 'name',
    'output_freq': 'timestep',
    'convention': 'convention'
}

VERT_AXIS_NODE_ATTRS = {
    'id': 'unique_id',
    'n_glo': 'size',
    'name': 'name',
    'positive': 'positive_direction',
    'unit': 'units',
    'value': 'level_definition'
}

NON_SPATIAL_AXIS_ATTRS = {
    'id': 'name',
    'n_glo': 'size',
    'name': 'name',
    'value': 'definition',
    'unit': 'unit',
}

LOGGER = getLogger("reconfigurator.presenter")


class CommentedTreeBuilder(ElementTree.TreeBuilder):
    """XML tree builder to enable parser to handle comments"""
    def comment(self, data):
        self.start(ElementTree.Comment, {})
        self.data(data)
        self.end(ElementTree.Comment)


class MetadataIodefGenerator:
    """
    Takes a Metadata object and writes the metadata to XML in the form of an
    XIOS configuration file (i.e an iodef.xml file)
    """

    def __init__(self, metadata: Metadata):
        LOGGER.debug("Initialising iodef presenter")
        self._metadata = metadata
        self._root_node = None
        self._context_node = None
        self._axis_def = None
        self._grid_def = None
        self._field_def = None
        self._file_def = None

    def _create_root_nodes(self):
        """
        Create a root <context> element with <field_definition>,
        <axis_definition>, <grid_definition> and <file_definition> elements
        to append metadata to
        """
        LOGGER.info("Creating root nodes")
        self._root_node = ElementTree.Element('simulation')
        self._context_node = ElementTree.SubElement(self._root_node,
                                                    'context',
                                                    id="diagnostics")
        self._axis_def = ElementTree.SubElement(self._context_node,
                                                'axis_definition')
        self._grid_def = ElementTree.SubElement(self._context_node,
                                                'grid_definition',
                                                prec="8")
        self._field_def = ElementTree.SubElement(self._context_node,
                                                 'field_definition',
                                                 prec="8")
        self._file_def = ElementTree.SubElement(self._context_node,
                                                'file_definition')
        self._file_def.set('par_access', "collective")
        self._file_def.set('time_counter', "none")
        self._file_def.set('type', "one_file")

    def _find_template_nodes(self, template_path: str):
        """
        Use ElementTree to get the root, axis_definition, grid_definition
        field_definition, and file_definition nodes in the given template XML
        file
        :param template_path: String path to template XML file
        """
        parser = ElementTree.XMLParser(target=CommentedTreeBuilder())
        xml_tree = ElementTree.parse(template_path, parser)

        LOGGER.info("Template iodef specified, finding template nodes")
        template_invalid = False

        self._root_node = xml_tree.getroot()
        context_node = xml_tree.getroot().find('context')

        self._axis_def = context_node.find('axis_definition')
        if self._axis_def is None:
            template_invalid = True
            LOGGER.error('No axis definition found in template iodef')

        self._grid_def = context_node.find('grid_definition')
        if self._grid_def is None:
            template_invalid = True
            LOGGER.error('No grid definition found in template iodef')

        self._field_def = context_node.find('field_definition')
        if self._field_def is None:
            template_invalid = True
            LOGGER.error('No field definition found in template iodef')

        self._file_def = context_node.find('file_definition')
        if self._file_def is None:
            template_invalid = True
            LOGGER.error('No file definition found in template iodef')

        if template_invalid:
            raise ValueError("Definition section(s) missing from iodef "
                             "template")

    def _create_vertical_axis_node(self, axis: VerticalDimension) \
            -> ElementTree.Element:
        """
        Create and populate <axis> elements within the axis_definition node

        :param axis: VerticalDimension axis to be created
        :return: The newly made node in the XML tree containing the axis
        """
        LOGGER.debug("Creating axis: %s", axis.name)
        axis_node = ElementTree.SubElement(self._axis_def, 'axis')

        # Loop through axis attributes and set them
        for iodef_attr, entity_attr in sorted(VERT_AXIS_NODE_ATTRS.items()):
            attribute_value = getattr(axis, entity_attr, None)
            axis_node.set(iodef_attr, str(attribute_value))

        return axis_node

    def _create_non_spatial_axis_node(self, axis: NonSpatialDimension) \
            -> ElementTree.Element:
        """
        Create and populate <axis> elements within the axis_definition node

        :param axis: NonSpatialDimension axis to be created
        :return: The newly made node in the XML tree containing the axis
        """
        LOGGER.debug("Creating axis: %s", axis.name)
        axis_node = ElementTree.SubElement(self._axis_def, 'axis')

        # Loop through axis attributes and set them
        for iodef_attr, entity_attr in sorted(NON_SPATIAL_AXIS_ATTRS.items()):
            attribute_value = getattr(axis, entity_attr, None)
            axis_node.set(iodef_attr, str(attribute_value))

        if not axis.is_numerical:
            axis_node.attrib["label"] = axis_node.attrib.pop("value")
        return axis_node

    def _create_grid_node(self, grid: Grid) -> ElementTree.Element:
        """
        Create and populate <grid> elements within the grid_definition node

        :param grid: Grid to be created
        :return: The newly made node in the XML tree containing the grid
        """
        LOGGER.debug("Creating grid: %s", grid.unique_id)
        grid_node = ElementTree.SubElement(self._grid_def, 'grid')
        grid_node.set("id", grid.unique_id)

        domain_node = ElementTree.SubElement(grid_node, 'domain')
        domain_node.set("domain_ref", grid.domain)

        for axis in grid.axes:
            axis_node = ElementTree.SubElement(grid_node, 'axis')
            axis_node.set("axis_ref", axis)
        return grid_node

    def _create_field_group_node(self, field_group: FieldGroup) \
            -> ElementTree.Element:
        """
        Create and populate <field_group> elements within the
        field_definition node

        :param field_group: FieldGroup field group to be created
        :return: The newly made node in the XML tree containing the field group
        """
        LOGGER.debug("Creating field group: %s", field_group.name)
        field_group_node = ElementTree.SubElement(self._field_def,
                                                  'field_group')
        field_group_node.set('enabled', '.TRUE.')
        field_group_node.set('id', field_group.name)

        # Loop through optional field group attributes to set them
        for iodef_attr, entity_attr in sorted(
                FIELD_GROUP_NODE_OPTIONAL_ATTRS.items()):
            attribute_value = getattr(field_group, entity_attr, None)
            if attribute_value:
                field_group_node.set(iodef_attr, attribute_value)

        return field_group_node

    @staticmethod
    def _create_field_node(field: Field, group_node: ElementTree.Element) \
            -> ElementTree.Element:
        """
        Create <field> nodes in the XML tree inside the given field_group node

        :param field: Field entity containing metadata for the new field
        :param group_node: <field_group> node in the XML tree to append the
                            <field> node to
        :return: Newly created <field> node
        """
        LOGGER.debug("Creating field: %s", field.unique_id)
        field_node = ElementTree.SubElement(group_node, 'field')

        # Set mandatory field node attributes
        for iodef_attr, entity_attr in sorted(
                FIELD_NODE_MANDATORY_ATTRS.items()):
            field_node.set(iodef_attr, getattr(field, entity_attr))

        # Set optional field node attributes
        for iodef_attr, entity_attr in sorted(
                FIELD_NODE_OPTIONAL_ATTRS.items()):
            attribute_value = getattr(field, entity_attr, None)
            if attribute_value:
                field_node.set(iodef_attr, attribute_value)

        # Add variables to field node
        # These are used to store additional information about the field
        # in the iodef file, to be used within LFRic
        for variable_name, entity_attr in sorted(FIELD_NODE_VARIABLES.items()):
            attribute_value = getattr(field, entity_attr, None)
            if attribute_value is not None:
                # FoX requires booleans to be lowercase when parsing them
                if isinstance(attribute_value, bool):
                    attribute_value = str(attribute_value).lower()
                ElementTree.SubElement(
                    field_node,
                    'variable',
                    name=variable_name,
                    type=FIELD_NODE_VARIABLE_TYPES[variable_name]
                ).text = str(attribute_value)

        return field_node

    def _create_file_node(self, output_stream: OutputStream) \
            -> ElementTree.Element:
        """
        Create a <file> node within the file_definition

        :param output_stream: Information about the output stream to be put
                              into the <file> node
        :return: Newly created <file> node within the file_definition
        """
        LOGGER.debug("Creating file: %s", output_stream.name)
        file_node = ElementTree.SubElement(self._file_def, 'file',
                                           convention='CF', enabled='.TRUE.')

        for iodef_attr, stream_attr in sorted(FILE_NODE_ATTRS.items()):
            file_node.set(iodef_attr, getattr(output_stream, stream_attr))

        # create the field node within this file node
        for stream_field in output_stream.get_fields():
            file_field_name = f"{stream_field.field_ref}__" \
                              f"{stream_field.temporal}_" \
                              f"{output_stream.timestep}"
            ElementTree.SubElement(file_node, 'field',
                                   field_ref=stream_field.field_ref,
                                   name=file_field_name,
                                   operation=stream_field.temporal)
            # can add extra variables to the stream field node here
        return file_node

    def write_tree_to_xml(self, output_path: str):
        """
        Writes tree to XML

        :param output_path: Path to output XML
        """
        LOGGER.info("Writing XML file")

        xml = minidom.parseString(ElementTree.tostring(self._root_node))
        xml_text = '\n'.join([line for line in xml.toprettyxml(
            indent=' ' * 4).split('\n') if line.strip()])

        with open(output_path, 'w') as output_file:
            output_file.write(xml_text)

    def generate(self, output_file_path: str, template_file_path: str = None,
                 namelist_file_path: str = None):
        """
        Write iodef.xml file using data from the Metadata object. This file
        is based on an input template iodef_xml file if one is supplied

        :param output_file_path: File path to output iodef.xml to
        :param template_file_path: Template iodef.xml file to build XML tree
                                   from
        :param namelist_file_path: LFRic configuration file containing
                                   extrusion namelist to be updated

        """
        if template_file_path is None:
            LOGGER.info("No template XML file specified")
            self._create_root_nodes()
        else:
            self._find_template_nodes(template_file_path)

        # Build the XML tree for the axis_definition
        for axis in self._metadata.get_vertical_dimensions():
            if axis.active:
                self._create_vertical_axis_node(axis)

                # Update namelist if this is the primary vertical axis
                if namelist_file_path and axis.primary_axis == 'true':
                    LOGGER.info("Updating %s", namelist_file_path)
                    update_extrusion_namelist(namelist_file_path, axis)
        axes = len(self._axis_def.findall("axis"))
        LOGGER.info("Created %s axes", axes)

        for axis in self._metadata.get_non_spatial_dimensions().values():
            self._create_non_spatial_axis_node(axis)
            num_nsd = len(self._metadata.get_non_spatial_dimensions())
            LOGGER.info("Created %s Immutable Non-spatial Dimensions", num_nsd)

        # Build the XML tree for the grid_definition
        for grid in self._metadata.get_grids():
            self._create_grid_node(grid)
        grids = len(self._grid_def.findall("grid"))
        LOGGER.info("Created %s grids", grids)

        # Build the XML tree for the field_definition
        fields = 0
        for field_group in self._metadata.get_field_groups():
            group_node = self._create_field_group_node(field_group)
            for field_id in field_group.get_fields():
                field = self._metadata.get_field(field_id)
                if field.active:
                    self._create_field_node(field, group_node)
            fields += len(group_node.findall("field"))
        groups = len(self._field_def.findall("field_group"))
        LOGGER.info("Created %s fields in %s groups", fields, groups)

        # Build the XML tree for the file_definition
        stream_fields = 0
        for stream in self._metadata.get_output_streams():
            file = self._create_file_node(stream)
            stream_fields += len(file.findall("field"))
        files = self._file_def.findall("file")
        LOGGER.info("Created %s fields in %s files", stream_fields, len(files))

        output_file_path = Path(output_file_path).expanduser()
        self.write_tree_to_xml(output_file_path)
