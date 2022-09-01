##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Functionality to read LFRic configuration.nml and update an existing valid
extrusion namelist with values parsed from Rose configuration"""
import re

from entities import VerticalDimension

EXTRUSION_REGEX = re.compile(
        r"&extrusion\s+"
        r"(?P<domain_top>domain_top\s*=\s*\d+\.?\d*),?\s+"
        r"(?P<method>method\s*=\s*'[\d\w]+'),?\s+"
        r"(?P<number_of_layers>number_of_layers\s*=\s*\d+),?\s+"
        r"/\s+"
)


def update_extrusion_namelist(namelist_path: str,
                              vertical_dimension: VerticalDimension) -> None:
    """Read in LFRic configuration.nml and update extrusion namelist with
    the given new values
    :param namelist_path: Path to LFRic configuration.nml. This should contain
                          a valid extrusion namelist
    :param vertical_dimension: Vertical dimension object containing domain_top,
                               extrusion_method and number_of_layers
                               attributes"""
    with open(namelist_path, 'r') as file:
        namelist = file.read()

    match = EXTRUSION_REGEX.search(namelist)
    if not match:
        raise ValueError(f"Can't find valid extrusion namelist in "
                         f"{namelist_path}")

    extrusion = f"&extrusion\n" \
                f"domain_top={vertical_dimension.domain_top},\n" \
                f"method='{vertical_dimension.extrusion_method}',\n" \
                f"number_of_layers={vertical_dimension.number_of_layers},\n" \
                f"/\n\n"
    namelist = EXTRUSION_REGEX.sub(extrusion, namelist)

    with open(namelist_path, 'w') as file:
        file.write(namelist)
