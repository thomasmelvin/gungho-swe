##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import filecmp
import pytest
from pathlib import Path

from entities import VerticalDimension
from extrusion_namelist_updater import update_extrusion_namelist

TEST_DIR = Path(__file__).parent
NAMELIST_PATH = TEST_DIR / Path('input/configuration.nml')
BAD_NAMELIST_PATH = TEST_DIR / Path('input/bad_configuration.nml')
NO_NAMELIST_PATH = TEST_DIR / Path('configuration.nml')
KGO_NAMELIST_PATH_1 = TEST_DIR / Path('kgos/configuration_1.nml')
KGO_NAMELIST_PATH_2 = TEST_DIR / Path('kgos/configuration_2.nml')

VERT_DIM_1 = VerticalDimension('model_vert_axis_1_full_levels',
                               name='axis_1_full_levels',
                               domain_top=42,
                               extrusion_method='uniform',
                               number_of_layers=10)
VERT_DIM_2 = VerticalDimension('model_vert_axis_2_full_levels',
                               name='axis_2_full_levels',
                               domain_top=1000,
                               extrusion_method='method',
                               number_of_layers=3)


def test_update_extrusion_namelist():
    """Test extrusion namelist correctly updated"""
    update_extrusion_namelist(NAMELIST_PATH, VERT_DIM_2)
    assert filecmp.cmp(NAMELIST_PATH, KGO_NAMELIST_PATH_2)

    update_extrusion_namelist(NAMELIST_PATH, VERT_DIM_1)
    assert filecmp.cmp(NAMELIST_PATH, KGO_NAMELIST_PATH_1)


def test_update_namelist_no_extrusion():
    """Test an error is raised if the file does not have a valid extrusion
    namelist"""
    with pytest.raises(ValueError) as excinfo:
        update_extrusion_namelist(BAD_NAMELIST_PATH, VERT_DIM_2)
        assert "Can't find valid extrusion namelist" in excinfo.value


def test_update_no_namelist():
    """Test an error is raised if the file does not exist"""
    with pytest.raises(FileNotFoundError):
        update_extrusion_namelist(NO_NAMELIST_PATH, VERT_DIM_2)
