##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import filecmp
import logging
from xml.etree import ElementTree
from pathlib import Path

from metadata_extractor import MetadataExtractor
from metadata_iodef_generator import MetadataIodefGenerator
from metadata_validator import validate_metadata

TEST_DIR = Path(__file__).parent
IMMUTABLE_DATA_PATH = TEST_DIR / 'input/LFRic_meta_data.JSON'
ROSE_APP_PATH = TEST_DIR / 'input/rose-suite/rose-app.conf'
TEST_IODEF_PATH1 = TEST_DIR / 'kgos/test_iodef.xml'
OUTPUT_PATH1 = TEST_DIR / 'output/iodef.xml'

TEMPLATE_PATH = TEST_DIR / 'templates/minimal_iodef.xml'
TEST_IODEF_PATH2 = TEST_DIR / 'kgos/test_iodef2.xml'
OUTPUT_PATH2 = TEST_DIR / 'output/iodef_from_template.xml'

LOGGER = logging.getLogger()


class CommentedTreeBuilder(ElementTree.TreeBuilder):
    """XML tree builder to enable parser to handle comments"""
    def comment(self, data):
        self.start(ElementTree.Comment, {})
        self.data(data)
        self.end(ElementTree.Comment)


def test_without_template():

    extractor = MetadataExtractor(ROSE_APP_PATH, IMMUTABLE_DATA_PATH)
    metadata = extractor.extract_metadata()
    validate_metadata(metadata)
    generator = MetadataIodefGenerator(metadata)
    generator.generate(OUTPUT_PATH1)

    assert filecmp.cmp(TEST_IODEF_PATH1, OUTPUT_PATH1)


def test_with_template():
    extractor = MetadataExtractor(ROSE_APP_PATH, IMMUTABLE_DATA_PATH)
    metadata = extractor.extract_metadata()
    validate_metadata(metadata)
    generator = MetadataIodefGenerator(metadata)
    generator.generate(OUTPUT_PATH2, TEMPLATE_PATH)

    assert filecmp.cmp(TEST_IODEF_PATH2, OUTPUT_PATH2)
