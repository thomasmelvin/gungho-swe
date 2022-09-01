##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Extract metadata from a Rose suite and an immutable metadata JSON file and
    output it in an XIOS iodef.xml file"""

import argparse
import logging
from pathlib import Path

from metadata_extractor import MetadataExtractor
from metadata_iodef_generator import MetadataIodefGenerator
from metadata_validator import validate_metadata

LOGGER = logging.getLogger('reconfigurator')


def get_args():
    """Get command line arguments

    :return: List of values given from arguments
    """
    parser = argparse.ArgumentParser(
        description="This script parses diagnostics metadata from a rose "
                    "suite and a JSON file containing the immutable metadata. "
                    "It then adds it to either a new or existing XIOS "
                    "iodef.xml file")
    parser.add_argument('diagnostic_config_path', type=str,
                        help="Path to a rose-app.conf file. This should be "
                             "the file describing the diagnostic "
                             "configuration")
    parser.add_argument('immutable_data_path', type=str,
                        help="Path to the JSON file containing immutable "
                             "metadata")
    parser.add_argument('output_file_path', type=str,
                        help="Path to output the iodef.xml file to, including "
                             "filename and extension")
    parser.add_argument('--iodef-template', dest='template_path',
                        help="Path to a template iodef.xml file. This is an "
                             "existing iodef file to which the diagnostics "
                             "information will be added. This file should "
                             "contain a context node containing axis, grid, "
                             "field and file definition nodes")
    parser.add_argument('--namelist', dest='namelist_path',
                        help="Path to an LFRic configuration.nml file "
                             "containing a valid extrusion namelist. This "
                             "namelist will be updated with the values parsed "
                             "from the diagnostic configuration")
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        dest='verbosity',
                        help="Increases logging verbosity. -v logs INFO and "
                             "above, -vv logs DEBUG and above. "
                             "Default is WARNING and above")
    parser.add_argument('--force', '-f', dest='force', action='store_true',
                        help="Forces the script to complete and produce an "
                             "iodef.xml output file, even if the metadata is "
                             "not valid")
    parser.add_argument('--log-stdout', dest='log_stdout', action='store_true',
                        help="Writes log messages to file as well as the "
                             "console. This file, log.txt, will be saved in "
                             "the same directory as the iodef.xml file")
    return parser.parse_args()


def setup_logging(verbosity: int, iodef_file_path: str, log_stdout: bool)\
            -> None:
    """Configure logging level and log file if one is being used"""
    logging_levels = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    verbosity = logging_levels[min(verbosity, 2)]
    handlers = [logging.StreamHandler()]

    if log_stdout:
        # Put log in same directory as output XML file
        log_path = Path(iodef_file_path).expanduser().parent / "log.txt"
        handlers.append(logging.FileHandler(log_path, mode='w'))

    logging.basicConfig(format="%(levelname)s:%(name)s: %(message)s",
                        level=verbosity,
                        handlers=handlers)


def reconfigure_metadata():
    """Extract metadata from a Rose suite and an immutable metadata JSON file
        and output it in an XIOS iodef.xml file"""
    args = get_args()
    setup_logging(args.verbosity, args.output_file_path, args.log_stdout)
    meta_extractor = MetadataExtractor(args.diagnostic_config_path,
                                       args.immutable_data_path)
    metadata = meta_extractor.extract_metadata()
    validate_metadata(metadata, args.force)
    generator = MetadataIodefGenerator(metadata)
    generator.generate(args.output_file_path, args.template_path,
                       args.namelist_path)


if __name__ == '__main__':
    reconfigure_metadata()
