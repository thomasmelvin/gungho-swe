Gets configured fields from a Rose suite and extracts their immutable metadata 
from a JSON file to be output in an XIOS iodef.xml file.
    - Sample Rose suite and immutable metadata file held in inputs/

To run:
    cd src/
    python metadata_reconfigurator.py <rose-suite-path> <immutable-data-file> 
                                <output-path> --iodef-template <template-path>

To run tests:
    cd test/
    pytest integration_test.py
