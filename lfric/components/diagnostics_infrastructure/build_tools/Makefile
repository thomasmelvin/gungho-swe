##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

export ROOT_DIR   ?= ../../..

.PHONY: test
test: export PYTHONPATH := ${PYTHONPATH}:./diagnostics_metadata_reconfigurator/src
test:
	python -m pytest
