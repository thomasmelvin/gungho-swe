##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Automate developer checks.
#
RUN_DIR = $(dir $(lastword $(MAKEFILE_LIST)))

Q=@

# TODO: These tools were written a long time ago before we decided to adopt
# PEP-8 or really knew what we were doing. We are styling files as we touch
# them so the untouched need to be ignored. Eventually the ability to ignore
# files should be removed.
#
STYLEABLE = parserender/mirrorer.py parserender/tests
I=-

.PHONY: all
all: pytest pylint pycodestyle

.PHONY: pytest
pytest:
	$Qecho === Unit testing
	$Qpytest $(RUN_DIR)

.PHONY: pylint
pylint:
	$Qecho === Linting
	$Q$Ipylint $(STYLEABLE)

.PHONY: pycodestyle
pycodestyle:
	$Qecho === Styling
	$Qpycodestyle $(STYLEABLE)
