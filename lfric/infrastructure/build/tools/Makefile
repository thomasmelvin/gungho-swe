##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

HERE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

export PYTHONPATH := ${PYTHONPATH}:$(HERE_DIR):$(HERE_DIR)/diagnostics_metadata_reconfigurator/src

ifdef VERBOSE
    VERBOSE_ARG = -v
endif
ifdef DEBUG
    DEBUG_ARG = -vv
endif

.PHONY: test
test: pytest

.PHONY: pytest
pytest: mypy
	pytest $(VERBOSE_ARG) $(DEBUG_ARG)

STYLED_FILES = configurator
STYLED_APPS = GenerateFeigns GenerateLoader GenerateNamelist

# We pass the applications to mypy separately as they fight with each other to
# declare a __main__. This means the package modules get extra tested.
#
.PHONY: mypy
mypy: $(addprefix mypy/,$(STYLED_APPS))

.PHONY: mypy/%
mypy/%:
	mypy $(VERBOSE_ARG) $(STYLED_FILES) $*

.PHONY: style
style: pylint

# We pass the applications to pylint separately as they fight with each other
# to declare a __main__. This means the package modules get extra tested.
#
.PHONY: pylint
pylint: pycodestyle $(addprefix pylint/,$(STYLED_APPS))
	pylint $(VERBOSE_ARG) $(STYLED_FILES)


.PHONY: pylint/%
pylint/%:
	pylint $(VERBOSE_ARG) $*

.PHONY: pycodestyle
pycodestyle:
	pycodestyle $(VERBOSE_ARG) $(STYLED_APPS) $(STYLED_FILES)

