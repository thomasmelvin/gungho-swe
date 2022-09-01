##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(ROOT_DIR)/components/lfric-xios/source
export IGNORE_DEPENDENCIES += xios mod_wait
export EXTERNAL_STATIC_LIBRARIES += xios
export PRE_PROCESS_MACROS += USE_XIOS

.PHONY: import-lfric-xios
import-lfric-xios:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
            SOURCE_DIR=$(PROJECT_SOURCE) \
            OPTIMISATION_PATH=$(OPTIMISATION_PATH)
