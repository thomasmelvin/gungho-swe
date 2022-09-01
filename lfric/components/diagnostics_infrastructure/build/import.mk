##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(ROOT_DIR)/components/diagnostics_infrastructure/source
export IGNORE_DEPENDENCIES += fox_sax fox_common
export EXTERNAL_STATIC_LIBRARIES += FoX_sax FoX_common FoX_utils FoX_fsys

.PHONY: import-diagnostics_infrastructure
import-diagnostics_infrastructure:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
