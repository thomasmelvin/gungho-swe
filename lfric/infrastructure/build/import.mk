##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export IGNORE_DEPENDENCIES += netcdf MPI yaxt pfunit_mod
export EXTERNAL_DYNAMIC_LIBRARIES += yaxt yaxt_c netcdff netcdf hdf5

.PHONY: import-infrastructure
import-infrastructure:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
