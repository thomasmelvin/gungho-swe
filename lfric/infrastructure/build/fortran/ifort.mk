##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
IFORT_VERSION := $(shell ifort -v 2>&1 \
                       | cut -d' ' -f3 \
                       | awk -F "." '/[0-9]\.[0-9]/ { yy = $$1 % 100; printf "%03i%02i%02i\n", yy,$$2,$$3}' )

$(info ** Chosen Intel Fortran compiler version $(IFORT_VERSION))

ifeq ($(shell test $(IFORT_VERSION) -lt 0150001; echo $$?), 0)
  $(error IFort is too old to build dynamo. Must be at least 15.0.1)
endif

F_MOD_DESTINATION_ARG = -module$(SPACE)
OPENMP_ARG            = -qopenmp
FORTRAN_RUNTIME       = 

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2 -fp-model strict
FFLAGS_RISKY_OPTIMISATION = -O3 -xhost
FFLAGS_DEBUG              = -g -traceback
FFLAGS_WARNINGS           = -warn all -warn errors
FFLAGS_UNIT_WARNINGS      = -warn all
FFLAGS_INIT               = -ftrapuv
FFLAGS_RUNTIME            = -check all -fpe0
# Option for checking code meets Fortran standard - currently 2008
FFLAGS_FORTRAN_STANDARD   = -stand f08

#########################################################################
# Application and file-specific options referenced in build/compile_options.mk files 
# 
# These variables need explanatory comments and need to be exported
#
# -qoverride-limits applied to PSy-layer code due to Intel compiler bug ref #1486
# When the Intel bug is fixed, this option will be removed by #1490
export FFLAGS_INTEL_FIX_ARG         = -qoverride-limits
########################################################################

# The "-assume realloc-lhs" switch causes Intel Fortran prior to v17 to
# actually implement the Fortran2003 standard. At version 17 it becomes the
# default behaviour.
ifeq ($(shell test "$(IFORT_VERSION)" -lt 0170000; echo $$?), 0)
  $(info ** Activating Intel "Make it work" switch for version earlier than 17)
  FFLAGS_COMPILER += -assume realloc-lhs
endif

LDFLAGS_COMPILER =

FPPFLAGS = -P
