##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Locate and compile all Fortran source in the current directory.
#
# This is done as a standalone make file so that your base makefile does not
# require the dependency analysis. This is annoying if you issue a
# "make clean" on an already clean working copy. In that case the analysis
# will be performed in order that it then be deleted.
#
# The following variables may be specified to modify the compile process...
#
# ROOT: Project directory
# BIN_DIR: Path to directory for resulting executables.
#          Default: $(ROOT)/bin
# CXX_LINK: Set this macro to have the C++ runtime library linked to the
#           executable.
# EXTERNAL_LIBRARIES: Libraries required by link stage. These should be
#                     specified in static link order, even if you are linking
#                     dynamically.
# FFLAG_GROUPS: Space separated list of FFLAG_<group name> variables to use in
#               building up the FFLAGS variable. Only the group name is
#               specified.
# LINK_TYPE: 'static' or 'dynamic' linking.
#            Default: dynamic, except on Crays where it's static
# PROGRAMS: Names of programs to compile.
#           Default: Everything listed in programs.mk
# PRE_PROCESS_MACROS: Macro definitions in the form NAME[=MACRO] to be passed
#                     to the compiler.
# PROJECT_MAKE_DIR: Used to locate project specific targets modifiers and
#                   such.
# COMPILE_OPTIONS: Name of an optional file that can be included to list
#                  project specific compile options
#
##############################################################################

.SECONDEXPANSION:

# Build a set of "-I" arguments to seach the whole object tree:
INCLUDE_ARGS := $(subst ./,-I,$(shell find . -mindepth 1 -type d -print))

# Build a set of "-D" argument for any pre-processor macros
#
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

include programs.mk

PROGRAMS ?= $(basename $(notdir $(PROG_OBJS)))

# Convert the program names to all caps, append "_OBJS" to the end and
# dereference that to get a list of all objects needed by all programs:
#
ALL_OBJECTS = $(foreach proj, $(shell echo $(PROGRAMS) | tr a-z A-Z), $($(proj)_OBJS))

-include $(COMPILE_OPTIONS)

.PHONY: applications
applications: FFLAGS += $(foreach group, $(FFLAG_GROUPS), $(FFLAGS_$(group)))
applications: $(addprefix $(BIN_DIR)/, $(PROGRAMS))

##############################################################################

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk
include $(LFRIC_BUILD)/cxx.mk
-include $(COMPILE_OPTIONS)

BIN_DIR ?= $(ROOT)/bin

# If the compiler produces module files, tell it where to put them
#
ifdef F_MOD_DESTINATION_ARG
  MODULE_DESTINATION_ARGUMENT = $(F_MOD_DESTINATION_ARG)$(dir $@)
endif

ifdef CXX_LINK
  EXTERNAL_DYNAMIC_LIBRARIES += $(CXX_RUNTIME_LIBRARY)
endif

ifdef CRAY_ENVIRONMENT
  $(warning Running on a Cray, selecting static linking)
  LINK_TYPE ?= static
else
  LINK_TYPE ?= dynamic
endif

# Work out what to do with external libraries.
#
ifeq "$(LINK_TYPE)" "static"
  override EXTERNAL_STATIC_LIBRARIES  := $(EXTERNAL_STATIC_LIBRARIES) $(EXTERNAL_DYNAMIC_LIBRARIES)
  override EXTERNAL_DYNAMIC_LIBRARIES :=
else ifeq "$(LINK_TYPE)" "dynamic"
  # Nothing further needs to be done.
else
  $(error Unrecognised LINK_TYPE. Must be either "static" or "dynamic")
endif

##############################################################################
$(BIN_DIR)/%: %.x | $(BIN_DIR)
	$(call MESSAGE,Installing,$*)
	$(Q)cp $< $(BIN_DIR)/$(notdir $*)

.PRECIOUS: %.x
%.x: $$($$(shell basename $$* | tr a-z A-Z)_OBJS)
	$(call MESSAGE,Linking,$@)
	$(Q)$(LDMPI) $(LDFLAGS) -o $@ \
	             $^ \
	             $(patsubst %,-l%,$(EXTERNAL_STATIC_LIBRARIES)) \
	             $(patsubst %,-l%,$(EXTERNAL_DYNAMIC_LIBRARIES))

.PRECIOUS: %.o
%.o: %.f90
	$(call MESSAGE,Compile,$<)
	$(Q)$(FC) $(FFLAGS) \
	          $(MODULE_DESTINATION_ARGUMENT) \
	          $(INCLUDE_ARGS) -c -o $(basename $@).o $<

%.o: %.F90
	$(call MESSAGE,Pre-process and compile,$<)
	$(Q)$(FC) $(FFLAGS) \
	          $(MODULE_DESTINATION_ARGUMENT) \
	          $(INCLUDE_ARGS) $(MACRO_ARGS) -c -o $(basename $@).o $<

#############################################################################
# Directories
#
$(BIN_DIR):
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@

#############################################################################
# Dependencies
#
include dependencies.mk
