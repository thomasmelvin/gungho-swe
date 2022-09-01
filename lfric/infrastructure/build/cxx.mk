##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# We have to concern ourselves with C++ compilers as well, due to external
# dependencies
#
# These environment variables are set up:
#
# CXX: C++ compiler command
# CXX_COMPILER: The compiler determined to be in use
# CXX_RUNTIME_LIBRARY: Name of the C++ runtime library used by this compiler
#
###############################################################################
# We assume GCC is in use unless the CXX environment variable is set.
#
CXX ?= g++

# Try to work out the compiler from CXX. Potentially this is a full path with
# version appended. i.e. [/absolute/path/]compiler[-version]
#
# Thus we take the leaf name only using "notdir" to discard the path. We then
# substitute hyphens for spaces and take only the first word to discard any
# version identifier.
#
CXX_COMPILER := $(firstword $(subst -, ,$(notdir $(CXX))))

ifdef CRAY_ENVIRONMENT
  ifeq '$(PE_ENV)' 'CRAY'
    CXX_COMPILER = craycc
  else ifeq '$(PE_ENV)' 'INTEL'
    CXX_COMPILER = icc
  else ifeq '$(PE_ENV)' 'GNU'
    CXX_COMPILER = g++
  else ifeq '$(PE_ENV)' 'PGI'
    CXX_COMPILER = pgc++
  else
    $(error Unrecognised Cray programming environment)
  endif
endif

include $(LFRIC_BUILD)/cxx/$(CXX_COMPILER).mk
export CXX_RUNTIME_LIBRARY
