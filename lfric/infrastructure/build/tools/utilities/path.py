#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# A library of filename tools.

from __future__ import print_function

from __future__ import absolute_import
import os.path

###############################################################################
def replaceExtension( filename, extension ):
    (base, rubbish) = os.path.splitext( filename )
    return '{}.{}'.format( base, extension)