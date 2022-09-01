#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to strip the resolutions info.
'''
from __future__ import absolute_import
from jinja2 import contextfilter


def unique_entries(inDict):
    ''' Returns unique entries from a dictionary of lists '''

    unique=set()
    for k,v in inDict.items():
        unique |= set(v)

    return list(unique)
