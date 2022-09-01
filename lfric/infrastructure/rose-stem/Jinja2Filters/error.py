#!/usr/bin/env python
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
A Jinja2 filter which raises an exception, thus halting execution of the
template.
'''


def error(message):
    '''
    Generate an exception using the filter text as a message.

    @param [in] message Error message
    '''
    raise Exception(message)
