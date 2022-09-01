#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter which takes a number of strings and joins them with
the filter text.

e.g.

' foo ' | command_list('comm1', ['comm2, 'comm3'])
=> comm1 foo comm2 foo comm3
'''


def command_list(separator, *fragments):
    '''
    Takes an arbitrary collection of strings and joins them with the filter
    text into a single string.

    @param [in] separator Filter text string.
    @param [in] fragments Multiple strings or lists of strings.
    @return String resulting from joining the fragments
    '''
    commands = []
    for fragment in fragments:
        if hasattr(fragment, '__iter__'):
            commands.extend(fragment)
        else:
            commands.append(fragment)
    return separator.join(commands)
