#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to convert a dictionary into assignment strings.
'''


def dict_to_assign(context, indent):
    '''
    Takes a dictionary and returns a string of assigments "key = value".

    @param [in] context Dictionary.
    @param [in] indent Indentation string
    @return String resulting from setting the environment.
    '''
    env_variables = []
    for key, value in context.items():
        if not isinstance(value, dict):
            env_variables.append('%s = %s' % (key, value))

    joining_str = '\n' + indent
    return joining_str.join(env_variables)
