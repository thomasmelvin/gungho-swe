#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to break appart a macro call and extract the
environment arguments.
'''
from __future__ import absolute_import
import ast
import re
from jinja2 import contextfilter
import utilities


@contextfilter
def get_env_macro(context, call):
    '''
    Takes a string and parses any instances of an env dictionary.

    @param [in,out] context Jinja2 instance to run macro against. Dummy var.
    @param [in]     call    Invocation string.
    @return Tuple of application name, Configuration name,
            Dictionary of environment arguments and Macro name.
    '''
    if call.find('(') == -1:
        macro_name = call
        arguments = []
    else:
        macro_name = call[:call.index('(')]
        macro_arguments = call[call.index('(')+1:call.rindex(')')]
        arguments = utilities.dictionary_from_arguments(macro_arguments)

    normal_arguments = [argument for argument in arguments
                        if argument.find('=') == -1]
    keyword_arguments = [argument for argument in arguments
                         if argument.find('=') != -1]

    argument_list = []
    for argument in normal_arguments:
        argument_list.append(ast.literal_eval(argument.strip()))

    argument_dictionary = {}
    for argument in keyword_arguments:
        key, value = re.split(' *= *', argument)
        argument_dictionary[key.strip()] = ast.literal_eval(value.strip())

    # We only do work on the 'env' dictionary
    if 'env' in argument_dictionary:
        environment_dictionary = argument_dictionary['env']
    else:
        environment_dictionary = {}

    if len(normal_arguments) >= 2:
        app_name = argument_list[0]
        key = app_name + '_' + argument_list[1]
        return_value = app_name, key, environment_dictionary, macro_name
    else:
        return_value = None, None, None, None

    return return_value
