#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to strip the resolutions info.
'''
from __future__ import absolute_import
import ast
import re
import utilities


def get_precision_macro(call):
    '''
    Identifies/Processes the precision argument from the string input.

    The precision argument lists the requested default precisions for
    real variables. If present, the argment is processed to extract
    requested real-variable precision levels for the given
    application/configuration.
    ...
    @param [in] call Invocation string.
    @return List resulting from retrieving the requested precisions.
    '''
    # pylint: disable=too-many-locals, too-many-branches

    if call.find('(') == -1:
        arguments = []
    else:
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

    if len(argument_list) > 0:
        app_name = argument_list[0]
        app_key = app_name

        if len(argument_list) > 1:
            configuration = argument_list[1]
            app_key += '_' + configuration

        argument_dictionary = {}
        for argument in keyword_arguments:
            key, value = re.split(' *= *', argument)
            argument_dictionary[key.strip()] = ast.literal_eval(value.strip())

        # Return info about the precision argument
        resource_precisions = []
        if 'precisions' in argument_dictionary.keys():
            resource_precisions = argument_dictionary['precisions']

        return_value = (app_key,
                        resource_precisions)

    else:

        return_value = None, None

    return return_value
