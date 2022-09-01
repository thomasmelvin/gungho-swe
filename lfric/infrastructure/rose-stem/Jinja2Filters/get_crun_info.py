#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to strip the crun info.
'''
import re
from jinja2 import contextfilter


@contextfilter
def get_crun_info(call):
    '''
    Takes a string return dictionary values related to crunning:
        crun: Number of runs to do in the crun.
        ...
    @param [in] context Jinja2 instance to run macro against.
    @param [in] call    Invocation string.
    @return String resulting from setting the environment.
    '''
    if call.find('(') == -1:
        arguments = []
    else:
        arguments = re.split(', *', call[call.index('(')+1:call.rindex(')')])

    normal_arguments = [argument for argument in arguments
                        if argument.find('=') == -1]
    keyword_arguments = [argument for argument in arguments
                         if argument.find('=') != -1]

    argument_list = []
    for argument in normal_arguments:
        if argument[0] == '"':
            argument_list.append(argument[1:-1])
        else:
            argument_list.append(argument)

    argument_dictionary = {}
    for argument in keyword_arguments:
        key, value = re.split(' *= *', argument)
        argument_dictionary[key] = value

    # Return info about the crun arguments
    return_value = {}
    if 'crun' in argument_dictionary:
        return_value.update({'crun': int(argument_dictionary['crun'])})

    return return_value
