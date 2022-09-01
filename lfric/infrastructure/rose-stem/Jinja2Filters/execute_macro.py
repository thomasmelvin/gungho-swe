#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja2 filter to run a macro specified by a string.
'''
from __future__ import absolute_import
import ast
import re
from jinja2 import contextfilter
import utilities


@contextfilter
def execute_macro(context, call):
    '''
    Takes a string and executes it as though it were a Jinja2 macro call.

    The call string has the syntax <macro name>([<argument>]...).

    Arguments can be either position or keyword based.

    @param [in,out] context Jinja2 instance to run macro against.
    @param [in]     call     Invocation string.
    @return String resulting from calling the macro.
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

    return context.vars[macro_name](*argument_list, **argument_dictionary)
