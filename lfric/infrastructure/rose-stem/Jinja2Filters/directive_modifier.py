#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Implements a Jinja 2 filter which inserts values into batch manager directive
strings.
'''

import math
import re


def directive_modifier(directive, cores, walltime, xios_nodes=0, memory='30GB'):
    '''
    Substitutes values into batch manager directives.

    @param directive  Directive string or list of directives
    @param cores      Total system cores
    @param walltime   Walltime
    @param xios_nodes Number of xios nodes

    @return List of directives with substitutions made
    '''
    def choose_replacement(name, arguments, xios_nodes, memory):
        '''
        Converts a substitution name into a value.
        '''
        if name == 'nodes':
            node_size = int(arguments[0])
            middlebit = str(int(math.ceil(float(cores)
                                          / float(node_size))) + xios_nodes)
        elif name == 'cores':
            middlebit = str(cores)
        elif name == 'time_hhmmss':
            middlebit = str(walltime)
        elif name == 'memory':
            middlebit = str(memory)
        else:
            raise Exception('Unrecognised function name "{}"'.format(name))
        return middlebit

    if isinstance(directive, str):
        directive = [directive]

    pattern = re.compile(r'<<(\S+?)\((\S*?)\)>>')

    plist = []
    for directive_i in directive:
        processed = directive_i
        for match in pattern.finditer(directive_i):
            name = match.group(1)
            arguments = match.group(2).split(',')
            match_string = '<<%s(%s)>>' % (match.group(1), match.group(2))
            processed = processed.replace(match_string,
                                          choose_replacement(name,
                                                             arguments,
                                                             xios_nodes, memory))
        plist.append(processed)

    return plist


if __name__ == '__main__':
    DIRECTIVE = '-l=select=<<nodes(36)>>,walltime=<<time_hhmmss()>>'
    print(directive_modifier(DIRECTIVE, 122, '00:05:00'))
