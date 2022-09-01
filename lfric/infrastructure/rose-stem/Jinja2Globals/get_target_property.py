##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
from __future__ import absolute_import, print_function

from jinja2 import contextfunction

@contextfunction
def get_target_property(context, target_dets, prop_path=None, compiler=None):
    if prop_path is not None:
        path_parts = prop_path.split('.')
    else:
        path_parts = []
    subtree = context.vars['target'][target_dets['platform']]
    while len(path_parts) > 0:
        part = path_parts.pop(0)
        if part == 'compiler':
            if compiler is None:
                path_parts.insert(0, target_dets['compiler'])
            else:
                path_parts.insert(0, compiler)
        subtree = subtree[part]
    return subtree
