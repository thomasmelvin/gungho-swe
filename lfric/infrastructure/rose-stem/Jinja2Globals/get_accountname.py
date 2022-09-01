##############################################################################
# (c) Crown copyright 2019 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
from __future__ import absolute_import, print_function
from os import environ
from jinja2 import contextfunction

@contextfunction
def get_accountname(context):
    if 'BATCH_ACCOUNT_NAME' in context.vars:
        return context.vars['BATCH_ACCOUNT_NAME']
    else:
        return environ['USER']
