#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Return a list of meshnames to be generated for the list of science resolutions
'''

def setMeshNames(resolutionDict, groups):
    '''
    @param [in] nested dictionarys of rose/science groups and their resolution dependencies
    @param [in] list of rose groups to be used
    @return list of all the different resolutions
    '''
    meshnames=set()
    for group in groups:
        for sgroup,groupList in resolutionDict[group].items():
            meshnames |= set(groupList)
        
    return list(meshnames)
