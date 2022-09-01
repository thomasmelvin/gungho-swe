#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Process previously analysed dependency database. For fun and profit!

import logging
import os.path

from utilities.path import replaceExtension


###############################################################################
# Process dependency database.
#
class FortranProcessor():
    ###########################################################################
    # Constructor.
    #
    # Arguments:
    #   database  - FortranDatabase object holding details.
    #   objectdir - The directory which holds .o files.
    #   moduledir - The directory which holds .mod files.
    #
    def __init__(self, database, objectDirectory, moduleDirectory):
        self._database = database
        self._objectDirectory = objectDirectory
        self._moduleDirectory = moduleDirectory

    ###########################################################################
    # Examine the program unit dependecies and work out the file dependencies.
    #
    # Arguments:
    #   fileStore - FileDependencies object to accept computed dependencies.
    #   moduleObjects - Boolean, whether the compiler stores module
    #   information in object files

    def determineCompileFileDependencies(self, fileStore, moduleObjects=False):

        logging.getLogger(__name__).info(
            'Removing old file compile dependencies')

        fileStore.removeAllFileDependencies()

        logging.getLogger(__name__).info(
            'Determining file compile dependencies...')

        # we have 2 types of dependency
        # 1) module file dependencies: a module's .mod file depends on it's
        # source files .o file (ommitted of course if module information is
        # stored in object files)
        # 2) object file dependencies: %.o files depend on .mod files for
        # modules used within %.fF90 (or on object files if module
        # information stored in object files)

        if not moduleObjects:
            for module, file in self._database.getModules():
                moduleFilename = module + '.mod'
                modulePathname = os.path.join(self._moduleDirectory,
                                              os.path.split(file)[0],
                                              moduleFilename)

                sourceObjectFilename = replaceExtension(file, 'o')
                sourceObjectPathname = os.path.join(self._objectDirectory,
                                                    sourceObjectFilename)

                fileStore.addFileDependency(modulePathname,
                                            sourceObjectPathname)

        for unit, unitFilename, unitType, prerequisite, \
            prerequisiteFilename, prerequisiteType in \
                self._database.getCompileDependencies():

            unitObjectFilename = replaceExtension(unitFilename, 'o')
            unitObjectPathname = os.path.join(self._objectDirectory,
                                              unitObjectFilename)

            prereqModuleFilename = prerequisite + '.mod'
            prereqModulePathname = os.path.join(self._moduleDirectory,
                                                os.path.split(
                                                    prerequisiteFilename)[0],
                                                prereqModuleFilename)

            prereqObjectFilename = replaceExtension(prerequisiteFilename,
                                                    'o')
            prereqObjectPathname = os.path.join(self._objectDirectory,
                                                prereqObjectFilename)

            if prerequisiteType == 'module':
                message = '{0} depends on module {1}'.format(unit,
                                                             prerequisite)
                logging.getLogger(__name__).info(message)

                if moduleObjects:
                    fileStore.addFileDependency(unitObjectPathname,
                                                prereqObjectPathname)
                else:
                    fileStore.addFileDependency(unitObjectPathname,
                                                prereqModulePathname)

            if prerequisiteType == 'procedure':
                message = '{0} depends on procedure {1}'.format(unit,
                                                                prerequisite)
                logging.getLogger(__name__).info(message)

                fileStore.addFileDependency(unitObjectPathname,
                                            prereqObjectPathname)

    ###########################################################################
    # Determine all program units needed to build each program.
    #
    # TODO: Once we have a more recent version of SQLite we could consider
    # doing this at the database level.
    #
    def determineLinkDependencies(self):
        for program in self._database.getPrograms():
            logging.getLogger(__name__).info('Program {0}'.format(program))

            prerequisites = set()
            self._descend(program, prerequisites)
            unit, unit_file, prereq, prereq_file = next(
                self._database.getLinkDependencies(program))
            program_object_file = os.path.join(self._objectDirectory,
                                               replaceExtension(unit_file,
                                                                'o'))
            prerequisites.add(program_object_file)
            yield os.path.join(self._objectDirectory, program), \
                  sorted(list(prerequisites))

    ##########################################################################
    def _descend(self, programUnit, prerequisites):
        logger = logging.getLogger(__name__)
        for unit, unit_file, prereq, prereq_file \
                in self._database.getLinkDependencies(programUnit):
            logger.info('  Requires {0}'.format(prereq))

            prereq_object_file = os.path.join(self._objectDirectory,
                                              replaceExtension(prereq_file,
                                                               'o'))

            if prereq_object_file in prerequisites:
                logger.info('    Seen already, stopping descent')
            else:
                prerequisites.add(prereq_object_file)
                self._descend(prereq, prerequisites)
