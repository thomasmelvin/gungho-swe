#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Manages a database of dependency information.

import logging
import sqlite3
from abc import ABCMeta, abstractmethod
from pathlib import Path
from time import time


##############################################################################
# Databases throw this exception.
#
class DatabaseException(Exception):
    pass


##############################################################################
# Basic backend database functionality.
#
class _Database(metaclass=ABCMeta):
    ##########################################################################
    # Makes sure a described table exists in the database.
    #
    # Arguments:
    #   name    - String by which the table will be known.
    #   columns - List of name/type/modifiers tuples.
    #
    @abstractmethod
    def ensureTable(self, name, columns):
        pass

    ##########################################################################
    # Passes an SQL query to the backend.
    #
    # \deprecated Clearly this breaks all encapsulation. It should be replaced
    #             with a series of objects from which to build an abstract
    #             query.
    #
    # Arguments:
    #   query - Potentially multi-line SQL query.
    #
    @abstractmethod
    def query(self, query):
        pass


##############################################################################
# Database backend.
#
class SQLiteDatabase(_Database):
    ##########################################################################
    # Default constructor.
    #
    # Arguments:
    #   filename - The filename of the database.
    #
    def __init__(self, filename: Path):
        super(SQLiteDatabase, self).__init__()

        start_time = time()
        self._database = sqlite3.connect(str(filename), timeout=5.0)
        self._database.row_factory = sqlite3.Row
        message = 'Time to initialise database: {0}'
        logging.getLogger(__name__).debug(message.format(time() - start_time))

    ###########################################################################
    # Destructor.
    #
    def __del__(self):
        start_time = time()
        self._database.commit()
        self._database.close()
        message = 'Time to finalise database: {0}'
        logging.getLogger(__name__).debug(message.format(time() - start_time))

    ###########################################################################
    # Creates a table if it does not already exist.
    #
    # Arguments:
    #   name    - String by which table will be identified.
    #   columns - List of lists of column definition clauses.
    #             e.g. [['foo', 'integer', 'primary key'],
    #                   ['bar', 'integer']]
    #
    def ensureTable(self, name, columns):
        columnDefinitions = []
        for columnDetails in columns:
            columnDefinitions.append(' '.join(columnDetails))
        query = 'CREATE TABLE IF NOT EXISTS {} ( {} )'
        start_time = time()
        with self._database:
            columnList = ', '.join(columnDefinitions)
            self._database.execute(query.format(name, columnList))
        message = 'Time to ensure database table: {0} [{1}]'
        logging.getLogger(__name__).debug(message.format(time() - start_time,
                                                         name))

    ###########################################################################
    # Execute an SQL query against the database.
    #
    # Arguments:
    #   query - Either a string continaing a single SQL instruction or a list
    #           of strings, each containing a single SQL instruction.
    #
    def query(self, query):
        if isinstance(query, list):
            query = '; '.join(query) + ';'
        query = ' '.join(query.split())  # This wheeze collapses whitespace
        start_time = time()
        try:
            cursor = self._database.cursor()
            if query.count(';') > 1:
                cursor.executescript(query)
            else:
                cursor.execute(query)
            return cursor.fetchall()
        except sqlite3.IntegrityError as ex:
            raise DatabaseException(ex)
        finally:
            message = 'Time to query database: {0} [{1}]'
            logging.getLogger(__name__).debug(
                message.format(time() - start_time,
                               query))


###############################################################################
# Basic file dependencies.
#
class FileDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__(self, database):
        self._database = database
        self._database.ensureTable('file_dependency',
                                   (('file', 'TEXT', 'NOT NULL'),
                                    ('prerequisite', 'TEXT', 'NOT NULL'),
                                    ))

    ###########################################################################
    # Remove a file and all its dependencies from the database.
    #
    # Arguments:
    #   filename - The filename as it appears in the database.
    #
    def removeFile(self, filename):
        query = 'DELETE FROM file_dependency WHERE file="{}"'.format(filename)
        self._database.query(query)

    ###########################################################################
    # Remove all files anddependencies from the database.
    #
    def removeAllFileDependencies(self):
        query = 'DELETE FROM file_dependency'
        self._database.query(query)

    ###########################################################################
    # Add a dependency relationship to the database.
    #
    # Arguments:
    #   filename     - A filename string.
    #   prerequisite - The filename string of a file which the first depends
    #                  on.
    #
    def addFileDependency(self, filename, prerequisite):
        query = "INSERT INTO file_dependency VALUES ('{}','{}')"
        self._database.query(
            query.format(filename, prerequisite))

    ###########################################################################
    # Get all the file dependency relationships.
    #
    # Arguments:
    # Return:
    #   A generator yielding (filename, filename) tuples.
    #
    def getDependencies(self):
        query = 'SELECT * FROM file_dependency ORDER BY file'
        result = self._database.query(query)

        lastFile = None
        prerequisites = []
        for row in result:
            if row['file'] != lastFile:
                if prerequisites:
                    yield Path(lastFile), prerequisites
                lastFile = row['file']
                prerequisites = []

            prerequisites.append(Path(row['prerequisite']))

        if prerequisites:
            yield Path(lastFile), prerequisites


###############################################################################
# Fortran dependencies.
#
class FortranDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__(self, database):
        self._database = database

        self._database.ensureTable('fortran_unit_type',
                                   [('type', 'TEXT', 'PRIMARY KEY')])
        self._database.query(
            'INSERT OR IGNORE INTO fortran_unit_type VALUES("program")')
        self._database.query(
            'INSERT OR IGNORE INTO fortran_unit_type VALUES("module")')
        self._database.query(
            'INSERT OR IGNORE INTO fortran_unit_type VALUES("procedure")')

        self._database.ensureTable('fortran_dependency_type',
                                   [('type', 'TEXT', 'PRIMARY KEY')])
        self._database.query(
            'INSERT OR IGNORE INTO fortran_dependency_type VALUES("compile")')
        self._database.query(
            'INSERT OR IGNORE INTO fortran_dependency_type VALUES("link")')

        self._database.ensureTable('fortran_program_unit',
                                   (('unit', 'TEXT', 'PRIMARY KEY'),
                                    ('file', 'TEXT', 'NOT NULL'),
                                    ('type',
                                     'REFERENCES fortran_unit_type(type)')))
        self._database.ensureTable('fortran_unit_dependency',
                                   (('unit', 'TEXT', 'NOT NULL'),
                                    ('prerequisite', 'TEXT', 'NOT NULL'),
                                    ('type',
                                     'REFERENCES fortran_dependency_type('
                                     'type)')))

    ###########################################################################
    # Remove a Fortran source file and all its dependencies.
    #
    # Arguments:
    #   filename - Filename string as it appears in the database.
    #
    def removeFile(self, filename):
        query = [
            '''
            CREATE TEMPORARY TABLE _units AS SELECT unit
            FROM fortran_program_unit
            WHERE file="{filename}"
            '''.format(filename=filename),
            'DELETE FROM fortran_unit_dependency WHERE unit IN _units',
            'DROP TABLE _units',
            '''
            DELETE FROM fortran_program_unit WHERE file="{filename}"
            '''.format(filename=filename)
        ]
        self._database.query(query)

    ###########################################################################
    # Add a program to the database.
    #
    # Arguments:
    #   name     - Name of the program's program unit.
    #   filename - Source file in which program was found.
    #
    def addProgram(self, name, filename):
        # Changes are transacted to ensure other processes can't find the
        # database with half a program.
        query = "INSERT INTO fortran_program_unit VALUES ( '{}', '{}', " \
                "'program' )"
        self._database.query(query.format(name, filename))

    ###########################################################################
    # Add a module to the database.
    #
    # Arguments:
    #   name     - The name of the module's program unit.
    #   filename - The source file in which the modules was found.
    #
    def addModule(self, name, filename):
        try:
            query = "INSERT INTO fortran_program_unit VALUES ( '{}', '{}', " \
                    "'module' )"
            self._database.query(query.format(name, filename))
        except DatabaseException as ex:
            message = 'Unable to add module "{name}" from "{filename}": {ex}'
            newException = DatabaseException(message.format(name=name,
                                                            filename=filename,
                                                            ex=ex))
            newException.module = name
            newException.filename = filename
            raise newException

    ##########################################################################
    # Add a program-unit procedure to the database.
    #
    # Arguments:
    #   name     - The name of procedure's program unit.
    #   filename - The source file in which the procedure was found.
    #
    def addProcedure(self, name, filename):
        try:
            query = "INSERT INTO fortran_program_unit VALUES ( '{name}', " \
                    "'{filename}', 'procedure' )"
            self._database.query(query.format(name=name, filename=filename))
        except DatabaseException as ex:
            message = 'Unable to add procedure "{name}" from "{filename}": {' \
                      'ex}'
            newException = DatabaseException(message.format(name=name,
                                                            filename=filename,
                                                            ex=ex))
            newException.module = name
            newException.filename = filename
            raise newException

    ##########################################################################
    # Returns a list of all program units in the database.
    #
    def get_program_units(self):
        query = 'SELECT * FROM fortran_program_unit'
        rows = self._database.query(query)
        return [(row['unit'], Path(row['file'])) for row in rows]

    ###########################################################################
    # Add a compile dependency to the database.
    #
    # Arguments:
    #   unit
    #   prerequisite
    #
    def addCompileDependency(self, unit, prerequisite):
        query = "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', " \
                "'compile' )"
        self._database.query(query.format(unit, prerequisite))

    ###########################################################################
    # Add a link dependency to the database.
    #
    # Arguments:
    #   unit
    #   prerequisite
    #
    def addLinkDependency(self, unit, prerequisite):
        query = "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', " \
                "'link' )"
        self._database.query(query.format(unit, prerequisite))

    ###########################################################################
    # Gets all the programs from the database
    #
    # Arguments:
    # Return:
    #   List of program names.
    #
    def getPrograms(self):
        query = "SELECT unit FROM fortran_program_unit " \
                + " WHERE type='program' ORDER BY unit DESC"
        rows = self._database.query(query)
        return [row['unit'] for row in rows]

    ###########################################################################
    # Gets all the modules from the database
    #
    # Arguments:
    # Return:
    #   List of tuples (module name, file)
    #
    def getModules(self):
        query = "SELECT unit, file FROM fortran_program_unit " \
                + " WHERE type='module'"
        rows = self._database.query(query)
        return [(row['unit'], row['file']) for row in rows]

    ###########################################################################
    # Gets program unit dependencies for a program unit.
    #
    # This could be done (and was for a while) as a single SELECT statement.
    # The problem with that is it becomes hard to know when nothing is
    # returned because there is nothing to return and when it happens because
    # a dependency is missing from the database.
    #
    # Arguments:
    #   programUnit - Discover prerequisites for this.
    #
    # Return:
    #   Yields (unit, unit filename, prerequisite, prerequisite filename)
    #   tuples
    #
    def getLinkDependencies(self, programUnit):
        query = 'SELECT * FROM fortran_program_unit WHERE unit=\'{}\''
        unit_rows = self._database.query(query.format(programUnit))
        programUnitFilename = unit_rows[0]['file']

        query = '''
              SELECT unit, prerequisite FROM fortran_unit_dependency
              WHERE unit=\'{}\' and type=\'link\'
              ORDER BY unit
              '''.format(programUnit)
        rows = self._database.query(query)
        for unit, prerequisite in rows:
            query = 'SELECT file FROM fortran_program_unit WHERE unit=\'{}\''
            unit_rows = self._database.query(query.format(prerequisite))
            if len(unit_rows) == 0:
                message = 'Program unit "{}" requires "{}" but it was not ' \
                          'found in the database'
                raise DatabaseException(message.format(unit, prerequisite))
            yield programUnit, Path(programUnitFilename), prerequisite, Path(
                unit_rows[0]['file'])

    ###########################################################################
    # Gets all program unit dependencies from database for compilation.
    #
    # Arguments:
    #
    # Return:
    #   Yields (unit, unit filename, prerequisite, prerequisite filename)
    #   tuples
    #
    def getCompileDependencies(self):
        query = '''
                SELECT d.unit,
                       df.file AS unit_file,
                       df.type AS unit_type,
                       d.prerequisite,
                       pf.file AS prerequisite_file,
                       pf.type AS prerequisite_type
                FROM fortran_unit_dependency AS d,
                     fortran_program_unit AS df,
                     fortran_program_unit as pf
                WHERE df.unit=d.unit
                      AND pf.unit=d.prerequisite
                      AND d.type='compile'
                ORDER BY d.unit
                '''.strip()
        rows = self._database.query(query)

        for unit, unit_file, unit_type, prerequisite, prerequisite_file, \
            prerequisite_type \
                in rows:
            yield unit, Path(unit_file), unit_type, prerequisite, Path(
                prerequisite_file), prerequisite_type
