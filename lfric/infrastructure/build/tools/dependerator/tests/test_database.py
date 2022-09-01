#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from pathlib import Path

import pytest

from dependerator.database import (DatabaseException,
                                   FileDependencies,
                                   FortranDependencies,
                                   SQLiteDatabase)


class TestDatabase():
    def test_all(self, tmp_path: Path):
        """
        Exercise basic database functions of the SQLite wrapper.
        """
        uut = SQLiteDatabase(str(tmp_path / 'test.db'))

        uut.ensureTable('test', ["'alpha' integer primary key"])

        result = uut.query('select * from test')
        assert 0 == len(result)

        result = uut.query("insert into 'test' values (1)")
        assert 0 == len(result)

        result = uut.query('select * from test')
        assert 1 == len(result)
        assert 1 == result[0][0]


class TestFileDependency():
    def test_all(self, tmp_path: Path):
        """
        Ensure file dependencies are correctly managed.
        """
        database = SQLiteDatabase(str(tmp_path / 'file.db'))
        uut = FileDependencies(database)

        result = uut.getDependencies()
        assert [] == list(result)

        uut.addFileDependency('foo.f90', 'bar')
        result = uut.getDependencies()
        assert [(Path('foo.f90'), [Path('bar')])] == list(result)

        uut.addFileDependency(Path('foo.f90'), Path('baz'))
        uut.addFileDependency(Path('qux.f90'), Path('bar'))
        result = uut.getDependencies()
        assert [(Path('foo.f90'), [Path('bar'), Path('baz')]),
                (Path('qux.f90'), [Path('bar')])] == list(result)

        uut.removeFile('foo.f90')
        result = uut.getDependencies()
        assert [(Path('qux.f90'), [Path('bar')])] == list(result)

        uut.removeAllFileDependencies()
        result = uut.getDependencies()
        assert [] == list(result)


class TestFortranDependency():
    def test_add_program(self, tmp_path: Path):
        """
        Ensure that a program can be added and correctly retrieved.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        uut.addProgram('foo', 'bar.f90')

        result = uut.getPrograms()
        assert [u'foo'] == result

        # Should test that filename is correctly stored.

    def test_remove_source_file(self, tmp_path: Path):
        """
        Ensure that all traces of a file are removed correctly.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        self._populateDB(uut)

        # Delete something from the middle of a dependency chain
        uut.removeFile('bar.f90')

        result = uut.getPrograms()
        assert [u'fred', u'foo'] == list(result)

        result = uut.getCompileDependencies()
        assert [(u'foo', Path('foo.f90'), u'program', u'baz', Path('baz.f90'), u'module'),
                (u'fred', Path('fred.f90'), u'program', u'wilma', Path('wilma.f90'), u'module')] \
            == list(result)

        assert [(u'baz', Path('baz.f90'),
                 u'foo', Path('foo.f90'))] \
            == list(uut.getLinkDependencies('baz'))
        assert [(u'wilma',
                 Path('wilma.f90'),
                 u'fred',
                 Path('fred.f90'))] == list(uut.getLinkDependencies('wilma'))

    def test_programs(self, tmp_path: Path):
        """
        Ensure the the list of programs is correctly retrieved.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        self._populateDB(uut)

        programs = uut.getPrograms()
        assert [u'fred', u'foo'] == list(programs)

    def test_modules(self, tmp_path: Path):
        """
        Ensure the the list of modules is correctly retrieved.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        self._populateDB(uut)

        modules = uut.getModules()
        assert [(u'bar', u'bar.f90'), (u'baz', u'baz.f90'), (u'qux', u'qux.f90'), (u'wilma', u'wilma.f90')] == list(modules)

    def test_get_all_file_dependencies(self, tmp_path: Path):
        """
        Ensure all file dependencies are correctly returned.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        self._populateDB(uut)

        dependencies = list(uut.getCompileDependencies())
        assert [(u'bar', Path('bar.f90'), u'module', u'qux', Path('qux.f90'), u'module'),
                (u'foo', Path('foo.f90'), u'program', u'bar', Path('bar.f90'), u'module'),
                (u'foo', Path('foo.f90'), u'program', u'baz', Path('baz.f90'), u'module'),
                (u'fred', Path('fred.f90'), u'program', u'wilma', Path('wilma.f90'), u'module')] \
            == dependencies

        assert [(u'bar', Path('bar.f90'),
                 u'foo', Path('foo.f90'))] == list(uut.getLinkDependencies('bar'))
        assert [(u'baz', Path('baz.f90'),
                 u'foo', Path('foo.f90'))] == list(uut.getLinkDependencies('baz'))
        assert [(u'qux', Path('qux.f90'),
                 u'bar', Path('bar.f90'))] == list(uut.getLinkDependencies('qux'))
        assert [(u'wilma',
                 Path('wilma.f90'),
                 u'fred',
                 Path('fred.f90'))] == list(uut.getLinkDependencies('wilma'))

    def test_duplicate_module(self, tmp_path: Path):
        """
        Ensure that an exception is thrown if an attempt is made to add a
        module which is already in the database.
        """
        database = SQLiteDatabase(str(tmp_path / 'fortran.db'))
        uut = FortranDependencies(database)
        self._populateDB(uut)

        with pytest.raises(DatabaseException) as caught:
            uut.addModule('qux', 'cheese/qux.f90')
        assert 'qux' == caught.value.module
        assert 'cheese/qux.f90' == caught.value.filename

    def _populateDB(self, database):
        """
        Helper which creates a database full of test values.
        """
        database.addProgram('foo', 'foo.f90')
        database.addModule('bar', 'bar.f90')
        database.addModule('baz', 'baz.f90')
        database.addModule('qux', 'qux.f90')
        database.addProgram('fred', 'fred.f90')
        database.addModule('wilma', 'wilma.f90')

        database.addCompileDependency('foo', 'bar')
        database.addCompileDependency('foo', 'baz')
        database.addCompileDependency('bar', 'qux')
        database.addCompileDependency('fred', 'wilma')

        database.addLinkDependency('bar', 'foo')
        database.addLinkDependency('baz', 'foo')
        database.addLinkDependency('qux', 'bar')
        database.addLinkDependency('wilma', 'fred')
