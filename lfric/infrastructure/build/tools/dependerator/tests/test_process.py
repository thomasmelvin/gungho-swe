#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from pathlib import Path
import pytest

from dependerator.database import (FileDependencies,
                                   FortranDependencies,
                                   SQLiteDatabase)
from dependerator.process import FortranProcessor


class TestFortranProcessor():
    @pytest.fixture
    def databases(self, tmp_path_factory):
        filename = tmp_path_factory.mktemp('db-', True) / 'fortran.db'
        database = SQLiteDatabase(filename)
        return FortranDependencies(database), FileDependencies(database)

    def test_compile_dependencies(self, databases):
        self._populate_database(databases[0])

        uut = FortranProcessor(databases[0], "objects", "modules")
        uut.determineCompileFileDependencies(databases[1])

        assert [(Path('modules/bits/bar.mod'), [Path('objects/bits/bar.o')]),
                 (Path('modules/bits/baz.mod'), [Path('objects/bits/baz.o')]),
                 (Path('modules/bobs/corge.mod'), [Path('objects/bobs/grault.o')]),
                 (Path('modules/bobs/qux.mod'), [Path('objects/bobs/qux.o')]),
                 (Path('objects/bits/bar.o'), [Path('modules/bits/baz.mod')]),
                 (Path('objects/bobs/grault.o'),[Path('modules/bits/bar.mod')]),
                 (Path('objects/bobs/qux.o'), [Path('modules/bits/baz.mod')]),
                 (Path('objects/foo.o'), [Path('modules/bits/bar.mod'), Path('objects/quux.o')])] \
        == list(databases[1].getDependencies())

    def test_compile_dependencies_module_objects(self, databases):
        self._populate_database(databases[0])

        uut = FortranProcessor(databases[0], "objects", "modules")
        uut.determineCompileFileDependencies(databases[1], moduleObjects=True)

        assert [(Path('objects/bits/bar.o'), [Path('objects/bits/baz.o')]),
                (Path('objects/bobs/grault.o'), [Path('objects/bits/bar.o')]),
                (Path('objects/bobs/qux.o'), [Path('objects/bits/baz.o')]),
                (Path('objects/foo.o'), [Path('objects/bits/bar.o'), Path('objects/quux.o')])] \
               == list(databases[1].getDependencies())

    def test_link_dependencies(self, databases):
        self._populate_database(databases[0])

        import sys
        uut = FortranProcessor(databases[0], "objects", "modules")
        result = list(uut.determineLinkDependencies())

        assert [(u'objects/foo', ['objects/bits/bar.o',
                                  'objects/bits/baz.o',
                                  'objects/bobs/qux.o',
                                  'objects/foo.o', 'objects/quux.o'])] == result

    def _populate_database( self, database: FortranDependencies ):
        database.addProgram( 'foo', 'foo.f90' )
        database.addModule( 'bar', 'bits/bar.f90' )
        database.addModule( 'baz', 'bits/baz.f90' )
        database.addModule( 'qux', 'bobs/qux.f90' )
        database.addModule( 'corge', 'bobs/grault.f90')
        database.addProcedure( 'quux', 'quux.f90')

        database.addCompileDependency( 'foo', 'bar' )
        database.addCompileDependency( 'bar', 'baz' )
        database.addCompileDependency( 'qux', 'baz' )
        database.addCompileDependency( 'corge', 'bar')
        database.addCompileDependency( 'foo', 'quux')

        database.addLinkDependency( 'foo', 'bar' )
        database.addLinkDependency( 'bar', 'baz' )
        database.addLinkDependency( 'baz', 'qux' )
        database.addLinkDependency('corge', 'bar')
        database.addLinkDependency('foo', 'quux')

if __name__ == '__main__':
    unittest.main()
