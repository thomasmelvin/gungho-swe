#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
from pathlib import Path
import pytest
from textwrap import dedent

from dependerator.analyser import (FortranAnalyser,
                                   NamelistDescriptionAnalyser)
from dependerator.database import (FileDependencies,
                                   FortranDependencies,
                                   SQLiteDatabase)


class TestNamelistDescriptionAnalyser():
    @pytest.fixture
    def database(self, tmp_path_factory):
        filename = tmp_path_factory.mktemp('db-', True) / 'test.db'
        database = SQLiteDatabase(filename)
        return FileDependencies(database)

    def testAnalysis(self, database, tmp_path: Path):
        """
        Presents a namelist description to the analyser and ensures the correct
        dependencies are returned.
        """
        test_filename = tmp_path / 'test.nld'
        test_filename.write_text(dedent('''
            namelist foo

            bar  : string(filename) !gumph
            baz  : enumeration[ thing1, thing2 ]
            qux  : real
            fred : real[ 'qux * 2' ]

            end namelist foo
            '''))

        uut = NamelistDescriptionAnalyser(database)
        uut.analyse(test_filename)
        dependencies = list(database.getDependencies())
        assert [(Path('foo_configuration_mod.f90'),
                [test_filename])] == dependencies


class TestFortranAnalyser():
    @pytest.fixture
    def database(self, tmp_path_factory):
        filename = tmp_path_factory.mktemp('db-', True) / 'test.db'
        database = SQLiteDatabase(filename)
        return FortranDependencies(database)

    def test_continuation_lines(self, database, tmp_path: Path):
        """
        Ensure continuation lines are handled correctly.
        """
        database.addProgram(u'stock', u'oxo.f90')
        database.addCompileDependency(u'stock', u'widgit')
        database.addLinkDependency(u'stock', u'widgit')
        database.addCompileDependency(u'stock', u'thingy')
        database.addLinkDependency(u'stock', u'thingy')
        database.addModule(u'beef', u'cow.f90')
        database.addModule(u'pork', u'pig.f90')

        test_filename = tmp_path / 'cont.f90'
        test_filename.write_text(dedent('''
            subroutine thingy(cheese, meat, &
                              teapot, fishslice, &
                             )

            end subroutine thingy

            subroutine widgit(grunk, &
                              ! This comment should be ignored
                              grank)
              use beef
              implicit none
            end subroutine widgit

            subroutine old_school(bibble, &
                                & bobble)
              use pork
              implicit none
            end subroutine old_school
            '''))

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        dependencies = list(database.getCompileDependencies())
        assert [(u'old_school', test_filename, u'procedure', u'pork', Path('pig.f90'), u'module'),
                (u'stock', Path('oxo.f90'), u'program', u'thingy', test_filename, u'procedure'),
                (u'stock', Path('oxo.f90'), u'program', u'widgit', test_filename, u'procedure'),
                (u'widgit', test_filename, u'procedure', u'beef', Path('cow.f90'), u'module')] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('widgit'))
        assert [('widgit', test_filename,
                 u'beef', Path('cow.f90'))] == dependencies

    def testProcedure(self, database, tmp_path: Path):
        """
        Procedure as program unit.
        """
        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
           subroutine empty_sub()
           end subroutine empty_sub

           subroutine one_sub(argument)
           end subroutine one_sub

           real function cosd(degree)
           end function cosd
           '''))

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        assert [(u'cosd',      test_filename),
                (u'empty_sub', test_filename),
                (u'one_sub',   test_filename)] \
            == sorted(database.get_program_units())

    def testAnalyseProgram(self, database, tmp_path: Path):
        """
        Includes disparate case to ensure case insensitivity.
        """
        database.addModule('constants_mod', 'constants_mod.f90')
        database.addModule('trumpton_mod', 'trumpton_mod.f90')

        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
            program fOo
              use constAnts_mod, only : i_def
              use trumpton_Mod, only : hew, pew, barney, mcgrey, &
                                       cuthbirt, dibble, grub
              implicit none
            end program fOo
            '''))

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        programs = list(database.getPrograms())
        assert [u'foo'] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'foo', test_filename, u'program',
                 u'constants_mod', Path('constants_mod.f90'), u'module'),
                (u'foo', test_filename, u'program',
                 u'trumpton_mod', Path('trumpton_mod.f90'), u'module')] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('foo'))
        assert [(u'foo', test_filename,
                 u'constants_mod', Path('constants_mod.f90')),
                (u'foo', test_filename,
                 u'trumpton_mod', Path('trumpton_mod.f90'))] \
            == sorted(dependencies)

    def testAnalyseModule(self, database, tmp_path: Path):
        """
        Includes disparate case to ensure case insensitivity.
        """
        uut = FortranAnalyser([], database)

        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
            module foO
              ! Ignore this
              use consTants_mod, only : i_def
              use trumPton_mod, only : hew, pew, barney, mcgrey, &
                                       cuthbirt, dibble, grub
              implicit none
              private
            contains
          end module foO
          module truMpton_mod
          end module truMpton_mod
          '''))
        uut.analyse(test_filename)

        other_filename = tmp_path / 'other.f90'
        other_filename.write_text(dedent('''
            module coNstants_mod
            end module coNstants_mod
            '''))
        uut.analyse(other_filename)

        programs = list(database.getPrograms())
        assert [] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'foo', test_filename, u'module', u'constants_mod', other_filename, u'module'),
                (u'foo', test_filename, u'module', u'trumpton_mod', test_filename, u'module')] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('foo'))
        assert [(u'foo', test_filename, u'constants_mod', other_filename),
                (u'foo', test_filename, u'trumpton_mod', test_filename)] \
            == sorted(dependencies)

    def testAnalyseSubModule(self, database, tmp_path):
        """
        This test also includes disparate case to ensure case insensitivity is
        enforced.
        """
        uut = FortranAnalyser([], database)

        parent_filename = tmp_path / 'parent.f90'
        parent_filename.write_text(dedent('''
            module Parent
              implicit none
              private
              type, public :: test_type
              contains
                procedure foo
                procedure bar
                procedure baz
              end type test_type
              interface
                module subroutine foo(this, cheese)
                  class(test_type), intent(inout) :: this
                  real,             intent(in)    :: cheese
                end subroutine foo
                module subroutine bar(this, teapot)
                  class(test_type), intent(inout) :: this
                  character(*),     intent(in)    :: teapot
                end subroutine bar
                type(something) module function baz(this)
                  class(test_type), intent(in) :: this
                end function baz
              end interface
            end module Parent
            submodule (pArent) chIld3
              implicit none
            contains
              type(something) module function baz(this)
                class(test_type), intent(in) :: this
              end function baz
            end submodule chIld3
            '''))
        uut.analyse(parent_filename)

        child1_filename = tmp_path / 'child1.f90'
        child1_filename.write_text(dedent('''
            submodule(paRent) Child1
              implicit none
              type :: secondary_type
              contains
                procedure baz
              end type secondary_type
              interface
                module subroutine baz(this, wibble)
                  class(secondary_type), intent(inout) :: this
                  real,                  intent(in)    :: wibble
                end subroutine baz
              end interface
            contains
              module subroutine foo(this, cheese)
                implicit none
                class(test_type), intent(inout) :: this
                real,             intent(in)    :: cheese
                type(secondary_type) :: thang
                thang = secondary_type()
                write(6, *) cheese
                call thang%baz(12.7)
              end subroutine foo
            end submodule Child1
            '''))
        uut.analyse(child1_filename)

        child2_filename = tmp_path / 'child2.f90'
        child2_filename.write_text(dedent('''
            submodule (parent) cHild2
              implicit none
            contains
              module procedure bar
                implicit none
                write(6, *) teapot
              end procedure bar
            end submodule cHild2
            '''))
        uut.analyse(child2_filename)

        child3_filename = tmp_path / 'child3.f90'
        child3_filename.write_text(dedent('''
            submodule (parEnt:chilD1) grandChild
              implicit none
            contains
              module subroutine baz(this, wibble)
                implicit none
                class(secondary_type), intent(inout) :: this
                real,                  intent(in)    :: wibble
                write(6, *) wibble
              end subroutine baz
            end submodule grandChild
            '''))
        uut.analyse(child3_filename)

        programs = list(database.getPrograms())
        assert [] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'child1', child1_filename, u'module', u'parent', parent_filename, u'module'),
                (u'child2', child2_filename, u'module', u'parent', parent_filename, u'module'),
                (u'grandchild', child3_filename, u'module', u'child1', child1_filename, u'module')] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('parent'))
        assert [(u'parent', parent_filename, u'child1', child1_filename),
                (u'parent', parent_filename, u'child2', child2_filename)] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('child1'))
        assert [(u'child1', child1_filename, u'grandchild', child3_filename)] \
            == sorted(dependencies)

    def testFunctionInModuleName(self, database, tmp_path: Path):
        """
        Ensure the analyser isn't tripped up by naked global level procedures
        as program units.
        """
        uut = FortranAnalyser([], database)

        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
            module function_thing_mod
              use constants_mod, only : i_def
              implicit none
              private
           contains
           end module function_thing_mod
           '''))
        uut.analyse(test_filename)

        other_filename = tmp_path / 'other.f90'
        other_filename.write_text(dedent('''
            module constants_mod
            end module constants_mod
            '''))
        uut.analyse(other_filename)

        depend_filename = tmp_path / 'dependson.f90'
        depend_filename.write_text(dedent('''
            subroutine dependson
            end dependson
            '''))

        programs = list(database.getPrograms())
        assert [] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'function_thing_mod', test_filename, u'module',
                 u'constants_mod', other_filename, u'module')] == sorted(dependencies)

        dependencies = list(database
                            .getLinkDependencies('function_thing_mod'))
        assert [(u'function_thing_mod', test_filename,
                 u'constants_mod', other_filename)] == sorted(dependencies)

    def testDependsOn(self, database, tmp_path: Path):
        '''
        The analyser has to be able to track dependencies using the deprecated
        "depends on:" comments of the UM.
        '''
        database.addProcedure(u'flibble', u'flibble.f90')
        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
            module function_thing_mod

              use constants_mod, only : i_def

              implicit none

              ! Add in an interface block - this will test to make sure
              ! we don't pick up a spurious subroutine call
              interface
                subroutine wooble ()
                end subroutine
              end interface

              private

              ! Comments before the "depends on" shouldn't upset it.

              ! depends on: flibble.o

              ! depends on: wooble

            contains

            end module function_thing_mod
            '''))

        other_filename = tmp_path / 'other.f90'
        other_filename.write_text(dedent('''
            module constants_mod
            contains
              subroutine wooble
              end subroutine wooble
            end module constants_mod
            '''))

        depend_filename = tmp_path / 'wooble.f90'
        depend_filename.write_text(dedent('''
            subroutine wooble
            end subroutine wooble
            '''))

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)
        uut.analyse(other_filename)
        uut.analyse(depend_filename)

        programs = list(database.getPrograms())
        assert [] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'function_thing_mod', test_filename, u'module',
                 u'constants_mod', other_filename, u'module'),
                (u'function_thing_mod', test_filename, u'module',
                 u'flibble', Path('flibble.f90'), u'procedure'),
                (u'function_thing_mod', test_filename, u'module',
                 u'wooble', depend_filename, u'procedure')] == sorted(dependencies)

        dependencies = list(database
                            .getLinkDependencies('function_thing_mod'))
        assert [(u'function_thing_mod', test_filename,
                 u'constants_mod', other_filename),
                (u'function_thing_mod', test_filename,
                 u'flibble', Path('flibble.f90')),
                (u'function_thing_mod', test_filename,
                 u'wooble', depend_filename)] == sorted(dependencies)

    def testAbstractInterface(self, database, tmp_path: Path):
        """
        The analyser must ignore abstract interface definitions. These are not
        program units.
        """
        first_filename = tmp_path / 'test.f90'
        first_filename.write_text(dedent('''
            module first_mod

              implicit none

              private

              abstract interface
                subroutine thing_face()
                  implicit none
                end subroutine thing_face
              end interface

            contains

            end module first_mod
            '''))

        second_filename = tmp_path / 'test2.f90'
        second_filename.write_text(dedent('''
            module second_mod

              implicit none

              private

              abstract interface
                subroutine thing_face()
                  implicit none
                end subroutine thing_face
              end interface

            contains

          end module second_mod
          '''))

        uut = FortranAnalyser([], database)
        uut.analyse(first_filename)
        uut.analyse(second_filename)

    def testExternal(self, database, tmp_path):
        """
        Ensure "external" works as a dependency marker.
        """
        database.addModule('wibble', 'wibble.f90')
        database.addModule('bibble', 'bibble.f90')
        database.addModule('ibble', 'ibble.f90')
        database.addModule('gribble', 'gribble.f90')
        test_filename = tmp_path / 'test.f90'
        test_filename.write_text(dedent('''
           program boo

             implicit none

             external ibble
             external wibble, bibble, gribble

             call wibble()
             call bibble()

           end program boo
           '''))

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        programs = list(database.getPrograms())
        assert [u'boo'] == programs

        dependencies = list(database.getCompileDependencies())
        assert [(u'boo', test_filename, u'program', u'bibble', Path('bibble.f90'), u'module'),
                (u'boo', test_filename, u'program', u'gribble', Path('gribble.f90'), u'module'),
                (u'boo', test_filename, u'program', u'ibble', Path('ibble.f90'), u'module'),
                (u'boo', test_filename, u'program', u'wibble', Path('wibble.f90'), u'module')] \
            == sorted(dependencies)

        dependencies = list(database.getLinkDependencies('boo'))
        assert [(u'boo', test_filename, u'bibble', Path('bibble.f90')),
                (u'boo', test_filename, u'gribble', Path('gribble.f90')),
                (u'boo', test_filename, u'ibble', Path('ibble.f90')),
                (u'boo', test_filename, u'wibble', Path('wibble.f90'))] \
            == sorted(dependencies)
