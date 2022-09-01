##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
System test of Fortran Dependency Analyser
"""
import os
import subprocess
from pathlib import Path


class TestFortranDependencyAnalyser():

    def test_dependencies(self, tmp_path: Path):

        database = tmp_path / 'dependencies.db'
        dependencies_file = tmp_path / 'dependencies.mk'
        programs_file = tmp_path / 'programs.mk'
        test_path = Path(__file__)
        source_dir = test_path.parent / 'source'
        tool_dir = test_path.parent.parent
        os.chdir(source_dir)

        source_files = []
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if file.endswith(".f90") or file.endswith(".F90"):
                    rel_dir = os.path.relpath(root, source_dir)
                    rel_file = os.path.join(rel_dir, file)
                    source_files.append(rel_file)

        def build():

            for file in source_files:
                command = ['python', tool_dir / 'DependencyAnalyser',
                           '-verbose',
                           database,
                           file]
                process = subprocess.run(command)
                assert process.returncode == 0

            command = ['python', tool_dir / 'DependencyRules', '-verbose',
                       '-database', database, '-moduledir', 'modules',
                       '-objectdir',
                       'objects', dependencies_file]
            process = subprocess.run(command)
            assert process.returncode == 0

            command = ['python', tool_dir / 'ProgramObjects', '-database',
                       database, '-objectdir', 'objects', programs_file]
            process = subprocess.run(command)
            assert process.returncode == 0

            assert Path('../expected.dependencies.mk').read_text() == \
                   dependencies_file.read_text()
            assert Path('../expected.programs.mk').read_text() == \
                   programs_file.read_text()

        build()

        # rebuild with one changed file
        source_files = [source_files[0]]
        build()
