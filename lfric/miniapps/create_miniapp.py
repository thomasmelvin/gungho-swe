#!/usr/bin/env python

##############################################################################
# (c) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

""" Create a new miniapp based on the existing skeleton miniapp.
This will replace instances of "skeleton" with the name provided
in the command line arguments."""
from __future__ import absolute_import
from __future__ import print_function
import argparse
import re
import os
from shutil import copyfile

skeleton_word = "skeleton"


def replace_keep_case(word, replacement, text):
    """ Wrapper to re.sub to replace a word while maintaining the case"""
    def func(match):
        """Function converts replacement string to have same case as match"""
        grp = match.group()
        if grp.islower():
            return replacement.lower()
        if grp.istitle():
            return replacement.title()
        if grp.isupper():
            return replacement.upper()
        return replacement
    return re.sub(word, func, text, flags=re.I)


def run(miniapp_name, miniapp_dir):
    """ Create a new miniapp based on the existing skeleton miniapp """

    miniapp_var = re.sub(pattern="-", repl="_", string=miniapp_name)
    miniapp_root = "{}/{}".format(miniapp_dir, miniapp_name)
    skeleton_root = "{}/{}".format(miniapp_dir, skeleton_word)

    if os.path.exists(miniapp_root):
        raise ValueError("Miniapp already exists: {}".format(miniapp_root))
    else:
        os.makedirs(miniapp_root)

    # Traverse skeleton root directory, and create the corresponding
    # files and directories for the newly named miniapp
    for root, dirs, files in os.walk(skeleton_root):

        relpath = os.path.relpath(root, skeleton_root)
        newpath = replace_keep_case(skeleton_word, miniapp_var, relpath)

        for each_dir in dirs:
            if os.path.islink(root+"/"+each_dir):
                # If dir is a symlink, copy it (assumes a relative symlink)
                olddir = "{}/{}".format(root, each_dir)
                newdir = "{}/{}/{}".format(miniapp_root, newpath, each_dir)
                copyfile(olddir, newdir, follow_symlinks=False)
            else:
                # Create new directories (substituting name where required)
                newdir = replace_keep_case(skeleton_word, miniapp_var,
                                           each_dir)
                dir_out = "{}/{}/{}".format(miniapp_root, newpath, newdir)
                os.makedirs(dir_out)

        for each_file in files:
            # Create new files (replacing name in text where required)
            file_out = replace_keep_case(skeleton_word, miniapp_var, each_file)
            oldfile = "{}/{}".format(root, each_file)
            newfile = "{}/{}/{}".format(miniapp_root, newpath, file_out)
            if re.search(r"\.nc", each_file):
              copyfile(oldfile, newfile)
            else:
              f_in = open(oldfile)
              f_out = open(newfile, "w")
              for line in f_in.readlines():
                  f_out.write(replace_keep_case(skeleton_word,
                                                miniapp_var, line))
              # Close the files
              f_out.close()
              f_in.close()

    print("""
Successfully created the new miniapp: {}
Please now look at the miniapp documentation and populate it with
some functionality.

Note that if you wish your miniapp to go into the repository, you
should now run
 > fcm add {}

(Do it now before building or adding other files in there!)
""".format(miniapp_name, miniapp_root))

if __name__ == "__main__":

    # Parse the command line arguments
    PARSER = argparse.ArgumentParser(
        description="""
description: Makes a copy of the skeleton miniapp by replacing
instances of 'skeleton' with the input 'name'.""")
    PARSER.add_argument("name", help="Name of the new miniapp")
    PARSER.add_argument("-d", "--miniapp-directory",
                        default=(os.path.dirname(os.path.realpath(__file__))),
                        help="""Directory in which the miniapps are stored.
Default is the location of this script""")
    ARGS = PARSER.parse_args()

    run(ARGS.name, ARGS.miniapp_directory)
