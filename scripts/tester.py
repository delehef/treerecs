#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Nov  9 10:06:45 2017
@author: Nicolas Comte

Treerecs – A tree reconciliation tool

Copyright (C) 2018  INRIA

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys

from src.RecSoftExec import *
from src.testWithData import *

if __name__ == "__main__":
    """
        If a gene tree is not specified in args, the program will use 
        examples/ensembl_compara_73/phyml_trees.txt and run several trees
        from this file to test results. 
    """
    # Parser =======================================
    parser = argparse.ArgumentParser(prog="RecSoftRunner")
    parser.add_argument("-g", "--genetree", default="",
                        help="Gene tree filename in newick format.")
    parser.add_argument("-s", "--speciestree", default="",
                        help="Species tree filename in newick format.")
    parser.add_argument("-t", "--threshold", type=float, default=0.0,
                        help="Contraction threshold.")
    parser.add_argument("-S", "--smap", type=str, default="",
                        help="Gene<>Species map file in Smap format.")
    parser.add_argument("-r", "--reroot",
                        action="store_true", help="Find best root.")
    parser.add_argument("-p", "--parallelize",
                        action="store_true", help="Parallelize treerecs and profilenj.")
    parser.add_argument("-N", "--repetition", type=int, default=1,
                        help="Number of repetition")
    parser.add_argument("-F", "--resfile", type=str, default="resfile.txt",
                        help="Name of the file containing ProfileNJ results")
    parser.add_argument("--tcost", type=float, default=-1, help="Allow hgt and set"
                                                             "a cost value.")
    parser.add_argument("-d", "--dcost", type=float, default=2.0,
                        help="Set the gene duplication cost.")
    parser.add_argument("-l", "--lcost", type=float, default=1.0,
                        help="Set the gene loss cost.")

    # Path parser ==================================
    parser.add_argument("--treerecs", type=str, default="",
                        help="Path to Treerecs.")
    parser.add_argument("--genetreeeditor", type=str, default="",
                        help="Path to geneTreeEditor")
    parser.add_argument("--profilenj", type=str, default="",
                        help="Path to profileNJ.")
    parser.add_argument("--notung", type=str, default="",
                        help="Set the path to Notung exec and use it.")
    parser.add_argument("--eccetera", type=str, default="",
                        help="Set the path to Eccetera exec and use it.")
    parser.add_argument("--ranger", type=str, default="",
                        help="Set the path to Ranger-DTL "
                             "(or main directory path) exec and use it.")
    parser.add_argument("--keepfiles",
                        action="store_true", help="Keep temporary files.")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="Verbose.")
    parser.add_argument("-q", "--quiet",
                        action="store_true",
                        help="Avoid prints of progression bars.")

    val_parse = vars(parser.parse_args())

    # Run ==========================================
    cpath = os.path.dirname(os.path.abspath(__file__))

    if ( val_parse["genetree"] != "" ):
        compareReconciliationSoftwaresWithOneTree_mode(val_parse)

    else:
        missing_program = False;
        val_parse_simplified = val_parse

        val_parse_simplified["genetree"] = "GENETREE"
        val_parse_simplified["speciestree"] = "SPECIESTREE"
        val_parse_simplified["smap"] = "SMAP"

        if ( val_parse["treerecs"] == "" ):
            sys.stdout.write("Please provide the Treerecs exec with --treerecs"
                             " option\n")
            missing_program = True
        else:
            print "Use commands like:"
            print "$ " + (createTreerecsCommand(val_parse_simplified)[0])

        if( val_parse["profilenj"] == ""):
            sys.stdout.write("Please provide the ProfileNJ exec with --profilenj"
                             " option\n")
            missing_program = True
        else:
            print "$ " + (createProfileNJCommand(val_parse_simplified)[0])

        if missing_program:
            sys.exit(1)

        compareTreerecsAndProfileNJWithPhymlTrees(val_parse)


    sys.exit(0)  # success !
