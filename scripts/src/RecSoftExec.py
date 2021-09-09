#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue June  5 10:16:00 2018
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

import os
import sys
from random import shuffle

from RecSoftComposer import*
from RecSoftOutParser import*
from CLExec import *

def compareReconciliationSoftwaresWithOneTree_mode(val_parse):
    """
    If a gene tree is given in command line argument, this program will execute
    different softwares in a random order.
    :param val_parse: command line arguments.
    :return:
    """

    # Test data
    if (not os.path.isfile(val_parse["genetree"])):
        print "genetree " + val_parse["genetree"] + " does_not exist."
        sys.exit(1)

    if (not os.path.isfile(val_parse["speciestree"])):
        print "speciestree " + val_parse["speciestree"] + " does_not exist."
        sys.exit(1)

    if (val_parse["smap"] != "" and not os.path.isfile(
            val_parse["genetree"])):
        print "Smap " + val_parse["smap"] + " does_not exist."
        sys.exit(1)

    # Create commands
    commands = []
    output_files = []

    if val_parse["treerecs"] != "":
        treerecs_command = createTreerecsCommand(val_parse)

        commands.append(treerecs_command[0])
        output_files += treerecs_command[1]

    if val_parse["profilenj"] != "":
        profilenj_command = createProfileNJCommand(val_parse)

        commands.append(profilenj_command[0])
        output_files += profilenj_command[1]

    if val_parse["notung"] != "":
        if val_parse["reroot"] and val_parse["threshold"] > 0.0 :
            sys.stdout.write("\rWarning: Notung cannot root and rearrange in the"
                             " same time. It will only root the tree.\n")
        if val_parse["smap"] != "":
            sys.stdout.write("\rWarning: Notung does not support SMAP files. Please use"
                             "gene names with species name in prefix separated with '_'.\n")

        notung_command = createNotungCommand(val_parse)
        commands.append(notung_command[0])
        output_files += notung_command[1]

    if(val_parse["eccetera"] != ""):
        eccetera_command = createEcceteraCommand(val_parse)
        commands.append(eccetera_command[0])
        output_files += eccetera_command[1]

    if(val_parse["ranger"] != ""):
        if val_parse["reroot"] and val_parse["threshold"] > 0.0 :
            sys.stdout.write("\rWarning: Ranger cannot root and rearrange in the"
                             " same time. It will only root the tree.\n")

        if val_parse["smap"] != "":
            sys.stdout.write("\rWarning: Ranger does not support SMAP files. Please use"
                             "gene names with species name in prefix separated with '_'.\n")

        ranger_command = createRangerCommand(val_parse)
        commands.append(ranger_command[0])
        output_files += ranger_command[1]

    total_exec_times = {}

    for command in commands:
        total_exec_times[command] = 0.0


    program_turn = commands * val_parse["repetition"]
    shuffle(program_turn)

    # Evaluate commands
    if(val_parse["verbose"]):
        sys.stdout.write("\rEvaluate commands:\n")

    iteration = 0
    for p in program_turn:
        if(not val_parse["quiet"]):
            sys.stdout.flush()
            sys.stdout.write("\rProgression: " + str(
                int(100 * (iteration) / len(program_turn))) + "%.")

        if(val_parse["verbose"]):
            sys.stdout.write("\r# " + p + "\n")

        command_res = evaluateCommand(p)
        total_exec_times[p] += command_res[1]

        iteration += 1

    if(not val_parse["quiet"]):
        sys.stdout.write("\rProgression: " + str(
            int(100 * (iteration) / len(program_turn))) + "%." + os.linesep)

    for command in commands:
        total_exec_times[command] /= val_parse["repetition"]
        sys.stdout.write("$ " + shortCommand(command)[0:43] + "[...]\n\texecution time average = "
                         + str(total_exec_times[command]))
        sys.stdout.write("\n")

    # Get results
    treerecs_total_cost = DEFAULT_ERROR_VALUE
    profilenj_total_cost = DEFAULT_ERROR_VALUE
    notung_total_cost = DEFAULT_ERROR_VALUE
    eccetera_total_cost = DEFAULT_ERROR_VALUE
    ranger_total_cost = DEFAULT_ERROR_VALUE

    if(val_parse["treerecs"] != ""):
        treerecs_total_cost = getTreerecsOutputScore(DEFAULT_TREERECS_OUTPUT_FILENAME + ".nwk")

    if(val_parse["profilenj"] != ""):
        profilenj_total_cost = getProfileNJOutputScore(DEFAULT_PROFILENJ_OUTPUT_FILENAME)

    if(val_parse["notung"] != ""):
        notung_total_cost = getNotungOutputScore(DEFAULT_NOTUNG_OUTPUT_FILENAME)

    if(val_parse["eccetera"] != ""):
        eccetera_total_cost = getEcceteraOutputScore(DEFAULT_ECCETERA_OUTPUT_FILENAME)

    if(val_parse["ranger"] != ""):
        ranger_total_cost = getRangerOutputScore(DEFAULT_RANGER_OUTPUT_FILENAME)

    # Print results
    file_to_keep = []
    if(treerecs_total_cost != DEFAULT_ERROR_VALUE):
        print "Treerecs result = " + str(treerecs_total_cost)
    elif val_parse["treerecs"] != "" :
        file_to_keep += treerecs_command[1]

    if(profilenj_total_cost != DEFAULT_ERROR_VALUE):
        print "ProfileNJ result = " + str(profilenj_total_cost)
    elif val_parse["profilenj"] != "" :
        file_to_keep += profilenj_command[1]

    if(notung_total_cost != DEFAULT_ERROR_VALUE):
        print "Notung result = " + str(notung_total_cost)
    elif val_parse["notung"] != "" :
        file_to_keep += notung_command[1]

    if(eccetera_total_cost != DEFAULT_ERROR_VALUE):
        print "EcceTERA result = " + str(eccetera_total_cost)
    elif val_parse["eccetera"] != "":
        file_to_keep += eccetera_command[1]

    if(ranger_total_cost != DEFAULT_ERROR_VALUE):
        print "Ranger-DTL result = " + str(ranger_total_cost)
    elif val_parse["ranger"] != "" :
        file_to_keep += ranger_command[1]

    # Then, delete temp files
    if not val_parse["keepfiles"]:
        for file_to_remove in output_files:
            if not file_to_remove in file_to_keep:
                evaluateCommand("rm " + file_to_remove)

    if len(file_to_keep) > 0 :
        print "Un-expected behaviour(s), several files are kept:\n"
        for f in file_to_keep:
            print "\t* " + f + "\n"

    print "Commands used:"

    for command in commands:
        print "$ " + command + "\n"
