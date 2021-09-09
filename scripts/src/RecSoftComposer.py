#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue June  5 10:09:00 2018
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
from CLExec import *

DEFAULT_TREERECS_OUTPUT_FILENAME = "temp_treerecs_output.txt"

DEFAULT_PROFILENJ_OUTPUT_FILENAME = "temp_profilenj_output.txt"

DEFAULT_NOTUNG_OUTPUT_FILENAME = "temp_notung_output.txt"

DEFAULT_ECCETERA_OUTPUT_FILENAME = "temp_eccetera_output.txt"

DEFAULT_RANGER_OUTPUT_FILENAME = "temp_ranger_output.txt"

def shortCommand(cmd_str):
    cmd_ll = cmd_str.split(" ")
    prog_path, prog_name = os.path.split(cmd_ll[0])
    return prog_name + " " + " ".join(cmd_ll[1:])

def createTreerecsCommand_(treerecs_path,
                          genetree_filename,
                          speciestree_filename,
                          contraction_threshold,
                          output_filename=DEFAULT_TREERECS_OUTPUT_FILENAME,
                          dupcost=2.0,
                          losscost=1.0,
                          reroot=False,
                          parallelize=False,
                          sample_size=1,
                          smap=""
                          ):
    """
    Creates Treerecs command. Returns a pair with first the command and as
    second, the list of created files -> (str, [str]).
    """
    treerecs_command = treerecs_path
    treerecs_command += " -g " + genetree_filename
    treerecs_command += " -s " + speciestree_filename
    treerecs_command += " -t " + str(contraction_threshold)
    treerecs_command += " -d " + str(dupcost)
    treerecs_command += " -l " + str(losscost)
    treerecs_command += " -n " + str(sample_size)
    treerecs_command += " -o " + output_filename
    treerecs_command += " --inferior-than-threshold-only -q"

    if (smap != ""):
        treerecs_command += " -S " + smap

    if (reroot):
        treerecs_command += " -r"

    if (parallelize):
        treerecs_command += " -P"

    return (treerecs_command, [output_filename + ".nwk"])

def createTreerecsCommand(val_parse):
    return createTreerecsCommand_(val_parse["treerecs"],
                                 val_parse["genetree"],
                                 val_parse["speciestree"],
                                 contraction_threshold=val_parse["threshold"],
                                 reroot=val_parse["reroot"],
                                 smap=val_parse["smap"],
                                 output_filename=DEFAULT_TREERECS_OUTPUT_FILENAME,
                                 parallelize=val_parse["parallelize"],
                                 dupcost=val_parse["dcost"],
                                 losscost=val_parse["lcost"])

def createProfileNJCommand_(profilenj_path,
                           genetree_filename,
                           speciestree_filename,
                           contraction_threshold,
                           output_filename=DEFAULT_PROFILENJ_OUTPUT_FILENAME,
                           dupcost=2.0,
                           losscost=1.0,
                           reroot=False,
                           parallelize=False,
                           sample_size=1,
                           smap=""
                           ):
    """
    Creates ProfileNJ command. Returns a pair with first the command and as
    second, the list of created files -> (str, [str]).
    """
    profilenj_command = profilenj_path
    profilenj_command += " -g " + genetree_filename
    profilenj_command += " -s " + speciestree_filename
    profilenj_command += " --seuil " + str(contraction_threshold)
    profilenj_command += " --cost " + str(dupcost)
    profilenj_command += " " + str(losscost)
    profilenj_command += " --slimit " + str(sample_size)
    profilenj_command += " -o " + output_filename

    if smap != "":
        profilenj_command += " -S " + smap
    else:
        profilenj_command += " --sep _ --spos prefix"

    if (reroot):
        profilenj_command += " -r best"

    if (parallelize):
        profilenj_command += " --parallelize"

    return (profilenj_command,[output_filename])

def createProfileNJCommand(val_parse):
    return createProfileNJCommand_(val_parse["profilenj"],
                                  val_parse["genetree"],
                                  val_parse["speciestree"],
                                  contraction_threshold=val_parse["threshold"],
                                  reroot=val_parse["reroot"],
                                  smap=val_parse["smap"],
                                  output_filename=DEFAULT_PROFILENJ_OUTPUT_FILENAME,
                                  parallelize=val_parse["parallelize"],
                                  dupcost=val_parse["dcost"],
                                  losscost=val_parse["lcost"])

def createRangerCommand_(ranger_path,
                           genetree_filename,
                           speciestree_filename,
                           contraction_threshold,
                           output_filename=DEFAULT_RANGER_OUTPUT_FILENAME,
                           dupcost=2.0,
                           losscost=1.0,
                           trancost=3.0,
                           reroot=False,
                           sample_size=1,
                           ):
    """
    Creates Ranger-DTL command. Returns a pair with first the command and as
    second, the list of created files -> (str, [str]).
    """

    path, ranger_name = os.path.split(ranger_path)

    if(ranger_name == "Ranger-DTL.linux"):
        ranger_command = path + "/../"
    else:
        ranger_command = ranger_path

    if reroot:
        ranger_command += "/CorePrograms/OptRoot.linux"
    else:
        ranger_command += "/SupplementaryPrograms/OptResolutions.linux"

    temp_ranger_input_name = output_filename + ".temp"

    evaluateCommand("cat " + speciestree_filename + " " + genetree_filename
                    + " > " + temp_ranger_input_name)

    ranger_command += " -i " + temp_ranger_input_name
    ranger_command += " -o " + output_filename

    if(not reroot):
        ranger_command += " -B " + str(contraction_threshold)

    ranger_command += " -D " + str(dupcost)
    ranger_command += " -L " + str(losscost)
    if(trancost >= 0):
        ranger_command += " -T " + str(trancost)
    else:
        ranger_command += " -T " + str(dupcost * 50)

    if(not reroot):
        ranger_command += " -N " + str(sample_size)

    #ranger_command += " -q"

    return (ranger_command, [output_filename, temp_ranger_input_name])

def createRangerCommand(val_parse):
    return createRangerCommand_(val_parse["ranger"],
                               val_parse["genetree"],
                               val_parse["speciestree"],
                               val_parse["threshold"],
                               DEFAULT_RANGER_OUTPUT_FILENAME,
                               reroot=val_parse["reroot"],
                               trancost=val_parse["tcost"],
                               dupcost=val_parse["dcost"],
                               losscost=val_parse["lcost"])

def createNotungCommand_(notung_path,
                        genetree_filename,
                        speciestree_filename,
                        contraction_threshold,
                        output_filename=DEFAULT_NOTUNG_OUTPUT_FILENAME,
                        dupcost=2.0,
                        losscost=1.0,
                        reroot=False,
                        sample_size=1
                        ):
    """
    Creates Notung command. Returns a pair with first the command and as
    second, the list of created files -> (str, [str]).
    """
    notung_command = "java -jar " + notung_path
    notung_command += " -g " + genetree_filename
    notung_command += " -s " + speciestree_filename
    notung_command += " --threshold " + str(contraction_threshold)
    notung_command += " --costdup " + str(dupcost)
    notung_command += " --costloss " + str(losscost)
    notung_command += " --maxtrees " + str(sample_size)

    other_output_file = os.path.split(genetree_filename)[1]

    if (reroot):
        notung_command += " --root"
        other_output_file += ".rooting.0"
    else:
        notung_command += " --rearrange"
        other_output_file += ".rearrange.0"

    notung_command += ""

    notung_command += " > " + output_filename

    return (notung_command, [output_filename, other_output_file])

def createNotungCommand(val_parse):
    return createNotungCommand_(val_parse["notung"],
                               val_parse["genetree"],
                               val_parse["speciestree"],
                               val_parse["threshold"],
                               DEFAULT_NOTUNG_OUTPUT_FILENAME,
                               reroot=val_parse["reroot"],
                               dupcost=val_parse["dcost"],
                               losscost=val_parse["lcost"]
                               )

def createEcceteraCommand_(eccetera_path,
                          genetree_filename,
                          speciestree_filename,
                          contraction_threshold,
                          output_filename=DEFAULT_ECCETERA_OUTPUT_FILENAME,
                          dupcost=2.0,
                          losscost=1.0,
                          hgtcost=3.0,
                          reroot=False,
                          smap=""
                          ):
    """
    Creates EcceTERA command. Returns a pair with first the command and as
    second, the list of created files -> (str, [str]).
    """
    eccetera_command = eccetera_path
    eccetera_command += " gene.file=" + genetree_filename
    eccetera_command += " species.file=" + speciestree_filename
    eccetera_command += " dated=0"
    eccetera_command += " collapse.mode=1"
    eccetera_command += " collapse.threshold=" + str(contraction_threshold)
    eccetera_command += " dupli.cost=" + str(dupcost)
    eccetera_command += " loss.cost=" + str(losscost)
    if hgtcost >= 0.0:
        eccetera_command += " HGT.cost=" + str(hgtcost)
    else:
        eccetera_command += " compute.T=false"

    if (smap != ""):
        eccetera_command += " gene.mapping.file=" + smap

    if (reroot):
        eccetera_command += " resolve.trees=1"
    else:
        eccetera_command += " resolve.trees=0"

    eccetera_command += " > " + output_filename

    return (eccetera_command, [output_filename])

def createEcceteraCommand(val_parse):
    return createEcceteraCommand_(val_parse["eccetera"],
                                 val_parse["genetree"],
                                 val_parse["speciestree"],
                                 val_parse["threshold"],
                                 DEFAULT_ECCETERA_OUTPUT_FILENAME,
                                 reroot=val_parse["reroot"],
                                 smap=val_parse["smap"],
                                 hgtcost=val_parse["tcost"],
                                 dupcost=val_parse["dcost"],
                                 losscost=val_parse["lcost"])
