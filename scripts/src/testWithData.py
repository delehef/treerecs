#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue June  5 11:00:00 2018
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

from RecSoftComposer import*
from RecSoftOutParser import*

def saveProfileNJRes(output_filename, genetree_filename, support_threshold, rerooted, result):
    """
    Append results in output_filename.

    :param output_filename: file to write data.
    :param genetree_filename: genetree input file used to save data.
    :param support_threshold: branch support threshold used to solve genetrees.
    :param rerooted: boolean value, true if the reroot option has been used, else false.
    :param result: total cost of the resulting reconciliation
    :return: nothing
    """
    file = open(output_filename, "a")
    rerooted_str = ""
    if (rerooted):
        rerooted_str += "true"
    else:
        rerooted_str += "false"

    file.write(str(genetree_filename) + "\t"
               + str(support_threshold) + "\t"
               + str(rerooted_str) + "\t"
               + str(result) + "\n")

    file.close()

def getSavedProfileNJRes(output_filename, genetree_filename, support_threshold, rerooted):
    """
    Return score saved in file according to parameters
    :param output_filename:
    :param genetree_filename:
    :param support_threshold:
    :param rerooted:
    :return:
    """
    if(os.path.exists(output_filename)):
        file = open(output_filename, "r")
        rerooted_str = ""
        if (rerooted):
            rerooted_str += "true"
        else:
            rerooted_str += "false"

        for line in file:
            line_genetree,line_support_threshold,line_rerooted_str,result = line.split("\t")
            if (line_genetree == genetree_filename
                and line_support_threshold == str(support_threshold)
                and line_rerooted_str == rerooted_str):
                return result

        file.close()

    return ""

def compareTreerecsAndProfileNJWithPhymlTrees(val_parse):
    """
    If there is no genetree given in command line argument, this script will
    test 13128 gene trees in phyml trees randomly according to parameters in command
    line.

    If ProfileNJ has already been ran with a particular genetree and parameters
    before, this script will use these saved results instead of run again
    ProfileNJ. As a consequence, there will be no execution time comparisons
    between Treerecs and ProfileNJ.

    :param val_parse: command line arguments.
    :return:
    """

    ensembl_compara_73_path = "../examples/ensembl_compara_73/"

    if not os.path.exists(ensembl_compara_73_path):
        sys.stdout.write("Error: " + ensembl_compara_73_path + " does not exist.\n")
        sys.exit(1)

    treerecs_output_filename = "temp_treerecs_output.txt"

    profilenj_output_filename = "temp_profilenj_output.txt"

    savefile_used = False

    speciestree_file = ensembl_compara_73_path + "species_tree.txt"
    smap_file = ensembl_compara_73_path + "phyml_trees.smap"

    # Split phyml_trees in multiple files and move them into a new folder "temp_phyml_trees_split"
    sys.stdout.flush()
    sys.stdout.write("\rPreparing files...")
    evaluateCommand("rm -rf temp_phyml_trees_split")
    evaluateCommand("mkdir temp_phyml_trees_split")
    evaluateCommand(
        "awk '{print $0 \";\"> \"temp_phyml_trees_split_\" NR}' RS=';' "
        + ensembl_compara_73_path + "phyml_trees.txt")
    evaluateCommand("mv temp_phyml_trees_split_* temp_phyml_trees_split")
    evaluateCommand("rm temp_phyml_trees_split/temp_phyml_trees_split_13129")

    recordedExecTime_treerecs = []
    recordedExecTime_profilenj = []

    recordedTotalCost_treerecs = []
    recordedTotalCost_profilenj = []

    # Then evaluate each genetree
    path, dirs, files = os.walk("temp_phyml_trees_split").next()

    sorted(files)

    iteration = 0
    number_of_trees = val_parse["repetition"]

    # Do a random use of files
    #shuffle(files)

    for genetree_file in files[0:number_of_trees]:
        genetree_file = "./temp_phyml_trees_split/" + genetree_file
        # Eliminate all endl character
        evaluateCommand(
            "cat " + genetree_file + " | tr -d \"" + os.linesep + "\" > " + genetree_file + "_strip")
        evaluateCommand("mv " + genetree_file + "_strip " + genetree_file)

    for genetree_file in files[0:number_of_trees]:
        genetree_file = "./temp_phyml_trees_split/" + genetree_file

        sys.stdout.flush()
        sys.stdout.write("\rProgression: " + str(
            int(100 * (iteration) / number_of_trees)) + "%.")

        # Generate commands
        treerecs_command = createTreerecsCommand_(val_parse["treerecs"],
                                                 genetree_file,
                                                 speciestree_file,
                                                 contraction_threshold=
                                                 val_parse["threshold"],
                                                 reroot=val_parse["reroot"],
                                                 smap=smap_file,
                                                 output_filename=treerecs_output_filename,
                                                 parallelize=val_parse[
                                                     "parallelize"])[0]

        # Evaluate Treerecs command
        treerecs_res = evaluateCommand(treerecs_command)

        # Save execution time
        recordedExecTime_treerecs.append(treerecs_res[1])

        # Get total costs computed with treerecs
        treerecs_total_cost = getTreerecsOutputScore(
            treerecs_output_filename + ".nwk")

        # Save score.
        recordedTotalCost_treerecs.append(treerecs_total_cost)

        sys.stdout.flush()
        sys.stdout.write("\rProgression: " + str(
            int(100 * (iteration + 0.5) / number_of_trees)) + "%.")

        # Evaluate with ProfileNJ
        # First try with a file which contains saved results
        profilenj_saved_res = getSavedProfileNJRes(val_parse["resfile"],
                                                   genetree_file,
                                                   val_parse["threshold"],
                                                   val_parse["reroot"])

        # if there is no data, run ProfileNJ
        if(profilenj_saved_res == ""):
            # If profileNJ has not been use and saved with its results before,
            # run it.
            profilenj_command = createProfileNJCommand_(val_parse["profilenj"],
                                                       genetree_file,
                                                       speciestree_file,
                                                       contraction_threshold=
                                                       val_parse["threshold"],
                                                       reroot=val_parse[
                                                           "reroot"],
                                                       smap=smap_file,
                                                       output_filename=profilenj_output_filename,
                                                       parallelize=val_parse[
                                                           "parallelize"])[0]
            profilenj_res = evaluateCommand(profilenj_command)

        if(profilenj_saved_res == ""):
            recordedExecTime_profilenj.append(profilenj_res[1])
        else :
            recordedExecTime_profilenj.append(0.0)
            savefile_used = True

        # If results have already been computed for these parameters, use them
        # or extract them from a ProfileNJ's run output.
        if(profilenj_saved_res != ""):
            profilenj_total_cost = profilenj_saved_res

        else :
            profilenj_total_cost = getProfileNJOutputScore(
                profilenj_output_filename)
            saveProfileNJRes(val_parse["resfile"],
                             genetree_file,
                             val_parse["threshold"],
                             val_parse["reroot"],
                             profilenj_total_cost)

        recordedTotalCost_profilenj.append(profilenj_total_cost)

        # Then, delete temp files
        evaluateCommand("rm " + treerecs_output_filename + ".newick")
        if(profilenj_saved_res == ""):
            evaluateCommand("rm " + profilenj_output_filename)

        iteration += 1

    sys.stdout.flush()
    sys.stdout.write("\rProgression: " + str(
        int(100 * (iteration) / number_of_trees)) + "%." + os.linesep)

    # Compute average execution time and control similarity in results.
    execTimeAverage_treerecs = 0.0
    execTimeAverage_profilenj = 0.0
    for i in range(len(recordedExecTime_treerecs)):
        execTimeAverage_treerecs += recordedExecTime_treerecs[i]
        execTimeAverage_profilenj += recordedExecTime_profilenj[i]

        if (float(recordedTotalCost_treerecs[i]) !=
                float(recordedTotalCost_profilenj[i])):
            print "In genetree " + str(i + 1) + "/" + str(
                len(recordedExecTime_treerecs)) + ":"
            print "Treerecs and ProfileNJ have a different result: Treerecs has " + str(
                recordedTotalCost_treerecs[
                    i]) + " and ProfileNJ has " + str(
                recordedTotalCost_profilenj[i])

    if(savefile_used == False): #if profileNJ has been ran for all trees (
        # without a use of a save of results) we can compare execution times
        # between ProfileNJ and Treerecs.
        execTimeAverage_treerecs /= len(recordedExecTime_treerecs)
        execTimeAverage_profilenj /= len(recordedExecTime_profilenj)

        print "Treerecs has an execution time average at " + str(
            execTimeAverage_treerecs) + "."
        print "ProfileNJ has an execution time average at " + str(
            execTimeAverage_profilenj) + "."
        print "Global time ratio = " + str(
            execTimeAverage_profilenj / execTimeAverage_treerecs) + "."

    # Then rm all files
    evaluateCommand("rm -r ./temp_phyml_trees_split")
