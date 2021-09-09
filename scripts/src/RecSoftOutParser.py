#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue June  5 10:18:00 2018
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

DEFAULT_ERROR_VALUE = -1

def extractFilename(str):
    """
    Returns the filename only from a path.
    :param str: path to file
    :return:
    """
    head, tail = os.path.split(str)
    return tail

def isConvertible(value):
    """
    Check if str can match with int or float.
    """
    if(isConvertibleToFloat(value)):
        return True

    try:
        int(value)
        return True
    except:
        return False


def isConvertibleToFloat(value):
    """
    Check if str matches with int or float.
    """
    try:
        float(value)
        return True
    except:
        return False

def strToNum(str):
    """
    Returns a number in string to numeric (float or int)
    :param str:
    :return: float or int
    """
    if not isConvertible(str):
        print "Invalid read: " + str
        return DEFAULT_ERROR_VALUE

    if (isConvertibleToFloat(str)):
        return float(str)
    else:
        return int(str)

def getTreerecsOutputScore(filename):
    """
        Returns total cost of a treerecs solution.
    """
    file = open(filename, "r")

    line = file.readline()

    file.close()

    # first: find the total cost field
    request = "total cost = "

    pos = line.find(request)

    pos += len(request)

    # then extract associated value
    totalCost_str_begin = pos

    while line[pos] != "," and pos < len(line):
        pos += 1

    totalCost_str_end = pos

    totalCost_str = line[totalCost_str_begin: totalCost_str_end]

    # and return value as float
    return strToNum(totalCost_str)


def getProfileNJOutputScore(filename):
    """
        Returns total cost of a profilenj solution.
    """
    file = open(filename, "r")

    line = file.readline()

    file.close()

    # first: find the total cost field
    request = "m_cost="

    pos = line.find(request)

    pos += len(request)

    # then extract associated value and return as float
    totalCost_str = line[pos: len(line)]
    return strToNum(totalCost_str)


def getNotungOutputScore(filename):
    """
        Returns total cost of a notung solution.
    """
    file = open(filename, "r")

    lines = file.readlines()
    line = lines[-1]
    if(len(line) == 1):
        # due to the new line in root mode with Notung, this must be take in
        # consideration. So we use the previous line before the end
        line = lines[-2]

    file.close()

    # first: find the total cost field
    request = "core is "

    pos = line.find(request)

    pos += len(request)

    # then extract associated value and return as float
    totalCost_str = line[pos: len(line)]
    return strToNum(totalCost_str)

def getEcceteraOutputScore(filename):
    """
        Returns total cost of a EcceTERA solution.
    """
    file = open(filename, "r")

    lines = file.readlines()
    line = lines[-1]

    file.close()

    # first: find the total cost field
    request = ": "

    pos = line.find(request)

    pos += len(request)

    # then extract associated value and return as float
    totalCost_str = line[pos: len(line)]

    return strToNum(totalCost_str)

def getRangerOutputScore(filename):
    """
        Returns total cost of a Ranger-DTL solution.
    """
    file = open(filename, "r")

    lines = file.readlines()
    line = ""
    line_temp = ""
    for line_temp in lines:
        if line_temp.find("inimum ") != -1:
            line = line_temp

    file.close()

    # first: find the total cost field
    request = ": "

    pos = line.find(request)

    pos += len(request)

    # then extract associated value and return as float
    totalCost_str = line[pos: len(line)]

    if totalCost_str == "":
        return DEFAULT_ERROR_VALUE

    print

    return strToNum(totalCost_str)

