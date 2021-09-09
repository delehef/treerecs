#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Tue June  5 17:35:00 2018
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

import sys
import time
import subprocess

def syscmd(cmd, encoding=''):
    """
    Runs a command on the system, waits for the command to finish, and then
    returns the text output of the command. If the command produces no text
    output, the command's return code will be returned instead.
    """

    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         close_fds=True)
    p.wait()
    if (p.returncode != 0):
        print p.stdout.read()
        #sys.exit(1)

    return p.returncode

def evaluateCommand(command):
    """
    Runs a command and returns [returned value, execution time].
    """
    start_time = time.time()

    returncode = syscmd(command)

    end_time = time.time()

    return [returncode, round(end_time - start_time, 4)]