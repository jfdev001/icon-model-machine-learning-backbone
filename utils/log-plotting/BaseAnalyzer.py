#!/usr/bin/env python3

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import os


class BaseAnalyzer_:
    """
    General parameters and functionalities to analyze icon log files and extract the data.

    Methods
    -------
    cut_ranks(line)
        Cuts the ranks from the beginning of the line
    analyze_line(line)
        Extracts the data of intrest from the lines of the log-file
    processing(exp_id, job_id)
        Creates a pandas dataframe out of the collected data
    save(f, output_dir)
        Saves the dataframe to a csv file
    """

    def __init__(self):
        """Parameters
        ----------
        df : pandas dataframe
            The dataframe to which analyzed data is stored to.
        filename : str
            The name of the csv file the df is stored in."""

        self.df = None
        self.filename = None

    def cut_ranks(self, line):
        """Cuts the ranks from the beginning of the line"""
        lsplit = line.split(":")
        try:
            int(lsplit[0])
            return ":".join(lsplit[1:]).strip()
        except ValueError:
            return line

    def analyze_line(self, line):
        """Extracts the data of intrest from the lines of the log-file"""
        raise NotImplementedError

    def processing(self, exp_id, job_id):
        """Creates a pandas dataframe out of the collected data"""
        raise NotImplementedError

    def save(self, f, output_dir):
        """Saves the dataframe to a csv file"""
        fn_base = os.path.basename(f).split(".log")[0]
        os.makedirs(output_dir, exist_ok=True)
        output_filename = f"{output_dir}/{fn_base}.{self.filename}.csv"
        self.df.to_csv(output_filename)
        return output_filename
