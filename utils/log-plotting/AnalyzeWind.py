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

import datetime
import logging
import re
import sys

import pandas as pd
from BaseAnalyzer import BaseAnalyzer_


class AnalyzeWind(BaseAnalyzer_):
    """Analyze lines in log file according to wind speed information."""

    def __init__(self):
        super().__init__()
        self.data = []
        self.dates = []
        self.ts = None
        self.filename = "wind_speed"

    def analyze_line(self, line):
        line = super().cut_ranks(line)
        if "Time step:" in line:
            self.ts = self.parse_ts(line)
            logging.debug(f"found {self.ts}")
        elif "MAXABS VN," in line:
            vn, lvn, w, lw = self.parse_vel(line)
            self.data.append((vn, lvn, w, lw))
            self.dates.append(self.ts)
        return self.data, self.dates

    def processing(self, exp_id, job_id):
        self.df = pd.DataFrame(
            self.data, columns=("vn", "level_vn", "w", "level_w")
        )
        self.df.insert(loc=0, column="dates", value=self.dates)
        self.df["exp_name"] = [exp_id] * len(self.df)
        self.df["job_id"] = [job_id] * len(self.df)
        if len(self.dates) == 0:
            logging.critical(
                f"Could not find time steps in log for {exp_id} {job_id}"
            )
            sys.exit(1)
        logging.info("Dates %s to %s" % (self.dates[0], self.dates[-1]))
        return self.df

    def parse_ts(self, line):
        """Extract datetime the last digits of a line"""
        return datetime.datetime.strptime(
            " ".join(line[:-4].split()[-2:]), "%Y-%m-%d %H:%M:%S"
        )

    def parse_vel(self, line):
        """the string is splitted into sections and the values are saved as floats/integers"""
        split = line.split()
        try:
            vn = float(split[6])
            lvn = int(split[9][:-1])
            w = float(split[10])
            lw = int(split[13][:-1])
        except ValueError:
            matches = re.findall("([-0-9.E+]+)", line)
            vn = matches[1]
            lvn = matches[3]
            w = matches[4]
            lw = matches[6]
        return vn, lvn, w, lw


custom_analyzers = [
    AnalyzeWind,
]
