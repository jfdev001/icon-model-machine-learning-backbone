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

import logging
import re

import pandas as pd
from BaseAnalyzer import BaseAnalyzer_


class AnalyzeTimer(BaseAnalyzer_):
    """Analyze lines of timer reports in log file"""

    def __init__(self):
        super().__init__()
        self.data = []
        self.ranks = {}
        self.filename = "timer_report"

    def timer_ranks(self, line):
        """Extract the timer report ranks."""
        rank = line.split(":")[0]
        logging.debug(f"{rank} for {line}")
        if line.endswith("Timer report"):
            self.ranks[rank] = rank
        if " Timer report, ranks" in line:
            self.ranks[rank] = line.split(" ")[-1]

    def fill_single_match(self, match, rank):
        if match:
            groups = match.groups()
            filled = [
                *groups[0:5],
                rank,
                *groups[5:7],
                rank,
                groups[7],
                rank,
                groups[8],
                rank,
                groups[9],
            ]
            return filled

    def analyze_line(self, line):
        self.timer_ranks(line)
        if len(self.ranks) == 0:
            return
        number = "[0-9.E\+hms]+"
        columns = f"\s*({number})\s*({number})\s*\[({number})\]\s*({number})\s*({number})\s*\[({number})\]\s*({number})\s*\[({number})\]\s*({number})\s*\[({number})\]\s*({number})\s*{number}$"
        columns_single = f"\s*({number})\s*({number})\s*({number})\s*({number})\s*({number})\s*({number})\s*({number})$"

        rank = line.split(":")[0]

        process = re.match(rf"([0-9\s]+):  ()(.+?){columns}", line)
        if process:
            process = process.groups()
        else:
            process = self.fill_single_match(
                re.match(rf"([0-9\s]+):  ()(.+?){columns_single}", line), rank
            )

        sub_process = re.match(rf"([0-9\s]+):  (\s*)L\s(.+?){columns}", line)
        if not sub_process:
            sub_process = self.fill_single_match(
                re.match(rf"([0-9\s]+):  (\s*)L\s(.+?){columns_single}", line),
                rank,
            )
        else:
            sub_process = sub_process.groups()

        special_sub_process = re.match(r"([0-9\s]+):  (\s*)L\s(.+±)", line)
        if sub_process:
            self.process_pattern(sub_process)
        elif process:
            self.process_pattern(process)
        elif special_sub_process:
            offset = len(special_sub_process.groups()[1])
            self.data.append(
                (
                    special_sub_process.groups()[0],
                    offset,
                )
                + special_sub_process.groups()[2:]
                + ("",) * 11
            )

        return self.data

    def processing(self, exp_id, job_id):
        column_names = (
            "rank",
            "offset",
            "name",
            "#calls",
            "t_min(s)",
            "min_rank",
            "t_avg(s)",
            "t_max(s)",
            "max_rank",
            "total_min(s)",
            "total_min_rank",
            "total_max(s)",
            "total_max_rank",
            "total_avg(s)",
        )

        self.df = pd.DataFrame(self.data, columns=column_names)

        self.df["ranks"] = [self.ranks[x] for x in self.df["rank"]]
        offset = self.df["offset"].values.tolist()
        names = self.df["name"].values.tolist()
        names = self.name(offset, names)
        self.df.name = names
        self.df = self.df.drop("offset", axis=1)
        self.df["exp_name"] = [exp_id] * len(offset)
        self.df["job_id"] = [job_id] * len(offset)
        return self.df

    def time(self, time):
        """Reformat time in units of seconds."""
        t_hm = re.match(r"([0-9]+)h([0-9]+)m", time)
        t_ms = re.match(r"([0-9]+)m([0-9]+)s", time)
        t_s = re.match(r"([0-9.]+)s", time)
        if t_hm:
            t_hm = t_hm.groups()
            t_s = int(t_hm[0]) * 3600 + int(t_hm[1]) * 60
        elif t_ms:
            t_ms = t_ms.groups()
            t_s = int(t_ms[0]) * 60 + int(t_ms[0])
        elif t_s:
            t_s = float(t_s.group(1))
        return t_s

    def process_pattern(self, line):
        line = list(line)
        for i, elm in enumerate(line):
            if i == 1:
                line[i] = len(elm)
            elif i == 3:
                line[i] = float(elm)
            elif i in [4, 6, 7]:
                line[i] = self.time(elm)
        return self.data.append(line)

    def name(self, offset, names):
        """Rewrite names of subprocesses as hirarchie 'directory'."""
        layer = [int((i - 1) / 3 + 1) if i > 0 else 0 for i in offset]
        directories = []
        for j in range(len(names)):
            directory = []
            for i in range(j + 1):
                if i == 0:
                    directory.append(names[i])
                elif layer[i] > layer[i - 1]:
                    directory.append(names[i])
                elif layer[i] <= layer[i - 1]:
                    directory = directory[: layer[i]]
                    directory.append(names[i])
            separator = "/"
            directories.append(separator.join(directory))
        return directories


custom_analyzers = [
    AnalyzeTimer,
]
