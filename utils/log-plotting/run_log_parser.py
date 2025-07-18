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

# This script takes an icon (atmosphere) log file and extracts the timer report information.

import importlib
import logging
import os
import re
import sys
from argparse import ArgumentParser


class RunLogParser:
    """Extract data from ICON log file."""

    def __init__(self, output_dir, job_id, exp_id, custom_modules):
        self.output_dir = output_dir
        self.job_id = job_id
        self.exp_id = exp_id
        self.custom_modules = custom_modules
        self.analyzers = []
        self.get_custom_analyzers()

    def process_file(self, f):
        self.parse_file(f)
        self.prepare_output_directory()
        for a in self.analyzers:
            a.save(f, self.output_dir)

    def prepare_output_directory(self):
        if self.output_dir != ".":
            os.makedirs(self.output_dir, exist_ok=True)

    def parse_file(self, infile):
        """Read log file and extract data."""
        logging.debug("parsing %s" % infile)
        with open(infile, "r") as file:
            for line in file:
                line = self.strip_line(line)
                for a in self.analyzers:
                    a.analyze_line(line)
        file.close
        for a in self.analyzers:
            a.processing(self.exp_id, self.job_id)

    def get_custom_analyzers(self):
        """Import custom analyzer scripts."""
        for script in self.custom_modules:
            if (spec := importlib.util.find_spec(script)) is not None:
                module = importlib.util.module_from_spec(spec)
                sys.modules[script] = module
                spec.loader.exec_module(module)
                logging.info(f"{script!r} has been imported")
                if module.custom_analyzers is not None:
                    for obj in module.custom_analyzers:
                        self.analyzers.append(obj())
                        logging.info(f"{obj!r} has been added to analyzers")
                else:
                    logging.error(
                        f"Can't find the analyzer(s) in module {script!r}"
                    )
                    sys.exit(1)
            else:
                logging.error(f"Can't find the module {script!r}")
                sys.exit(1)
        if not self.custom_modules:
            logging.info("No custom analyzer specified")

    def strip_line(self, line):
        line = line.strip()
        if mo := re.match(
            r"([-0-9:.]+T[-0-9:.]+)", line
        ):  # T to make sure this is a timer, not a rank.
            line = line[mo.end() :]
        return line


options = {}


def parse_args():
    """Parses the command line arguments"""
    parser = ArgumentParser()
    parser.description = "extract timer report from icon log file"
    parser.add_argument(
        "FILES", nargs="+", help="names of (multiple) log file(s)"
    )
    parser.add_argument(
        "-o", "--output_dir", help="give path to output directory"
    )
    parser.add_argument(
        "--job_id", help="give job id to be processed", required=True
    )
    parser.add_argument(
        "--exp_id", help="give name of experiment", required=True
    )
    parser.add_argument(
        "--custom_modules",
        nargs="*",
        help="python script to analyze the log customly",
    )
    options = parser.parse_args()
    if options.output_dir is None:
        logging.warning(
            """
        Output directory not given.
        File will be saved in the current directory.
        Path can be supplied with --output_dir.
        """
        )
        options.output_dir = "."
    return options


if __name__ == "__main__":
    options = parse_args()

    rlp = RunLogParser(
        output_dir=options.output_dir,
        job_id=options.job_id,
        exp_id=options.exp_id,
        custom_modules=options.custom_modules,
    )
    for f in options.FILES:
        rlp.process_file(f)
