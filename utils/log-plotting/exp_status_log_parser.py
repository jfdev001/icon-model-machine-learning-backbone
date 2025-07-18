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

# This script extracts information from an icon cycle time log.

import datetime
import logging
from argparse import ArgumentParser

import pandas as pd


def main():
    options = parse_args()
    process_file(options.FILE, options.output_dir, options.exp_id)


options = {}


def parse_args():
    """Parses the command line arguments"""
    parser = ArgumentParser()
    parser.description = "extract cycle times from experiment log file"
    parser.add_argument("FILE", help="name of exp_id.log file")
    parser.add_argument(
        "-o",
        "--output_dir",
        help="give path to output directory",
    )
    parser.add_argument("--exp_id", help="give name of experiment")
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
    if options.exp_id is None:
        logging.critical(
            "Could not determine experiment name, please supply with --exp_id"
        )
    return options


def process_file(infile, output_dir, exp_id):
    """read file and get run times of one cycle"""
    logging.debug("parsing %s" % infile)
    data = []
    with open(infile, "r") as file:
        for line in file:
            data.append(line.rstrip().split(" "))
    file.close
    df = sort_data(data)
    df["exp_id"] = exp_id
    output_filename = output_dir + "/" + str(exp_id) + ".exp_status.csv"
    df.to_csv(output_filename)
    return df


def sort_data(data):
    data_sorted = []
    for i, elm in enumerate(data):
        if elm[-1] == "end":
            sdpd = compute_sdpd(elm[1], elm[2], data[i - 1][0], elm[0])
            data_sorted.append([sdpd] + [data[i - 1][0]] + elm[:-1])

    df = pd.DataFrame(
        data_sorted,
        columns=(
            "sdpd",
            "run_start",
            "run_end",
            "sim_start",
            "sim_end",
            "job_id",
        ),
    )
    return df


def compute_sdpd(sim_start, sim_end, run_start, run_end):
    reformat = datetime.datetime.strptime
    f1 = "%Y-%m-%dT%H:%M"
    f2 = "%Y-%m-%dT%H:%M:%SZ"
    sdpd = (reformat(sim_end, f1) - reformat(sim_start, f1)) / (
        reformat(run_end, f2) - reformat(run_start, f2)
    )
    return sdpd


if __name__ == "__main__":
    main()
