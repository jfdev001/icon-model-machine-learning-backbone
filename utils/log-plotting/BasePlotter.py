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

import glob
import logging
import os
from argparse import ArgumentParser

import pandas as pd


def get_parser():
    """Create and return a new argument parser."""
    parser = ArgumentParser()

    parser.add_argument("csv_dir", help="name of csv file directory")
    parser.add_argument(
        "-o",
        "--output_dir",
        default=".",
        help="give path to output directory",
    )
    parser.add_argument("--plot_format", default="png")

    return parser


def parse_args(parser):
    """Parse command line arguments and return them as a dictionary."""
    op = parser.parse_args()
    if op.output_dir is None:
        logging.warning(
            """
        Output directory not given.
        File will be saved in the current directory.
        Path can be supplied with -d.
        """
        )
        op.output_dir = "."
    elif not os.path.exists(op.output_dir):
        logging.debug(f"Creating output directory {op.output_dir}")
    logging.info(f"Plots are saved to {op.output_dir}")

    options = vars(op)
    return options


class BasePlotter_:
    """
    General parameters and functionalities for icon log plotting.

    Parameters
    ----------
    files : str
        The directory of csv files containing the timer data. These csv files can be produced by AnalyzeTimer.
    output_dir : str
        The directory where the plots are going to be stored.
    plot_format : str
        The image file format of the generated plots.

    Methods
    -------
    extract_data()
        Extracts data from csv files into pandas dataframe.
    plot()
        Plots the extracted data.
    """

    def __init__(self, csv_dir, csv_name, output_dir, plot_format):
        """Parameters
        ----------
        files : str
            The directory of csv files containing the timer data. These csv files can be produced by AnalyzeTimer.
        output_dir : str
            The directory where the plots are going to be stored.
        plot_format : str
            The image file format of the generated plots."""

        self.files = sorted(glob.glob(f"{csv_dir}/*{csv_name}.csv"))
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.plot_format = plot_format

    def extract_data(self):
        """Extracts data from csv files into pandas dataframe."""
        dataframes = []
        if isinstance(self.files, list):
            for infile in self.files:
                df = pd.read_csv(infile)
                dataframes.append(df)
            results = pd.concat(dataframes)
        else:
            results = pd.read_csv(self.files)
        return results

    def plot(self):
        """Plots the extracted data."""
        raise NotImplementedError
