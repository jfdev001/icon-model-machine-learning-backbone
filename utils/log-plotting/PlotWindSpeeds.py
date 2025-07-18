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

# This script takes an icon (atmosphere) log file and extracts the wind speed information.

from __future__ import print_function

import logging
import sys
from datetime import datetime

import dateutil
import matplotlib.pyplot as plt

logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)

try:
    import seaborn as sns

    sns.set()
except ModuleNotFoundError:
    pass

import matplotlib.dates as mdates

myFmt = mdates.DateFormatter("%Y-%m-%d %H:%M")

from BasePlotter import BasePlotter_, get_parser, parse_args

options = {}


def main():
    # parse_args()
    # logging.debug(str(options))
    parser = get_parser()
    parser.description = "Plot icon log file wind speed information"
    parser.add_argument(
        "-x", "--xlim", nargs=2, help="min and max time values, e.g. 2020-01-17"
    )
    options = parse_args(parser)
    if options.get("xlim", False):
        xl = options.get("xlim", False)
        print(dir(datetime))
        options["xlim"] = (
            dateutil.parser.parse(xl[0]),
            dateutil.parser.parse(xl[1]),
        )
    logging.info("options:", options)

    pltws = PlotWindSpeeds(
        options["csv_dir"], options["output_dir"], options["plot_format"]
    )
    pltws.plot()


class PlotWindSpeeds(BasePlotter_):
    def __init__(self, csv_dir, output_dir, plot_format):
        csv_name = "wind_speed"
        super().__init__(csv_dir, csv_name, output_dir, plot_format)

    def plot(self):
        """Plots extracted data depending on experiment/job ID."""

        data_exp = self.extract_data()

        # extract exp_id and plot all jobs
        unique_df = data_exp["exp_name"].drop_duplicates()
        if unique_df.size == 1:
            exp_id = unique_df.to_list()[0]
            self.plot_data(data_exp, exp_id)
        else:
            exp_ids = unique_df.to_list()
            logging.critical(f"More than one experiment ID detected. {exp_ids}")
            sys.exit(1)

        # extract last job_id and plot only last job_id
        current_job_id = data_exp["job_id"].max()
        if current_job_id is not None:
            data_job = data_exp.loc[
                data_exp["job_id"] == int(current_job_id), :
            ]
            self.plot_data(data_job, f"{exp_id}.last_job", current_job_id)
        else:
            logging.critical(f"Could not find a finished job in {self.files}.")
            sys.exit(1)

    def plot_data(self, data, filename_base, job_id=""):
        """Plots multiple graphs from extracted wind speed data."""

        self.plot_opts = {"linewidth": 0.05, "markersize": 1, "style": ".-"}
        if options.get("xlim", False):
            self.plot_opts["xlim"] = options.get("xlim", False)

        self.data = data
        self.filename_base = filename_base
        self.job_id = job_id

        self.plt_lineplot(
            title="Maximum wind speeds",
            ylabel="Wind speed in m/s",
            filename="max_VN+W_log",
            limit=(10, 400),
            y1="vn",
            color1="C0",
            label1="Max(VN)",
            y2="w",
            color2="C1",
            label2="Max(W)",
            logy=True,
            legend=True,
        )

        self.plt_lineplot(
            title="Maximum horizontal wind speed",
            ylabel="Wind speed in m/s",
            filename="max_V",
            limit=(100, 500),
            y1="vn",
            color1="C0",
            label1="Max(VN)",
        )

        self.plt_lineplot(
            title="Maximum vertical wind speed",
            ylabel="Wind speed in m/s",
            filename="max_W",
            limit=(10, 100),
            y1="w",
            color1="C1",
            label1="Max(W)",
        )

        self.plt_lineplot(
            title="Levels of maximum horizontal and vertical wind speed",
            ylabel="Level (counting down from the top)",
            filename="max_VN+W_level",
            limit=None,
            y1="level_vn",
            color1="C2",
            label1="Level of max VN",
            y2="level_w",
            color2="C3",
            label2="Level of max W",
            legend=True,
            inverty=True,
        )

        self.plt_lineplot(
            title="Maximum of vertical wind speed and level",
            ylabel="Wind speed in m/s, level",
            filename="max_W+levelW",
            limit=(0, 100),
            y1="w",
            color1="C1",
            label1="Max(W)",
            y2="level_w",
            color2="C3",
            label2="Level of max W",
        )

        self.plt_dist("w")
        self.plt_dist("vn")

    def plt_lineplot(
        self,
        title,
        ylabel,
        filename,
        limit,
        y1,
        color1,
        label1,
        y2=None,
        color2=None,
        label2=None,
        logy=False,
        legend=False,
        inverty=False,
    ):
        """Plots line plots with two axis"""

        plt.figure()
        self.data.plot(
            x="dates",
            y=y1,
            logy=logy,
            color=color1,
            label=label1,
            **self.plot_opts,
        )
        if y2 is not None and y2:
            self.data.plot(
                x="dates",
                y=y2,
                ax=plt.gca(),
                logy=logy,
                color=color2,
                label=label2,
                **self.plot_opts,
            )
        if inverty:
            ax = plt.gca()
            ax.invert_yaxis()

        plt.ylabel(ylabel)
        plt.title(f"{title}\n{self.job_id}", wrap=True)
        if legend:
            plt.legend(loc="upper right")
        plt.ylim(limit)
        plt.xticks(rotation=70, ha="right")
        plt.subplots_adjust(bottom=0.4)
        # plt.gca().xaxis.set_major_formatter(myFmt)
        plt.savefig(
            f"{self.output_dir}/wind_speeds.{self.filename_base}.{filename}.{self.plot_format}"
        )
        plt.close("all")

    def plt_dist(self, variable):
        """Plots wind speed distributions.

        variable : 'vn' or 'w'
            Indicates to either plot distribution on vertical ('w') or horizontal ('vn') wind speed.
        """

        plt.figure()
        scatter_opts = self.plot_opts.copy()
        scatter_opts["linewidth"] = 0
        scatter_opts["marker"] = "."
        del scatter_opts["style"]
        plt.plot(
            variable,
            f"level_{variable}",
            data=self.data,
            color="C4",
            label=f"Max({variable.upper()})",
            **scatter_opts,
        )
        plt.ylabel("Level (counting from the top)")
        plt.xlabel("Speed (m/s)")

        def title(direction, job_id):
            plt.title(
                f"Distribution of maximum {direction} wind speeds and levels of occurence\n{job_id}"
            )

        if variable == "w":
            plt.ylim((70, 0))
            plt.xlim((0, 100))
            title("vertical", self.job_id)
        elif variable == "vn":
            plt.ylim((50, 0))
            plt.xlim((100, 500))
            title("horizontal", self.job_id)
        plt.savefig(
            f"{self.output_dir}/wind_speeds.{self.filename_base}.level{variable.upper()}_vs_{variable.upper()}.{self.plot_format}"
        )


custom_plotter = [
    PlotWindSpeeds,
]

if __name__ == "__main__":
    main()
