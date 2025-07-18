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

#'This script plots timers of icon logfiles.'

import logging

import matplotlib.pyplot as plt
import numpy as np

try:
    import seaborn as sns

    sns.set()
except ModuleNotFoundError:
    pass

from BasePlotter import BasePlotter_, get_parser, parse_args


def main():
    parser = get_parser()
    parser.description = "plots timer report data from icon log file"
    parser.add_argument(
        "-p", "--process", help="name of the compute process", default="all"
    )
    parser.add_argument(
        "-i",
        "--individual",
        help="plotting of timers for individual calls, otherwise plot total timers",
        action="store_true",
    )
    parser.add_argument("--figsize", nargs=2, type=float, default=[8, 6])

    options = parse_args(parser)
    logging.info("options:", options)

    pltt = PlotTimer(
        csv_dir=options["csv_dir"],
        output_dir=options["output_dir"],
        plot_format=options["plot_format"],
        figsize=options["figsize"],
    )
    pltt.plot(
        individual=options.get("individual", False), process=options["process"]
    )


class PlotTimer(BasePlotter_):
    """Plot timer of the icon log file.

    Parameters
    ----------
    files : str
        The directory of csv files containing the timer data. These csv files can be produced by AnalyzeTimer.
    output_dir : str
        The directory where the plots are going to be stored.
    plot_format : str
        The image file format of the generated plots.
    figsize : touple, optional
        The figure size of the plots.

    Methods
    -------

    """

    def __init__(self, csv_dir, output_dir, plot_format, figsize=[8, 6]):
        """Parameters
        ----------
        files : str
            The directory of csv files containing the timer data. These csv files can be produced by AnalyzeTimer.
        output_dir : str
            The directory where the plots are going to be stored.
        plot_format : str
            The image file format of the generated plots.
        figsize : touple, optional
            The figure size of the plots."""
        csv_name = "timer_report"
        super().__init__(csv_dir, csv_name, output_dir, plot_format)
        self.figsize = figsize

    def plot_single_process(self, individual, process):
        """Plots a single figure of a compute process.

        This plot can either contain the data of a single call or averaged the over total amount of calls.

        Parameters
        ----------
        individual : bool
            Indicates single call (True) or total (False).
        process : str
            The name of the compute process.
        """

        df_all = self.extract_data()
        df = df_all[df_all["name"] == process].copy()
        df["x"] = self.xaxis(df["job_id"].values.tolist())

        self.color_palette(df)

        fig, ax = plt.subplots()
        if individual:
            df["yerr_min"] = df["t_avg(s)"] - df["t_min(s)"]
            df["yerr_max"] = df["t_max(s)"] - df["t_avg(s)"]
            self.scatter_plot("t_avg(s)", df, ax)
            specification = "individual"
        else:
            df["yerr_min"] = df["total_avg(s)"] - df["total_min(s)"]
            df["yerr_max"] = df["total_max(s)"] - df["total_avg(s)"]
            self.scatter_plot("total_avg(s)", df, ax)
            specification = "total"

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())
        plt.ylabel(specification + " time in s")
        plt.xlabel("job_id")
        plt.xticks(
            df.x,
            [
                (
                    str(df.job_id.values[i])
                    if df.job_id.values[i] != df.job_id.values[i - 1]
                    else None
                )
                for i in range(len(df))
            ],
            rotation=70,
        )
        plt.subplots_adjust(bottom=0.27)
        plt.title("time of compute process " + process, wrap=True)

        replace = str.maketrans("/.", "-_")
        plt.savefig(
            f"{self.output_dir}/timer_report.{specification}."
            + process.translate(replace)
            + "."
            + self.plot_format
        )
        plt.close(fig)

    def color_palette(self, df):
        """Generates a color palette depending on the ranks."""

        translations = {b: a for a, b in enumerate(sorted(set(df.ranks)))}
        try:
            pal = sns.color_palette("husl", len(translations))
        except NameError as e:
            logging.info(f"Seaborn color palette creation failed: {e}")
            pal = [plt.cm.PiYG(i) for i in np.linspace(0, 1, len(translations))]
        self.color = [pal[translations[x]] for x in df.ranks]
        return self.color

    def scatter_plot(self, t_avg, df, ax):
        """Makes a scatter plot."""

        [
            df.iloc[[i]].plot.scatter(
                "x",
                t_avg,
                ax=ax,
                s=50,
                label=label,
                zorder=1,
                color=self.color[i % len(self.color)],
            )
            for i, label in enumerate(df.ranks)
        ]
        try:
            ax.errorbar(
                "x",
                t_avg,
                yerr=(df.yerr_min, df.yerr_max),
                data=df,
                fmt="none",
                ecolor=self.color,
                zorder=2,
                label=None,
            )
        except ValueError as e:
            logging.critical(f"Can't process \n{str(df)}")
            raise e

    def xaxis(self, job_id):
        """Generates a list out of job-IDs to use as an X-axis."""

        x = [0]
        for i in range(len(job_id) - 1):
            if job_id[i] == job_id[i + 1]:
                x.append(x[i] + 1)
            else:
                x.append(x[i] + 2)
        return x

    def plot(self, individual=False, process="all"):
        """Plots the time of a compute process against its job-ID.

        This is done either for one process or all processes."""
        df_all = self.extract_data()
        if process == "all":
            logging.info("plotting all processes")
            names = sorted(set(df_all["name"].values.tolist()))
            for process in names:
                self.plot_single_process(individual, process)
        else:
            logging.info("plotting", process)
            self.plot_single_process(individual, process)


custom_plotter = [
    PlotTimer,
]

if __name__ == "__main__":
    main()
