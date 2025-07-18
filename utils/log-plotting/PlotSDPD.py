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

# This script plots the simulated day per day (SDPD) of one experiment.

import matplotlib.pyplot as plt
import numpy as np
from BasePlotter import BasePlotter_, get_parser, parse_args


def main():
    parser = get_parser()
    parser.description = (
        "Plots simulated days per day (SDPD) of one experiment."
    )

    options = parse_args(parser)
    print("options:", options)

    pltSDPD = PlotSDPD(
        options["csv_dir"], options["output_dir"], options["plot_format"]
    )
    pltSDPD.plot()


class PlotSDPD(BasePlotter_):
    def __init__(self, csv_dir, output_dir, plot_format, figsize=[10, 6]):
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
        csv_name = "exp_status"
        super().__init__(csv_dir, csv_name, output_dir, plot_format)
        self.figsize = figsize

    def plot(self):
        """Plots simulated days per day (SDPD) of one experiment."""
        df = self.extract_data()
        x = df["run_start"]
        y = df["sdpd"]

        sdpd_mean = np.mean(y)
        sdpd_std = np.std(y)

        fig = plt.figure(figsize=self.figsize)
        plt.plot(x, y, ".-")

        plt.ylabel("Throughput [SPDP]")
        plt.xlabel("start time of the run")
        plt.xticks(x, [date[:10] for date in x], rotation=80, ha="right")
        plt.margins(x=0)
        plt.subplots_adjust(bottom=0.27)

        unique_df = df["exp_id"].drop_duplicates()
        if unique_df.size == 1:
            exp_id = unique_df.to_list()[0]
            plt.title(f"Throughput of {exp_id}")
            plt.figtext(
                0.65,
                0.3,
                f"Mean: ({round(sdpd_mean, 2)} +/- {round(sdpd_std, 2)}) SDPD",
            )
            plt.savefig(f"{self.output_dir}/sdpd.{exp_id}.{self.plot_format}")
        else:
            exp_ids = unique_df.to_list()
            plt.title(f"Throughput of {exp_ids[0]} ff.")
            plt.figtext(
                0.65,
                0.35,
                f"Experiments: {exp_ids}\n"
                + f"Mean: ({round(sdpd_mean, 2)} +/- {round(sdpd_std, 2)}) SDPD",
            )
            plt.savefig(
                f"{self.output_dir}/{exp_ids[0]}_ff_sdpd_plot.{self.plot_format}"
            )
        plt.close(fig)


if __name__ == "__main__":
    main()
