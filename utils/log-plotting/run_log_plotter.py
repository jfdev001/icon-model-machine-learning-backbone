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

import importlib
import logging
import sys

from BasePlotter import get_parser, parse_args


class RunLogPlotter:
    """Plot data from ICON log file."""

    def __init__(self, csv_dir, output_dir, plot_format, custom_plotter):
        self.csv_dir = csv_dir
        self.output_dir = output_dir
        self.plot_format = plot_format
        self.custom_plotter = custom_plotter
        self.plotters = []
        self.get_custom_plotter()

    def get_custom_plotter(self):
        """Import custom analyzer scripts."""
        for script in self.custom_plotter:
            if (spec := importlib.util.find_spec(script)) is not None:
                module = importlib.util.module_from_spec(spec)
                sys.modules[script] = module
                spec.loader.exec_module(module)
                logging.info(f"{script!r} has been imported")
                if module.custom_plotter is not None:
                    for obj in module.custom_plotter:
                        self.plotters.append(obj)
                        logging.info(f"{obj!r} has been added to plotters")
                else:
                    logging.error(
                        f"Can't find the plotter(s) in module {script!r}"
                    )
                    sys.exit(1)
            else:
                logging.error(f"Can't find the module {script!r}")
                sys.exit(1)
        if not self.custom_plotter:
            logging.info("No custom plotter specified")

    def plot(self):
        for P in self.plotters:
            p = P(
                csv_dir=self.csv_dir,
                output_dir=self.output_dir,
                plot_format=self.plot_format,
            )
            p.plot()


if __name__ == "__main__":
    parser = get_parser()
    parser.description = "Plot data from icon log file"
    parser.add_argument(
        "--custom_modules",
        nargs="*",
        help="python script to analyze the log customly",
    )
    options = parse_args(parser)

    rlp = RunLogPlotter(
        csv_dir=options["csv_dir"],
        output_dir=options["output_dir"],
        plot_format=options["plot_format"],
        custom_plotter=options["custom_modules"],
    )
    rlp.plot()
