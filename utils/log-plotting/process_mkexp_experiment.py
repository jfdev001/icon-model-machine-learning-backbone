#! %{JOB.python3} #%# -*- mode: python -*- vi: set ft=python :

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

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


import fileinput
import logging
import os
import sys

import build_index_html
import exp_status_log_parser
import run_log_parser
import run_log_plotter

try:
    import dateutil

    time_guesser = dateutil.parser.parse
except ModuleNotFoundError:
    import mtime

    time_guesser = mtime.DateTime


def main(args):
    plot_format = args.plot_format
    custom_analyzer = args.custom_analyzer
    custom_plotter = args.custom_plotter
    exp_id = args.exp_id
    if args.job_id is not None:
        job_id = args.job_id
    else:
        job_id = get_job_id(args.script_dir, exp_id, args.start_date)
    log_prefix = f"{exp_id}.run.{job_id:0>8}"
    log_file = get_log_file_name(args.log_dir, log_prefix)
    csv_dir = f"{args.mon_dir}/timer_data"
    output_dir_plots = f"{args.mon_dir}/log_plotting/"

    apply_run_log_parser(
        log_prefix, log_file, csv_dir, job_id, exp_id, custom_analyzer
    )
    exp_status_log = f"{args.script_dir}/{exp_id}.log"
    apply_exp_status_log_parser(exp_status_log, csv_dir, exp_id)
    apply_run_log_plotter(
        csv_dir, output_dir_plots, plot_format, custom_plotter
    )
    run_sdpd_plotter(csv_dir, output_dir_plots, plot_format)
    apply_build_index_html(
        output_dir_plots, exp_id, args.prioritized_plots, plot_format
    )


def get_job_id(script_dir, exp_id, start_date):
    current_job_id = None
    log_file = f"{script_dir}/{exp_id}.log"
    for line in fileinput.input(log_file):
        logging.debug(line)
        (timestamp, startdate, enddate, job_id, state) = line.rstrip().split(
            " "
        )
        logging.debug(str((timestamp, startdate, enddate, job_id, state)))
        if state == "end" and time_guesser(startdate) == time_guesser(
            start_date
        ):
            current_job_id = job_id
    if current_job_id is None:
        logging.critical(
            f"Could not find a finished job for start_date {start_date} in {log_file}."
        )
        sys.exit(1)
    logging.debug(f"current job id: {current_job_id}")
    return current_job_id


def get_log_file_name(log_dir, log_prefix):
    log_file = f"{log_dir}/{log_prefix}.log"
    logging.debug(f"log file: {log_file}")
    return log_file


def apply_run_log_parser(
    log_prefix, log_file, csv_dir, job_id, exp_id, custom_analyzer
):
    timer_csv_file = f"{csv_dir}/{log_prefix}.timer_report.csv"
    wind_csv_file = f"{csv_dir}/{log_prefix}.wind_speed.csv"
    if not (os.path.exists(timer_csv_file) and os.path.exists(wind_csv_file)):
        rlp = run_log_parser.RunLogParser(
            output_dir=csv_dir,
            job_id=job_id,
            exp_id=exp_id,
            custom_modules=custom_analyzer,
        )
        rlp.process_file(log_file)


def apply_exp_status_log_parser(exp_status_log, csv_dir, exp_id):
    exp_status_log_parser.process_file(exp_status_log, csv_dir, exp_id)


def apply_run_log_plotter(
    csv_dir, output_dir_plots, plot_format, custom_plotter
):
    rlp = run_log_plotter.RunLogPlotter(
        csv_dir=csv_dir,
        output_dir=output_dir_plots,
        plot_format=plot_format,
        custom_plotter=custom_plotter,
    )
    rlp.plot()


def run_sdpd_plotter(csv_dir, output_dir_plots, plot_format):
    from PlotSDPD import PlotSDPD

    pltSDPD = PlotSDPD(
        csv_dir=csv_dir, output_dir=output_dir_plots, plot_format=plot_format
    )
    pltSDPD.plot()


def apply_build_index_html(
    output_dir_plots, exp_id, prioritized_plots, plot_format
):
    os.makedirs(output_dir_plots, exist_ok=True)
    build_index_html.build(
        output_dir_plots, exp_id, prioritized_plots, plot_format=plot_format
    )
