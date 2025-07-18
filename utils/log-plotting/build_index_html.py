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

# This script buildes the index html page which overviews all figures from icon log file data.
import os
import shutil
import textwrap
from argparse import ArgumentParser
from pathlib import Path


def main():
    options = parse_args()
    build(options.PATH, options.exp_id, options.prio_plots, options.plot_format)


def parse_args():
    """Parses the command line arguments"""
    global options
    parser = ArgumentParser()
    parser.description = "displayes timer plots in html page"
    parser.add_argument("PATH", help="name of directory containing all plots")
    parser.add_argument(
        "--exp_id", help="give name of experiment", required=True
    )
    parser.add_argument("--plot_format", default="png")
    parser.add_argument(
        "--prio_plots",
        nargs="*",
        help="give list of plots which should appear at the top",
        default=[],
    )
    options = parser.parse_args()
    return options


def div(txt, class_):
    """Format text as html div."""
    return f"<div class='{class_}'>{txt}</div>"


def title_parser(path):
    """Generate title out of keyword in path."""

    if path.split(".")[0] == "timer_report":
        title = f'Time of compute process {path.split(".")[-2]}'
    if path.split(".")[0] == "wind_speeds":
        discription = path.split(".")[-2]
        if discription == "max_VN+W_log":
            title = "Maximum wind speeds"
        elif discription == "max_V":
            title = "Maximum horizontal wind speed"
        elif discription == "max_W":
            title = "Maximum vertical wind speed"
        elif discription == "max_VN+W_level":
            title = "Levels of maximum horizontal and vertical wind speed"
        elif discription == "max_W+levelW":
            title = "Maximum of vertical wind speed and level"
        elif discription == "levelW_vs_W":
            title = "Distribution of maximum vertical wind speeds and levels of occurence"
        elif discription == "levelVN_vs_VN":
            title = "Distribution of maximum horizontal wind speeds and levels of occurence"
        else:
            title = ""
    if path.split(".")[0] == "sdpd":
        title = "Simulated days per day"
    return title


def img_id(path):
    """Get image ID."""
    return path.split(".")[-3] + path.split(".")[-2]


def img(path):
    """Generate the lines for one image card."""
    image_id = img_id(path)
    title = title_parser(path)
    entry_txt = div(
        div(
            "\n".join(
                [
                    f"<img id='{image_id}' src='{path}' "
                    "class='card-img-top'/>",
                    div(
                        "\n".join(
                            [
                                f'<h5 class="card-title">{title}</h5>',
                            ]
                        ),
                        "card-body",
                    ),
                ]
            ),
            "card",
        ),
        "col",
    )
    return entry_txt


def prioritized_plots(prio_plots, lines, avail, plot_format):
    """Add lines of prioritized plots."""
    lines_prio_plots = []
    count = 0
    if prio_plots is not None:
        for f in prio_plots:
            if f in avail:
                lines_prio_plots.append(img(f + "." + plot_format))
                count += 1
    if count != 0:
        lines.append("<h4>Priority performance plots</h4>")
        lines.append('<div class="row row-cols-1 row-cols-md-3 g-4">')
        lines.extend(lines_prio_plots)
        lines.append("</div>")


def wind_speed(lines, plots_exp, plots_job, plot_format):
    """Generate the lines for one wind speed image cards."""
    plot_tuples = zip(plots_exp, plots_job)
    for plot_exp, plot_job in plot_tuples:
        plot_exp = plot_exp + "." + plot_format
        plot_job = plot_job + "." + plot_format
        image_id_exp = img_id(plot_exp)
        img_text = div(
            div(
                "\n".join(
                    [
                        f"<img id='{image_id_exp}' src='{plot_exp}' "
                        "class='card-img-top'/>",
                        "\n".join(
                            [
                                f"<img src='{plot_job}' "
                                "class='card-img-bottom'/>",
                            ]
                        ),
                        div(
                            "\n".join(
                                [
                                    f'<h5 class="card-title">{title_parser(plot_exp)}</h5>',
                                ]
                            ),
                            "card-body",
                        ),
                    ]
                ),
                "card",
            ),
            "col",
        )
        lines.append('<div class="card-group">')
        lines.append(img_text)
        lines.append("</div>")


def timer_report(lines, timer, plot_format):
    """Add lines of timer report plots."""
    for f in timer:
        lines.append(img(f + "." + plot_format))


def index(lines, plots_exp, timer, plot_format):
    """Generate lines for the index of all figuers located at the bottom page."""
    all_plots = plots_exp + timer
    lines.append('<div class="sticky-bottom"><nav class="nav flex-column">')
    for path in all_plots:
        path = path + "." + plot_format
        # image_id = path.split(".")[-3] + path.split(".")[-2]
        image_id = img_id(path)
        title = title_parser(path)
        lines.append(
            f'<a class="nav-link" href="#{image_id}">{title}</a><br />'
        )
    lines.append("</nav></div>")


def main_text(path, exp_id, prio_plots, plot_format):
    """Generate the lines of text for the overview page."""
    plots_job = []
    plots_exp = []
    timer = []
    lines = []

    for plot in sorted(os.listdir(path)):
        if plot.split(".")[0] == "wind_speeds":
            if plot.split(".")[2] == "last_job":
                plots_job.append(plot.rsplit(".", 1)[0])
                # plots_job.append(plot)
            else:
                plots_exp.append(plot.rsplit(".", 1)[0])
                # plots_exp.append(plot)
        elif plot.split(".")[-1] != plot_format:
            None
        else:
            timer.append(plot.rsplit(".", 1)[0])

    lines.append(f"<h3>Experiment {exp_id}</h3>")
    prioritized_plots(
        prio_plots, lines, plots_job + plots_exp + timer, plot_format
    )
    lines.append("<h4>Wind speed plots</h4>")
    lines.append('<div class="row row-cols-1 row-cols-md-3 g-4">')
    wind_speed(lines, plots_exp, plots_job, plot_format)
    lines.append("</div>")
    lines.append("<h4>All performance plots</h4>")
    lines.append('<div class="row row-cols-1 row-cols-md-3 g-4">')
    timer_report(lines, timer, plot_format)
    lines.append("</div>")
    index(lines, plots_exp, timer, plot_format)
    return lines


def build(path, exp_id, prio_plots, plot_format):
    """Write lines to index.html."""
    shutil.copyfile(
        Path(__file__).parent / "bootstrap.min.css",
        Path(path) / "bootstrap.min.css",
    )

    header = textwrap.dedent(
        """
    <!doctype html>
    <html lang="en">
      <head>
        <!-- Required meta tags -->
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">

        <!-- Bootstrap CSS -->
        <link href="bootstrap.min.css" rel="stylesheet">
        <title>Log Monitoring</title>
      </head>
      <body>
        <div class="container-fluid">
          <h1>
            Log Monitoring
          </h1>
          <p>
            This page provides an overview of figures plotted with data which was extracted from the icon log file.
          <p>
          <input class="form-control searchbox-input" type="text" placeholder="Type something here to search...">
          <br>
    """
    )

    footer = textwrap.dedent(
        """
        </div>
        <script>
          $(document).ready(function(){
            $('.searchbox-input').on("keyup", function() {
              var value = $(this).val().toLowerCase();
              $(".col").filter(function() {
                $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
              });
            });
          });
        </script>
      </body>
    </html>
    """
    )

    with open(f"{path}/index.html", "w") as html_file:
        text = (
            header
            + "\n".join(main_text(path, exp_id, prio_plots, plot_format))
            + footer
        )
        html_file.write(text)


if __name__ == "__main__":
    main()
