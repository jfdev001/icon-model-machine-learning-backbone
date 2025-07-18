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

import argparse
import os
import sys

needle = "Memory Size Used (MB) Min"


def error(message, exit=1):
    sys.stdout.flush()
    sys.stderr.write(os.path.basename(sys.argv[0]) + ": " + message + "\n")
    sys.stderr.flush()
    sys.exit(exit)


def scan_stdin():
    mem = -1
    for line in sys.stdin:
        sys.stdout.write(line)  # forward line to Job log

        # Look for
        # Memory Size Used (MB) Min               :    48517.996
        if line.startswith(needle):
            first, second = line.split(":")
            second = second.strip().split(" ")[0]
            mem = max(mem, float(second))

    if mem == -1:
        error("Could not find 'Memory Size Used'", exit=2)

    return mem


message_template = """
###############################################################################
{error}
###############################################################################

Possible reasons for increased memory consumption
=================================================
* You enabled additional output variables.
* You allocated additional variables.
* You reduced the total number of VEs used.

Possible reasons for unexpected low consumption
===============================================
* You disabled additional output variables.
* You removed variables.
* You increased the total number of VEs used.

General reasons
===============
* You changed the hardware platform (e.g. SX30 instead of SX10)
* Update of external dependency
* Bad luck with the sampling of memory usage. In this case, please try again or
  consider increasing the limit for this experiment slightly.

The limits are defined in the experiment under `run/checksuite.nwp/nwpexp...`
Please keep the limits tight but practical. The goal of this test is to detect
additional 3D-variables which use about 180 to 220 MB.

###############################################################################
{error}
###############################################################################
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check the memory usage reported by NEC mpi"
    )
    parser.add_argument(
        "--upper_mem_limit",
        default=0,
        help="Upper limit of the " + needle,
        type=int,
    )
    parser.add_argument(
        "--lower_mem_limit",
        default=0,
        help="Lower limit of the " + needle,
        type=int,
    )

    args = parser.parse_args()

    mem = scan_stdin()

    if args.upper_mem_limit != 0:
        if mem > args.upper_mem_limit:
            error(
                message_template.format(
                    error="'%s' is higher than expected. (actual: %g; limit: %g)"
                    % (needle, mem, args.upper_mem_limit)
                ),
                exit=21,
            )
    if args.lower_mem_limit != 0:
        if mem < args.lower_mem_limit:
            error(
                message_template.format(
                    error="'%s' is lower than expected. (actual: %g; limit: %g)"
                    % (needle, mem, args.lower_mem_limit)
                ),
                exit=21,
            )
