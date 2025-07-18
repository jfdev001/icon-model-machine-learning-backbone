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
import fnmatch
import multiprocessing
import os
import re
import sys

FORTRAN_GLOB_PATTERNS = ["*.F90", "*.f90"]


def list_files(dirs_or_files):
    for dir_or_file in dirs_or_files:
        if not os.path.exists(dir_or_file):
            raise Exception(f"ERROR: '{dir_or_file}' does not exist")
        elif os.path.isfile(dir_or_file):
            yield dir_or_file
        elif os.path.isdir(dir_or_file):
            for subdir, _, filenames in os.walk(dir_or_file):
                for filename in filenames:
                    filepath = os.path.join(subdir, filename)
                    if any(
                        fnmatch.fnmatch(filepath, pattern)
                        for pattern in FORTRAN_GLOB_PATTERNS
                    ):
                        yield filepath
        else:
            print(
                f"WARNING: input argument '{dir_or_file}' "
                f"is neither a directory nor a file",
                file=sys.stderr,
            )


def check_file(filepath):
    with open(filepath, "rb") as f:
        raw = f.read(-1)
        txt = raw.decode("utf-8", errors="replace")

        # OpenMP sentinels (!$) are not allowed:
        omp_sentinel_match = re.search(r"(?m)^\s*!\$\s.*$", txt)

        # Fortran USE statements must not be interlined with the preprocessor
        # directives:
        interlined_use_match = re.search(
            r"(?im)^\n*(\s*use(?:\s*&\n|\s)(?:.*&\n)+\s*#.*)$", txt
        )

        return (
            filepath,
            omp_sentinel_match.group() if omp_sentinel_match else None,
            interlined_use_match.group(1) if interlined_use_match else None,
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Checks Fortran source files for discouraged patterns"
    )

    parser.add_argument(
        "files_or_dirs",
        nargs="+",
        metavar="FILE_OR_DIRECTORY",
        help="path to a file or directory to check",
    )

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    files_with_omp_sentinels = []
    files_with_interlined_uses = []
    with multiprocessing.Pool() as pool:
        for (
            filepath,
            omp_sentinel_string,
            interlined_use_string,
        ) in pool.imap_unordered(check_file, list_files(args.files_or_dirs)):
            if omp_sentinel_string:
                files_with_omp_sentinels.append((filepath, omp_sentinel_string))
            if interlined_use_string:
                files_with_interlined_uses.append(
                    (filepath, interlined_use_string)
                )

    exit_code = 0

    if files_with_omp_sentinels:
        exit_code = 1
        print(
            "ERROR: the following files contain OpenMP conditional compilation "
            "sentinels:\n"
            "\t{0}\n"
            "Replace the sentinels (!$) in the files above with the macro "
            "'#ifdef _OPENMP/#endif' directives.\n".format(
                "\n\t".join(
                    "{0}:\n{1}".format(n, s)
                    for n, s in sorted(files_with_omp_sentinels)
                )
            ),
            file=sys.stderr,
        )

    if files_with_interlined_uses:
        exit_code = 1
        print(
            "ERROR: the following files contain multi-line Fortran USE statements "
            "interlined with the preprocessor directives:\n"
            "\t{0}\n"
            "Avoid the interlining as described in "
            "https://gitlab.dkrz.de/icon/icon/-/merge_requests/471\n".format(
                "\n\t\n".join(
                    "{0}:\n{1}".format(n, s)
                    for n, s in sorted(files_with_interlined_uses)
                )
            ),
            file=sys.stderr,
        )

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
