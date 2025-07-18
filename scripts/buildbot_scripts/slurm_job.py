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

import re
import subprocess
import sys

from batch_job import BatchJob


class SlurmJob(BatchJob):
    def __init__(self, cmd, cwd):
        super().__init__(cmd, cwd)
        self.system = "Slurm"

    def submit(self, script):
        # --wait so that the subprocess waits for slurm completion
        # hence there no extra poll method needed like for PBS. the associated process can be used
        submit_cmd = self.cmd.split() + ["--wait"]

        if len(self.parents) > 0:
            parent_ids = [p.jobid for p in self.parents]
            submit_cmd.append(
                "--dependency=afterany:{}".format(",".join(parent_ids))
            )

        submit_cmd.append(script)
        print("submitting slurm job: '{}'".format(" ".join(submit_cmd)))
        sp = subprocess.Popen(
            submit_cmd,
            shell=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=self.cwd,
            encoding="UTF-8",
        )
        try:
            self.jobid = re.findall(r"\d+", sp.stdout.readline())[0]
        except:
            for line in sp.stderr.readlines():
                print(line)

        if not self.jobid or not self.jobid.isnumeric():
            print(
                "Parsing jobid from slurm job failed, got {}".format(self.jobid)
            )
            sys.exit(1)
        self.job = sp

    def poll(self, timeout):
        """Check if task is still running.

        Waits up to specified timeout in seconds for job to finish. If job
        finishes in time or has finished before, return True and set returncode.
        """
        try:
            returncode = self.job.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            return False
        else:
            self.returncode = returncode
            return True

    def cancel(self):
        if None is not self.jobid:
            qdel = subprocess.Popen(
                ["scancel", "-v", self.jobid],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.cwd,
                encoding="UTF-8",
            )
            print(qdel.stderr.readlines()[0].rstrip("\r\n"))
        else:
            print("Cannot find jobid to cancel job!")

    # check if a job was canceled
    def wasCanceled(self):
        if None is not self.jobid:
            checkState = subprocess.Popen(
                f"sacct -j{self.jobid} -Pn",
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=self.cwd,
                encoding="UTF-8",
            )
            out, err = checkState.communicate(timeout=2)
            jobState = out.split("|")[-2].split(" ")[0]
            return "CANCELLED" == jobState
        else:
            print("Cannot find jobid!")
            return False
