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
import subprocess
from abc import ABC, abstractmethod


class BatchJob(ABC):
    def __init__(self, cmd, cwd):
        self.system = "undefined"
        self.cmd = cmd
        self.cwd = cwd
        self.jobid = None
        self.job = None
        self.returncode = None
        self.parents = []

    def add_parent(self, parent):
        self.parents.append(parent)

    @abstractmethod
    def submit(self, script):
        pass

    @abstractmethod
    def poll(self, timeout):
        pass

    @abstractmethod
    def cancel(self):
        pass

    @abstractmethod
    def wasCanceled(self):
        pass

    def failed(self):
        # this only makes sense IF there is a returncode at all
        if None == self.returncode:
            print("This process has not yet returned any code!")
            return None

        _returncode = self.returncode
        _canceled = self.wasCanceled()
        return (0 != _returncode) or _canceled
