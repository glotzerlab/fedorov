# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.


def json_key_to_int(x):
    if isinstance(x, dict):
        return {int(k): v for k, v in x.items()}
    return x
