# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

import axom.sidre as sidre
import numpy as np


def main():
    datastore = sidre.DataStore()
    root = datastore.getRoot()
    fields = root.createGroup("fields")

    values = np.arange(8, dtype=np.float64)
    fields.createView("values", values).apply(sidre.TypeID.DOUBLE_ID, len(values))

    view_data = fields.getView("values").getDataArray()
    assert np.array_equal(view_data, values)

    print(f"Using installed axom.sidre {sidre.__version__}")


if __name__ == "__main__":
    main()
