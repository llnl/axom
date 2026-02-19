#!/bin/bash

# Copyright (c) Lawrence Livermore National Security, LLC and other
# Axom Project Contributors. See top-level LICENSE and COPYRIGHT
# files for dates and other details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#UNCRUSTIFY_EXEC=$DIR/../src/TPL/uncrustify/src/uncrustify

for file in $(find . -type f -regex ".*\.\(hpp\|cpp\)" -follow -print0 | xargs -0); do
#    $UNCRUSTIFY_EXEC -c uncrustify.cfg --no-backup $file
     $1 -c $2 --no-backup $file
    echo $file
done

    
#echo $FILES
