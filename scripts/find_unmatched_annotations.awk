# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

# This utility script can help find unmatched pairs of AXOM_ANNOTATE_BEGIN/AXOM_ANNOTATE_END.
# To use, run the following from an axom source directory
#  > git grep AXOM_ANNOTATE | awk -f scripts/find_unmatched_annotations.awk


# Extract strings inside AXOM_ANNOTATE_BEGIN("") and AXOM_ANNOTATE_END("")
/AXOM_ANNOTATE_BEGIN\(/ {
    match($0, /AXOM_ANNOTATE_BEGIN\("([^"]+)"\)/, arr);
    if (arr[1] != "") {
        begin[arr[1]]++;
    }
}

/AXOM_ANNOTATE_END\(/ {
    match($0, /AXOM_ANNOTATE_END\("([^"]+)"\)/, arr);
    if (arr[1] != "") {
        end[arr[1]]++;
    }
}

# At the end of processing, print results
END {
    print "Matched Annotations and Counts:"
    for (key in begin) {
        print "  '" key "' (" begin[key] ")";
    }

    print "\nUnmatched Annotations:"
    unmatched = 0;
    for (key in begin) {
        if (begin[key] != end[key]) {
            print "  '" key "':\tBEGIN=", begin[key], "END=", (key in end ? end[key] : 0);
            unmatched++;
        }
    }
    for (key in end) {
        if (!(key in begin)) {
            print "  '" key "':'\tBEGIN=0", "END=", end[key];
            unmatched++;
        }
    }

    if (unmatched == 0) {
        print "All annotations have matching BEGIN and END.";
    }
}