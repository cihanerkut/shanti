#!/usr/bin/awk -f

# merge_featureCounts.awk
# Copyright (C) Cihan Erkut 2021
#
# This script merges 3 separate featureCounts output files.
#
# FeatureCounts is run 3 times for each BAM file for unstranded, forwardly-stranded
# and reversely-stranded counting. Since in the long-read mode featureCounts is
# single-threaded, running one thread per strand per BAM is more efficient. However,
# each instance outputs a separate TSV file with the counts as well as data from GTF
# and other values. This script merges those 3 files for each BAM file into one TSV file.
#
# It was originally written in Python using pandas, however, I realized later that it
# is overkill for such a simple task. Therefore I rewrote it in AWK. Unlike pandas, which
# depends on numpy, etc., AWK comes preinstalled in most Linux systems.
#

BEGIN {
    OFS = "\t"
    print "Geneid", "Chr", "Strand", "Length", "gene_name", "gene_type", "num", "num_fwd", "num_rev"
}

# We assume featureCounts inserts a commented first line with the command line. Here we skip
# that line as well the header line, which we already created in BEGIN block.
FNR < 3 {
    next
}

# Input file names must contain s0, s1 or s2 according to the -s parameter passed to
# featureCounts.
# Gene annotation is taken from the unstranded (s0) featureCounts file, assuming it is
# the same in other files as well.
FILENAME ~ /s0/ {
    n = FNR - 3
    split($2, Chr_array, ";")
    split($5, Strand_array, ";")

    Geneid[n] = $1
    Chr[n] = Chr_array[1]  # We assume that all features are on the same chromosome
    Strand[n] = Strand_array[1]  # We assume that all features are on the same strand
    Length[n] = $6
    gene_name[n] = $7
    gene_type[n] = $8
    num[n] = $9
}

# Forwardly-stranded
FILENAME ~ /s1/ {
    n = FNR - 3
    num_fwd[n] = $9
}

# Reversely-stranded
FILENAME ~ /s2/ {
    n = FNR - 3
    num_rev[n] = $9
}

# Build merged table
END {
    for (i = 0; i < n; i++)
        print Geneid[i], Chr[i], Strand[i], Length[i], gene_name[i], gene_type[i], num[i], num_fwd[i], num_rev[i]
}