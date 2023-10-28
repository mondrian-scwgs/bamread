# from __future__ import print_function

# import sys

import cython
import numpy as np
import pysam

from libc.stdint cimport (int8_t, int16_t, int32_t, int64_t, uint16_t,
                          uint32_t, uint64_t)

# from pysam.libcalignedsegment cimport AlignedSegment

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef _bamread(filename, uint32_t mapq=0, uint64_t required_flag=0, uint64_t filter_flag=1540):

    cdef:
        uint16_t flag
        int64_t start
        int64_t end
        uint32_t count = 0
        uint32_t nfound = 0
        long [::1] starts
        long [::1] ends
        int32_t [::1] chromosomes
        int8_t [::1] strands
        uint16_t [::1] flags
        uint32_t [::1] cbids

    samfile = pysam.AlignmentFile(filename, "rb")

    count = samfile.mapped + samfile.unmapped

    flags_arr = np.zeros(count, dtype=np.uint16)
    flags = flags_arr

    starts_arr = np.zeros(count, dtype=long)
    starts = starts_arr

    ends_arr = np.zeros(count, dtype=long)
    ends = ends_arr

    chromosomes_arr = np.zeros(count, dtype=np.int32)
    chromosomes = chromosomes_arr

    strands_arr = np.zeros(count, dtype=np.int8)
    strands = strands_arr

    cbids_arr = np.zeros(count, dtype=np.uint32)
    cbids = cbids_arr

    cbs = []
    for tag in samfile.header['CO']:
        if tag.startswith('CB:'):
            cbs.append(tag[3:])
    
    cbmap = {}
    for cbid, cb in enumerate(cbs):
        cbmap[cb] = cbid

    for a in samfile:
        flag = a.flag

        # https://broadinstitute.github.io/picard/explain-flags.html

        if a.mapping_quality < mapq:
            continue
        if (flag & required_flag) != required_flag:
            continue
        if (flag & filter_flag) != 0:
            continue
        if not a.has_tag('CB'):
            continue
        if a.get_tag('CB') not in cbmap:
            continue

        start = a.reference_start
        end = a.reference_end

        if start < 0 or end < 0:
            continue

        flags[nfound] = flag
        strands[nfound] = flag & 0x10
        chromosomes[nfound] = a.reference_id
        starts[nfound] = start
        ends[nfound] = end
        cbids[nfound] = cbmap[a.get_tag('CB')]

        nfound += 1

    chrs = samfile.references

    samfile.close()

    return (chromosomes_arr[:nfound], starts_arr[:nfound], ends_arr[:nfound], strands_arr[:nfound],
            flags_arr[:nfound], cbids[:nfound], chrs, cbs)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef _bamread_all(filename, uint32_t mapq=0, uint64_t required_flag=0, uint64_t filter_flag=1540):

    cdef:
        uint16_t flag
        int64_t start
        int64_t end
        uint32_t count = 0
        uint32_t nfound = 0
        long [::1] starts
        long [::1] ends
        int32_t [::1] chromosomes
        int8_t [::1] strands
        uint16_t [::1] flags
        uint32_t [::1] cbids
        # cdef AlignedSegment a


    samfile = pysam.AlignmentFile(filename, "rb")

    count = bamfile.mapped + bamfile.unmapped

    flags_arr = np.zeros(count, dtype=np.uint16)
    flags = flags_arr

    starts_arr = np.zeros(count, dtype=long)
    starts = starts_arr

    ends_arr = np.zeros(count, dtype=long)
    ends = ends_arr

    qstarts_arr = np.zeros(count, dtype=long)
    qstarts = qstarts_arr

    qends_arr = np.zeros(count, dtype=long)
    qends = qends_arr

    chromosomes_arr = np.zeros(count, dtype=np.int32)
    chromosomes = chromosomes_arr

    strands_arr = np.zeros(count, dtype=np.int8)
    strands = strands_arr

    cbids_arr = np.zeros(count, dtype=np.uint32)
    cbids = cbids_arr

    cbs = []
    for tag in samfile.header['CO']:
        if tag.startswith('CB:'):
            cbs.append(tag[3:])
    
    cbmap = {}
    for cbid, cb in enumerate(cbs):
        cbmap[cb] = cbid

    query_names, query_sequences, cigarstrings, query_qualities = [], [], [], []

    for a in samfile:
        flag = a.flag

        # https://broadinstitute.github.io/picard/explain-flags.html

        if a.mapping_quality < mapq:
            continue
        if (flag & required_flag) != required_flag:
            continue
        if (flag & filter_flag) != 0:
            continue
        if not a.has_tag('CB'):
            continue
        if a.get_tag('CB') not in cbmap:
            continue

        start = a.reference_start
        end = a.reference_end

        if start < 0 or end < 0:
            continue

        flags[nfound] = flag
        strands[nfound] = flag & 0x10
        chromosomes[nfound] = a.reference_id

        starts[nfound] = start
        ends[nfound] = end

        qstarts[nfound] = a.query_alignment_start
        qends[nfound] = a.query_alignment_end

        cbids[nfound] = cbmap[a.get_tag('CB')]

        query_names.append(a.query_name)
        query_sequences.append(a.query_sequence)
        query_qualities.append(a.query_qualities)
        cigarstrings.append(a.cigarstring)

        nfound += 1

    chrs = samfile.references

    samfile.close()

    return (chromosomes_arr[:nfound], starts_arr[:nfound], ends_arr[:nfound], strands_arr[:nfound],
            flags_arr[:nfound], cbids[:nfound], chrs, cbs, qstarts[:nfound], qends[:nfound], query_names, query_sequences, cigarstrings, query_qualities)
