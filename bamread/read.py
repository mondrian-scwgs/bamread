from pathlib import Path
from typing import Tuple

import pandas as pd

from bamread.src.bamread import _bamread, _bamread_all  # type: ignore


def read_bam(
    f: Path, mapq: int = 0, required_flag: int = 0, filter_flag: int = 1540
) -> pd.DataFrame:
    chromosomes, starts, ends, strands, flags, cbids, chrmap, cbmap = _bamread(
        f, mapq, required_flag, filter_flag
    )

    chromosomes, ends, flags, starts, strands, cbs = _create_series(
        chrmap, chromosomes, ends, flags, starts, strands, cbmap, cbids,
    )

    return pd.DataFrame(
        {
            "Chromosome": chromosomes,
            "Start": starts,
            "End": ends,
            "Strand": strands,
            "Flag": flags,
            "CB": cbs,
        }
    )


def _create_series(
    chrmap, chromosomes, ends, flags, starts, strands, cbmap, cbids,
) -> Tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    chromosomes = pd.Categorical.from_codes(chromosomes, categories=chrmap)
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).replace({16: "-", 0: "+"}).astype("category")
    flags = pd.Series(flags)
    cbs = pd.Categorical.from_codes(cbids, categories=cbmap)
    return chromosomes, ends, flags, starts, strands, cbs


def read_bam_full(
    f: Path, mapq: int = 0, required_flag: int = 0, filter_flag: int = 1540
) -> pd.DataFrame:
    (
        chromosomes,
        starts,
        ends,
        strands,
        flags,
        cbids,
        chrmap,
        cbmap,
        qstarts,
        qends,
        query_names,
        query_sequences,
        cigarstrings,
        query_qualities,
    ) = _bamread_all(f, mapq, required_flag, filter_flag)

    chromosomes, ends, flags, starts, strands, cbs = _create_series(
        chrmap, chromosomes, ends, flags, starts, strands, cbmap, cbids,
    )
    qstarts = pd.Series(qstarts)
    qends = pd.Series(qends)
    query_names = pd.Series(query_names)
    query_sequences = pd.Series(query_sequences)
    cigarstrings = pd.Series(cigarstrings)
    query_qualities = pd.Series(query_qualities)

    return pd.DataFrame(
        {
            "Chromosome": chromosomes,
            "Start": starts,
            "End": ends,
            "Strand": strands,
            "Flag": flags,
            "QueryStart": qstarts,
            "QueryEnd": qends,
            "QuerySequence": query_sequences,
            "Name": query_names,
            "Cigar": cigarstrings,
            "Quality": query_qualities,
            "CB": cbs,
        }
    )
