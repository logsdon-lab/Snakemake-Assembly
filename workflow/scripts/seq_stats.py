#!/usr/bin/env python
import argparse
import os
import sys
import multiprocessing
from functools import partial
import pysam
import re
import numpy as np

ALN_PAT = re.compile(r".*\.(bam|sam|sam.gz|cram)")

# global var for inputs
args = None


def mp(func, *inargs, threads=8, chunksize=1000, msg="", **kwargs):
    with multiprocessing.Pool(threads) as pool:
        # convert to a tuple
        if len(inargs) > 1:
            iterator = []
            for arg in inargs:
                iterator.append(list(arg))
            iterator = zip(*iterator)
        else:
            iterator = inargs[0]
        iterator = list(iterator)
        print("Total to do:\t{}\t{}".format(len(iterator), func.__name__))

        # redifne function to pass kwargs
        new_func = partial(func, **kwargs)
        res = []
        for i, rtn in enumerate(pool.imap(new_func, iterator, chunksize=chunksize)):
            sys.stderr.write("\rDone with {} {}".format(i + 1, func.__name__))
            res.append(rtn)
        sys.stderr.write("\n")
        return res


def read_bam(f):
    try:
        names = []
        lengths = []
        bam = pysam.AlignmentFile(f, check_sq=False)
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or (not rec.is_secondary and not rec.is_supplementary):
                names.append(rec.query_name)
                lengths.append(len(rec.query_sequence))
        bam.close()
        sys.stderr.write("SAM/BAM read: {}\n".format(f))
        return (names, lengths)
    except ValueError:
        return None


def read_index(f):
    if os.path.exists(f + ".fai") and os.path.getmtime(f) > os.path.getmtime(
        f + ".fai"
    ):
        sys.stderr.write(f"Warning, index is older than {f}\n")
        # return(None)
    try:
        index = pysam.FastaFile(f)
        names = index.references
        lengths = index.lengths
        sys.stderr.write("Index read: {}\n".format(f))
        return (names, lengths)
    except (IOError, ValueError):
        return None


def read_fastx(f):
    try:
        names = []
        lengths = []
        with pysam.FastxFile(f, persist=False) as fastx:
            for rec in fastx:
                names.append(rec.name)
                lengths.append(len(rec.sequence))
        sys.stderr.write("Fastx read: {}\n".format(f))
        return (names, lengths)
    except IOError:
        return None


def read_bed(f):
    h = open(f)
    names = []
    lengths = []
    for line in h:
        if line[0] == "#":
            continue
        t = line.strip().split()
        lengths.append(int(t[2]) - int(t[1]))
        names.append(t[0])
    h.close()
    return (names, lengths)


def get_lengths(f):
    if re.match(ALN_PAT, f):
        rtn = read_bam(f)
        return rtn
    elif re.match(r".*\.bed", f):
        return read_bed(f)

    # try reading by index
    rtn = read_index(f)
    # try reading by fastx
    if rtn is None:
        rtn = read_fastx(f)
    return rtn


def calc_stats(lengths, qs, x=50, gsize=None):
    n = len(lengths)
    if gsize is None:
        total = sum(lengths)
    else:
        total = gsize
    lens = sorted(lengths)[::-1]
    mmax = lens[0]
    mmin = lens[-1]
    mean = total / n

    auN = sum([x * x for x in lens]) / total
    quantiles = np.quantile(lens, qs)

    # N50 stats
    count = 0
    for length in lens:
        count += length
        if count >= total * (x / 100):
            return (sum(lengths), n, mean, quantiles, mmin, mmax, length, auN)


def h_fmt(num):
    for unit in ["", "Kbp", "Mbp"]:
        if num < 1000.0:
            return "{:.2f}{}".format(num, unit)
        num /= 1000.0
    return "{:.2f}{}".format(num, "Gbp")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "infiles", nargs="+", help="fast{a,q}(.gz), sam, or bam inputs as a list."
    )
    parser.add_argument(
        "-t", "--threads", help="Number of threads to use", type=int, default=4
    )
    parser.add_argument(
        "-r", "--human", help="print human readable", action="store_true", default=False
    )
    parser.add_argument(
        "-q",
        "--quantiles",
        nargs="+",
        help="quantile(s) to calcaulte",
        type=float,
        default=[0.5],
    )
    parser.add_argument(
        "-g",
        help="calculate NG50, provide genome size, for hg use 3098794149",
        type=int,
        default=None,
    )
    parser.add_argument(
        "-d", help="store args.d as true if -d", action="store_true", default=False
    )
    args = parser.parse_args()

    infiles = []
    for path in args.infiles:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            infiles.append(path)
        else:
            sys.stderr.write("Skipping, because missing or empty: {}\n".format(path))

    threads = min(args.threads, len(args.infiles))
    str_qs = "\t".join(["{:g}%".format(x * 100) for x in args.quantiles])
    out = f"file\ttotalBp\tnSeqs\tmean\t{str_qs}\tmin\tmax\tN50\tauN\n"
    with multiprocessing.Pool(threads) as pool:
        for i, rtn in enumerate(pool.imap(get_lengths, infiles)):
            total, nseqs, mean, quantiles, mmin, mmax, N50, auN = calc_stats(
                rtn[1], args.quantiles, gsize=args.g
            )
            f = args.infiles[i]
            if args.human:
                str_out_qs = "\t".join([h_fmt(q) for q in quantiles])
                out_fmt = f"{f}\t{h_fmt(total)}\t{nseqs}\t{h_fmt(mean)}\t{str_out_qs}\t{h_fmt(mmin)}\t{h_fmt(mmax)}\t{h_fmt(N50)}\t{h_fmt(auN)}\n"
            else:
                str_out_qs = "\t".join(["{:g}".format(q) for q in quantiles])
                out_fmt = f"{f}\t{total}\t{nseqs}\t{mean}\t{str_out_qs}\t{mmin}\t{mmax}\t{N50}\t{auN}\n"
            out += out_fmt

    sys.stdout.write(out)
    sys.stdout.close()
