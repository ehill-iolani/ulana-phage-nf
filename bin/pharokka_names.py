#!/usr/bin/env python3
"""
Keeps Pharokka from giving two genes the same name, by renaming the contigs
it sees and putting the names back in its outputs afterwards.

Pharokka 1.10.1 names each protein in its working FASTA by joining the
contig name to a gene number with nothing in between, and the number counts
genes across the whole assembly, not per contig (processes.py:
`str(contig_list[i]) + str(i)`). With Flye's names, gene 825 on contig_8 and
gene 25 on contig_88 are both "contig_8825", and post-processing stops with
`ValueError: Duplicate key 'contig_8825'`. That needs more than one contig and
a hundred-odd genes, which is what a host-contaminated plaque prep gives
(WGS_MaterMagnus, PRJNA1502461 validation 2026-11-12: 181 contigs, 1.1 Mb).
Flye's names change from run to run, so the same reads can pass one run and
fail the next.

mark: every contig name that ends in a digit gets an "x" on the end
(contig_8 -> contig_8x). The digits at the end of a joined name are then
exactly the gene number, so no two can match ("contig_8x825" vs
"contig_88x25"). The names it changed go to a list, one per line.

restore: takes the "x" back off those names, and only those, wherever one
stands as a whole name in Pharokka's outputs (GFF, GenBank, tables, FASTA
headers), so every file uses the assembly's own contig names again. The files
are rewritten in place.

Pure standard library. Usage:

  pharokka_names.py mark --in assembly.fasta --out pharokka_in.fasta --names renamed.txt
  pharokka_names.py restore --names renamed.txt FILE [FILE ...]
"""
import argparse
import os
import re
import sys

SUFFIX = "x"


def mark(src, dst, names_path):
    renamed = []
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                name, sep, rest = header.partition(" ")
                if name[-1:].isdigit():
                    renamed.append(name)
                    name += SUFFIX
                line = f">{name}{sep}{rest}\n"
            fout.write(line)
    with open(names_path, "w") as f:
        f.writelines(n + "\n" for n in renamed)
    return renamed


def restore_pattern(names):
    # longest first, so contig_88x is never read as contig_8 + "8x"; a name
    # only counts when it isn't part of a longer word on either side
    alts = "|".join(re.escape(n) for n in sorted(names, key=len, reverse=True))
    return re.compile(
        rf"(?<![A-Za-z0-9_.])({alts}){re.escape(SUFFIX)}(?![A-Za-z0-9])")


def restore(names_path, files):
    with open(names_path) as f:
        names = [n.strip() for n in f if n.strip()]
    if not names:
        return 0
    pattern = restore_pattern(names)
    changed = 0
    for path in files:
        # the pharokka* glob also matches the database folder staged into the
        # task dir (pharokka_v1.8.0_databases); only plain files are outputs
        if os.path.isdir(path):
            continue
        try:
            with open(path, encoding="utf-8") as f:
                text = f.read()
        except FileNotFoundError:
            continue  # an unmatched shell glob, e.g. terL.f* when no terL was found
        except UnicodeDecodeError:
            continue  # not a text file; Pharokka writes none today
        new = pattern.sub(r"\1", text)
        if new != text:
            with open(path, "w", encoding="utf-8") as f:
                f.write(new)
            changed += 1
    return changed


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mark", help="rename contigs before Pharokka")
    m.add_argument("--in", dest="src", required=True)
    m.add_argument("--out", dest="dst", required=True)
    m.add_argument("--names", required=True)
    r = sub.add_parser("restore", help="put the names back in Pharokka's outputs")
    r.add_argument("--names", required=True)
    r.add_argument("files", nargs="*")
    a = p.parse_args(argv)
    if a.cmd == "mark":
        n = len(mark(a.src, a.dst, a.names))
        print(f"pharokka_names: {n} contig name(s) marked", file=sys.stderr)
    else:
        n = restore(a.names, a.files)
        print(f"pharokka_names: names restored in {n} file(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
