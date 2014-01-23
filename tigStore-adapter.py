#!/usr/bin/env python

import os
import sys
import subprocess
import shlex
import string
from pbcore.io import FastaIO
from multiprocessing import Pool

store = sys.argv[1]
utgid = int(sys.argv[2])

gkp_store = "%s.gkpStore" % store
tig_store = "%s.tigStore" % store

rc = string.maketrans('actgACTG', 'tgacTGAC')


class FragRec(object):
    def __init__(self, frgid, start, end):
        self.frgid = int(frgid)
        self.start = int(min(start, end))
        self.end = int(max(start, end))
        self.rev = start - end > 0


def process_unitig(unitig_id):

    # more self-descriptive id
    utgid = "unitig_%d" % unitig_id

    tigStore_args = shlex.split("tigStore -g %s -t %s 1 -d fr -u %d"
        % (gkp_store, tig_store, unitig_id))
    frags = []
    max_coor = 0
    out = subprocess.check_output(tigStore_args)
    out = out.split("\n")

    # save fragment ids to a file
    idfh = open("%s.frgids" % utgid, "w")
    for frag in out:
        """FRG 1453 179419,182165"""
        frag = frag.strip()
        if len(frag) == 0:
            continue
        frag = frag.replace(",", " ")
        frag = frag.strip().split()
        frag_id = frag[1]
        idfh.write("%s\n" % frag_id)

        frags.append(FragRec(frag_id, int(frag[2]), int(frag[3])))
        max_coor = max(int(frag[2]), int(frag[3]), max_coor)

    idfh.close()

    # get the fragment sequences
    gkcall = "gatekeeper -dumpfasta %s_frgs -iid %s.frgids %s"
    subprocess.call(shlex.split(gkcall % (utgid, utgid, gkp_store)))

    # load fragment seqs into memory
    seq_db = {}
    seqF = FastaIO.FastaReader("%s_frgs.fasta" % utgid)
    for r in seqF:
        frgid = int(r.name.split()[0].split(',')[1])
        seq_db[frgid] = r.sequence

    seq_array = ["N"] * max_coor

    # build the unitig sequence
    utgseq = ''
    for f in frags:
        sequence = seq_db[f.frgid]
        if f.rev:
            sequence = sequence.translate(rc)[::-1]

        for p in range(f.start, f.end):
            if seq_array[p] == "N":
                try:
                    seq_array[p] = sequence[p - f.start]
                except:
                    pass

    utgseq = "".join(seq_array)

    # create a simple 'layout' file for consumption by pbutgcns
    with open('%s.lay' % utgid, 'w') as fh:
        # first line is the initial consensus
        fh.write("%s %s\n" % (utgid, utgseq))
        # the rest are the fragments and their offsets
        for f in frags:
            qseq = seq_db[f.frgid]
            if f.rev:
                qseq = qseq.translate(rc)[::-1]

            fh.write('%d %d %d %s\n' % (f.frgid, f.start, f.end, qseq))

process_unitig(utgid)

"""
args = shlex.split("tigStore -g %s -t %s 1 -D unitiglist"
    % (gkp_store, tig_store))
out = subprocess.check_output(args)
out = out.split("\n")

unitig_id_list = []
for l in out:
    l = l.strip().split()
    if len(l) == 0:
        continue
    if l[0] == "maID":
        continue
    process_unitig(int(l[0]))
"""
