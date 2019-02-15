import re
import sys
import os
import math
import numpy as np
from collections import defaultdict
from pyfaidx import Fasta

import gapCalc

GAP_PARAM_FILE = os.path.join(os.path.dirname(
    __file__), "axtChain_lastalgapCost.txt")

# A small wrapper to print to stderr


def eprint(*args, **kwargs):
    if verbose:
        print(*args, file=sys.stderr, **kwargs)


def ssign(strand):
    return 1 if strand=='+' else -1


def chromNum(chrom):
    m = re.search(r"[cC]hr(?P<num>[\dWXYZ]+)", chrom)
    if m is not None:
        if m.group('num') == 'X' or m.group('num') == 'W':
            num = 100
        elif m.group('num') == 'Y' or m.group('num') == 'Z':
            num = 101
        else:
            num = int(m.group('num'))
        return num
    else:
        return chrom


def get_scaffold_N_cumsum(scaffold, scaffolds_fasta):
    scaffolds = Fasta(scaffolds_fasta)
    return get_N_cumsum(scaffolds[scaffold][:].seq)


def get_N_cumsum(sequence):
    cum_sum = []
    current = 0
    for s in sequence:
        current += 1 if s == 'N' else 0
        cum_sum.append(current)
    return cum_sum


class Alignment(object):
    """
        An alignment as provided by the last tab output format
        According to the specification all positions are 0-based
    """

    def __init__(self, score, subject, spos, slen, sstrand, ssize,
                 query, qpos, qlen, qstrand, qsize):
        self._score = int(score)
        self._subject = subject
        self._spos = int(spos)
        self._slen = int(slen)
        self._sstrand = sstrand
        self._ssize = int(ssize)
        self._query = query
        self._qpos = int(qpos)
        self._qlen = int(qlen)
        self._qstrand = qstrand
        self._qsize = int(qsize)
        self._sindex = None
        self._qindex = None
        self._selected = True

    @classmethod
    def frombedpe(align, qchrom, qstart, qend, schrom, sstart, send,
                  name, score, qstrand, sstrand, qsize, ssize):
        "Initialize alignment from bedpe fields"
        qlen = int(qend) - int(qstart)
        slen = int(send) - int(sstart)
        if qstrand == "+":
            qpos = int(qstart)
        else:
            qpos = int(qsize) - int(qend)
        if sstrand == "+":
            spos = int(sstart)
        else:
            spos = int(ssize) - int(send)
        return align(score, schrom, spos, slen, sstrand, ssize,
                     qchrom, qpos, qlen, qstrand, qsize)
    # @classmethod
    # def fromlasttab( align, score, subject, spos, slen, sstrand, ssize,
    #                               query, qpos, qlen, qstrand, qsize):

    @property
    def name(self):
        return self._name

    @property
    def qname(self):
        return self._query

    @property
    def sname(self):
        return self._subject

    @property
    def qsize(self):
        return self._qsize

    @property
    def qlen(self):
        return self._qlen

    @property
    def slen(self):
        return self._slen

    @property
    def score_simple(self):
        return self._qlen + self._slen

    @property
    def score(self):
        return 5 * self.qlen
        # return 2*self._score

    @property
    def strand(self):
        return self._qstrand

    @property
    def sstrand(self):
        return self._sstrand

    @property
    def qsize(self):
        return self._qsize

    @property
    def ssize(self):
        return self._ssize

    @property
    def sstart(self):
        if self._sstrand == "+":
            return self._spos
        else:
            return self._ssize - self._spos - self._slen + 1

    @property
    def send(self):
        if self._sstrand == "+":
            return self._spos + self._slen
        else:
            return self._ssize - self._spos

    @property
    def qstart(self):
        if self._qstrand == "+":
            return self._qpos
        else:
            return self._qsize - self._qpos - self._qlen + 1

    @property
    def qend(self):
        if self._qstrand == "+":
            return self._qpos + self._qlen
        else:
            return self._qsize - self._qpos

    def __str__(self):
        values = [self.qname, self.qstart, self.qend,
                  self._qlen, self._qindex,
                  self.sname, self.sstart, self.send,
                  self._slen, self.strand, self._sindex,
                  self.score]
        # values = [self.qname, self.qstart, self.qend,
        #          self.strand,
        #          self.sname, self.sstart, self.send],self.score]
        return "\t".join(map(str, values))

    def shrtstr(self):
        values = [self.qname, self.qstart, self.qend,
                  self._qindex,
                  self.sname, self.sstart, self.send,
                  self.strand, self._sindex, self.score]
        return "\t".join(map(str, values))

    def tinystr(self):
        return "%s:%d-%d <--> %s:%d-%d %d" % (self.qname, self.qstart, self.qend,
                                              self.sname, self.sstart, self.send,
                                              self.score)

    def lasttab(self):
        values = [self.score,
                  self._subject, self._spos, self._slen,
                  self._sstrand, self._ssize,
                  self._query, self._qpos, self._qlen,
                  self._qstrand, self._qsize]
        return "\t".join(map(str, values))

    def bedpe(self):
        values = [
            self.qname, self.qstart, self.qend,
            self.sname, self.sstart, self.send,
            ".", self.score, self.strand, self.sstrand,
            self.qsize, self.ssize]
        return "\t".join(map(str, values))

    def hsp2string(self, format="bedpe"):
        if format == "bedpe":
            return self.bedpe()
        elif format == "lasttab":
            return self.lasttab()
        else:
            return self.__str__()

    @property
    def selected(self):
        return self._selected

    @property
    def qindex(self):
        return self._qindex

    @qindex.setter
    def qindex(self, index):
        self._qindex = index

    @property
    def sindex(self):
        return self._sindex

    @sindex.setter
    def sindex(self, index):
        self._sindex = index


def qdistance(a, b):
    """
       distance between two hsps on the query side
    """
    # hsps on different chromosomes are separated by an infinite distance
    if a.qname != b.qname:
        return math.inf
    # qindex must have been defined
    if a.qindex is None or b.qindex is None:
        eprint("qindex not defined")
        exit(1)
    # we first sort them to set a as the lowest
    a, b = (a, b) if a.qindex < b.qindex else (b, a)
    # now that everything is set we just return the distance
    return b.qstart - a.qend - 1


def q_real_nuc_distance(a, b, cumsum):
    """
       distance between two hsps on the query side
    """
    # hsps on different chromosomes are separated by an infinite distance
    if a.qname != b.qname:
        return math.inf
    # qindex must have been defined
    if a.qindex is None or b.qindex is None:
        eprint("qindex not defined")
        exit(1)
    # we first sort them to set a as the lowest
    a, b = (a, b) if a.qindex < b.qindex else (b, a)
    # now that everything is set we just return the distance
    distance = b.qstart - a.qend - 1
    if distance < 0:
        return 0
    if cumsum:
        # the number of 'N's in the interval
        num_N = cumsum[b.qstart - 2] - cumsum[a.qend]
        assert (distance >= num_N), "More Ns than nucleotides %d > %d" % (
            num_N, distance)
        distance -= num_N
    return distance


def sdistance(a, b):
    """
       distance between two hsps on the subject side
    """
    # hsps on different chromosomes are separated by an infinite distance
    if a.sname != b.sname:
        return math.inf
    # qindex must have been defined
    if a.sindex is None or b.sindex is None:
        eprint("qindex not defined")
        exit(1)
    # we first sort them to set a as the lowest
    a, b = (a, b) if a.sindex < b.sindex else (b, a)
    # now that everything is set we just return the distance
    distance = b.sstart - a.send - 1
    # in the case of overlap the distance is negative, we set it as 0
    if distance < 0:
        return 0
    return distance


class Chain(Alignment):
    """
       All coordinates are 0-based
    """
    def __init__(self, num, hsps,  score=-1, mixed=False):
        sorted_hsps = sorted(hsps, key=lambda x: x.qindex)
        first, last = (sorted_hsps[0], sorted_hsps[-1])
        self._hsps = sorted_hsps
        # an unstranded chain is a mixed chain
        if mixed:
            # When chained in an unstranded manner the chain has
            # no natural strand
            # we therefore define the strand as the majority strand
            # The order is not any more monotonic sstart and send have to
            # be computed in an "exhaustive" manner
            score_minus, score_plus = (0, 0)
            sstart, send = (math.inf, 0)
            qstart, qend = (math.inf, 0)
            for h in hsps:
                sstart = h.sstart if h.sstart < sstart else sstart
                send = h.send if h.send > send else send
                qstart = h.qstart if h.qstart < qstart else qstart
                qend = h.qend if h.qend > qend else qend
                if h.strand == "+":
                    score_plus += h.score
                else:
                    score_minus += h.score
            qstrand = "+" if score_plus > score_minus else "-"
            qsize = first.qsize
            if qstrand == "+":
                qpos = qstart
            else:
                qpos = qsize - qend
            spos = sstart
            qlen = qend - qstart
            slen = send - sstart
            sstrand = "+"
        else:
            if first.strand == "+":
                qpos = first.qstart
                spos = first.sstart
                slen = last.send - first.sstart
            else:
                qpos = last.qsize - last.qend
                spos = last.sstart
                slen = first.send - last.sstart
            qstrand = first.strand
            sstrand = first.sstrand
            qlen = last.qend - first.qstart

        qsize = first.qsize
        ssize = first.ssize
        query = first.qname
        subject = first.sname
        name = num

        super().__init__(name,
                         subject, spos, slen, sstrand, ssize,
                         query, qpos, qlen, qstrand, qsize)
        self._score = score
    # @property
    # def selected(self):
    #    return self._selected

    @property
    def score(self):
        return self._score

    def select(self):
        self._selected = True

    def deselect(self):
        self._selected = False

    def len(self):
        return len(self._hsps)

    def elements(self):
        return self._hsps

    def display(self):
        for a in self._hsps:
            print(a)


class AlignSet(object):
    def __init__(self, hsps, qstart_cutoff=None,
                 qend_cutoff=None,
                 size_cutoff=100,
                 score_cutoff=200,
                 qgap_cutoff=5000,
                 sgap_cutoff=5000,
                 cum_N_sum=None,
                 stranded=True,
                 debug=False
                 ):
        self._hsps = [hsp for hsp in hsps if (
            hsp.qlen > size_cutoff and hsp.score > score_cutoff)]
        if qstart_cutoff:
            self._hsps = [
                hsp for hsp in self._hsps if hsp.qstart > qstart_cutoff]
        if qend_cutoff:
            self._hsps = [
                hsp for hsp in self._hsps if hsp.qend < qend_cutoff]
        self.index()
        self._debug = debug
        self._qgap_cutoff = qgap_cutoff
        self._sgap_cutoff = sgap_cutoff
        self._stranded = stranded
        self._cum_N_sum = cum_N_sum
        self._gapCalc = gapCalc.gapCalc(GAP_PARAM_FILE)

    def set_chaining_parameters(self, qgap_cutoff=5000,
                                sgap_cutoff=5000,
                                cum_N_sum=None,
                                stranded=True):
        self._qgap_cutoff = qgap_cutoff
        self._sgap_cutoff = sgap_cutoff
        self._stranded = stranded
        self._cum_N_sum = cum_N_sum
        if self._debug > 2:
            print("qgap_cutoff = %d, sgap_cutoff = %d" %
                  (qgap_cutoff, sgap_cutoff))

    def debug(self, level):
        self._debug = level

    def index(self):
        qindex = []
        sindex = []
        for i, hsp in enumerate(sorted(self._hsps,
                                       key=lambda x: (chromNum(x.qname), x.qstart))):
            hsp.qindex = i
            qindex.append(hsp)
        for i, hsp in enumerate(sorted(self._hsps,
                                       key=lambda x: (chromNum(x.sname), x.sstart))):
            hsp.sindex = i
            sindex.append(hsp)
        self._qindex = qindex
        self._sindex = sindex

    def __getitem__(self, i):
        return self.qindex(i)

    def __iter__(self):
        return (hsp for hsp in self._hsps if hsp.selected)

    def len(self):
        return len([hsp for hsp in self._hsps if hsp.selected])

    def qindex(self, i):
        if i > len(self._qindex):
            eprint("Index out of range")
            exit(1)
        return self._qindex[i]

    def sindex(self, i):
        if i > len(self._sindex):
            eprint("Index out of range")
            exit(1)
        return self._sindex[i]

    def elements(self, filter=True, sorted=False):
        if filter:
            return [hsp for hsp in self._hsps if hsp.selected]
        else:
            return self._hsps

    def hsp2string(self, format="bedpe", filter=False, file=sys.stdout):
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            if not filter or a.selected:
                print(a.hsp2string(format), file=file)

    def qdisplay(self, format="bedpe", filter=False, file=sys.stdout):
        self.hsp2string(format, filter, file)

    def qprint(self):
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            print(a)

    def sdisplay(self):
        for a in sorted(self._hsps, key=lambda x: x.sindex):
            print(a)

    def qindices2elements(self, qindices):
        return [self.qindex(i) for i in qindices]

    def filter(self, scaffold_size, coverage=0.95, score=10000,
               length_cutoff=15000):
        for a in sorted(self._hsps, key=lambda x: x.score, reverse=True):
            a.deselect()
        cumsum = 0
        num_chain = 0
        for a in sorted(self._hsps, key=lambda x: x.score, reverse=True):
            if a.qlen < length_cutoff:
                continue
            a.select()
            num_chain += 1
            cumsum += a.qlen
            if a.score < score or cumsum / scaffold_size > coverage:
                break
        return num_chain, cumsum / scaffold_size

    def overlaps(self):
        intervals = []
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            intervals.append(create_interval_from_list(
                list(map(str, [a.qname, a.qstart - 1, a.qend, a.qname,
                               ".", a.strand]))))
        pybed_intervals = BedTool(intervals).sort()
        return

    def strict_chaining(self, qgap_cutoff=5000,
                        sgap_cutoff=5000,
                        cum_N_sum=None,
                        stranded=True):
        """
           Chain the hsps using a consecutive criteria
              - strict consecutive stranded hsps
              - consecutive unstranded hsps

        """
        self.set_chaining_parameters(qgap_cutoff,
                                     sgap_cutoff,
                                     cum_N_sum,
                                     stranded)
        if self._stranded:
            return self._strict_chaining_stranded()
        else:
            return self._strict_chaining_unstranded()

    def _strict_chaining_stranded(self):
        """
           Chain the hsps using strict consecutive criteria
              - hsps are sorted along both genomes
              - two hsps are chained if
                  - they are consecutive in both genomes
                  - they are in the same strand in both genomes
                  - the distance separating them is less than
                     qgap_cutoff (resp sgap_cutoff)
        """
        prev = None
        chains = []
        num = 0
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            if not prev:
                current_chain = [a]
            elif a.sname != prev.sname:
                score = self.score_lis(current_chain, self._cum_N_sum)
                chains.append(Chain(num, current_chain, score=score))
                current_chain = [a]
                num += 1
            else:
                if a.strand == "+":
                    if (a.sindex == prev.sindex + 1
                        and prev.strand == "+"
                        and a.qstart - prev.qend < self._qgap_cutoff
                            and a.sstart - prev.send < self._sgap_cutoff):
                        current_chain.append(a)
                    else:
                        score = self.score_lis(current_chain, self._cum_N_sum)
                        chains.append(Chain(num, current_chain, score=score))
                        current_chain = [a]
                        num += 1
                else:
                    if (a.sindex == prev.sindex - 1
                        and prev.strand == "-"
                        and a.qstart - prev.qend < self._qgap_cutoff
                            and prev.sstart - a.send < self._sgap_cutoff):
                        current_chain.append(a)
                    else:
                        score = self.score_lis(current_chain, self._cum_N_sum)
                        chains.append(Chain(num, current_chain, score=score))
                        current_chain = [a]
                        num += 1

            prev = a
        if current_chain and len(current_chain):
            score = self.score_lis(current_chain, self._cum_N_sum)
            # if self._debug:
            #     print(current_chain)
            chains.append(Chain(num, current_chain, score=score))
        return chains

    def _strict_chaining_unstranded(self):
        """
           Chain consecutive hsps
              - two hsps are merged if they pair the same chromosomes
                and if they are consecutive on qindex as well as sindex
        """
        prev = None
        chains = []
        current_chain = None
        num = 0
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            if not prev:
                current_chain = [a]
            elif a.sname != prev.sname:
                score = self.score_lis(current_chain, self._cum_N_sum)
                chains.append(Chain(num, current_chain,
                                    score=score, mixed=True))
                current_chain = [a]
                num += 1
            else:
                if self._debug:
                    self.gapscore(prev, a)
                if (abs(a.sindex - prev.sindex) == 1
                    and q_real_nuc_distance(a, prev,
                                            self._cum_N_sum) < self._qgap_cutoff
                        and sdistance(a, prev) < self._sgap_cutoff):
                    current_chain.append(a)
                    if self._debug:
                        print("chained")
                else:
                    score = self.score_lis(current_chain, self._cum_N_sum)
                    chains.append(Chain(num, current_chain,
                                        score=score, mixed=True))
                    current_chain = [a]
                    num += 1
            prev = a

        if current_chain and len(current_chain):
            score = self.score_lis(current_chain, self._cum_N_sum)
            chains.append(Chain(num, current_chain, score=score, mixed=True))
        return chains

    def backtrack(self, lis, back):
        # extract all subchain roots (sorting in descending order)
        candidates = [i for i in reversed(lis.argsort().tolist())]
        # The corresponding chains
        chain_indices = []
        while len(candidates) > 0:
            current = candidates[0]
            chain = []
            while current is not None:
                # the current node should not be part of a previous chain
                if current not in candidates:
                    break
                chain.append(current)
                candidates.remove(current)
                if self._debug > 2:
                    print(current, back[current])
                current = back[current]
            chain_indices.append(chain)
        return chain_indices

    def single_backtrack(self, endpoint, back):
        current = endpoint
        chain = []
        while current is not None:
            chain.append(current)
            if self._debug > 1:
                print(current, back[current])
            current = back[current]
        return chain

    def is_a_successor(self, a_i, a_j):
        """
            Returns true if a_i is a sucessor of a_j
               we must always have a_i.qindex > q_j.qindex
        """
        assert a_i.qindex > a_j.qindex, " i.qindex <= j.qindex here"
        if a_i.sname != a_j.sname:
            return False

        if self._stranded:
            # if stranded hsps should be both increasing if strand = "+"
            # and increasing on x and decreasing on y if strand = "-"
            if a_i.strand == a_j.strand:
                if ((a_i.strand == "+" and
                     a_i.sindex > a_j.sindex)
                    or
                    (a_i.strand == "-" and
                        a_i.sindex < a_j.sindex)):
                    return True
        else:
            # if not stranded chaining we just check for increasing
            # indices on both axes
            if a_i.sindex > a_j.sindex:
                return True
        return False

    def debug_lis(self, lis, i, j, gap_ij, s_ij):
        print("  l[%d]=%3.1f\tgap_ij=%3.1f\ts(%d,%d)=%3.1f l[%d]+s(%d,%d)=%3.1f and l[%d]=%3.1f" %
              (j, lis[j],
               gap_ij,
               j, i, s_ij,
               j, j, i, lis[j] + s_ij,
               i, lis[i]))

    def chaining(self, qgap_cutoff=5000,
                 sgap_cutoff=5000,
                 cum_N_sum=None,
                 stranded=True,
                 dropoff=-15000,
                 score_cutoff=10000):

        self.set_chaining_parameters(qgap_cutoff,
                                     sgap_cutoff,
                                     cum_N_sum,
                                     stranded)

        hsps = list(self._qindex)

        chains = []
        num = 0
        while len(hsps) > 0:
            elements, score = self.dropoff_lis_score(hsps,
                                                     dropoff=dropoff)
            if score < score_cutoff:
                break
            chain = Chain(num, elements, score=score)
            if self._debug:
                print("chain: ", chain)
                chain.display()
            chains.append(chain)
            # we remove the chained alignments from the list
            for elt in elements:
                hsps.remove(elt)
            num += 1

        return chains

    def dropoff_chain(self, qgap_cutoff=5000,
                      sgap_cutoff=5000,
                      cum_N_sum=None,
                      stranded=True,
                      dropoff=-15000):

        self.set_chaining_parameters(qgap_cutoff,
                                     sgap_cutoff,
                                     cum_N_sum,
                                     stranded)

        hsps = list(self._qindex)
        elements, score = self.dropoff_lis_score(hsps, dropoff=dropoff)

        return elements, score

    def dropoff_lis_score(self, hsps, dropoff=-15000):

        a = sorted(hsps, key=lambda x: x.qindex)
        n = len(a)

        # Declare the list (array) for LIS and initialize LIS
        # values for all indexes
        # We initialize the lis array with the segment score
        lis = np.array([0] * n)
        for i in range(0, n):
            lis[i] = a[i].score

        # Declare the list (array) for backtrack and initialize to None
        back = [None] * n

        # we start at the highest score
        candidates = [i for i in reversed(lis.argsort().tolist())]
        start = candidates[0]

        if self._debug > 0:
            print("start=%d %s" % (start, a[start]))

        blacklist = defaultdict()

        # we extend to the right
        # Compute optimized LIS values extending to the right
        best_score = lis[start]
        best_end = start
        skipped = 0
        max_i = 0
        for i in range(start + 1, n):
            max_sij = -math.inf
            if self.is_a_successor(a[i], a[start]):
                skipped = 0
                for j in range(start, i):
                    if j in blacklist:
                        continue
                    gap_ij = self.gapscore(a[j], a[i])
                    s_ij = a[i].score - gap_ij
                    if self._debug > 1:
                        self.debug_lis(lis, i, j, gap_ij, s_ij)
                    if lis[j] + s_ij > lis[i]:
                        lis[i] = lis[j] + s_ij
                        back[i] = j
                if self._debug > 0:
                    print("right l[%d]=%3.1f gap_ij = %3.1f prev = %r hsp=%s"
                          % (i, lis[i], gap_ij, back[i], a[i].tinystr()))
                if back[i] is None:
                    blacklist[i] = 1
                elif lis[i] > max_i:
                    max_i = lis[i]
                    best_end = i
            else:
                blacklist[i] = 1
                skipped += 1
            if skipped > 100:
                break

        # we need a first backtracking here
        qindices_right = self.single_backtrack(best_end, back)
        # last element of the chain should be the start
        #  assert qindices_right[-1] == start

        # now we extend to the left
        # Compute optimized LIS values extending to the left
        # we proceed here in decreasing i values (hence decreasing
        # qindex values)
        best_score = lis[start]
        best_end = start
        skipped = 0
        max_i = 0
        for i in range(start - 1, -1, -1):
            if self.is_a_successor(a[start], a[i]):
                skipped = 0
                for j in range(start, i, -1):
                    if j in blacklist:
                        continue
                    gap_ij = self.gapscore(a[i], a[j])
                    s_ij = a[i].score - gap_ij
                    if self._debug > 2:
                        self.debug_lis(lis, i, j, gap_ij, s_ij)
                    if lis[j] + s_ij > lis[i]:
                        lis[i] = lis[j] + s_ij
                        back[i] = j
                if back[i] is None:
                    blacklist[i] = 1
                elif lis[i] > max_i:
                    max_i = lis[i]
                    best_end = i
                if self._debug > 0:
                    print("left l[%d]=%3.1f gap_ij = %3.1f prev = %r hsp=%s"
                          % (i, lis[i], gap_ij, back[i], a[i].tinystr()))
            else:
                blacklist[i] = 1
                skipped += 1
            if skipped > 100:
                break

        # we need a second backtracking here
        qindices_left = self.single_backtrack(best_end, back)

        # last element of the chain should be the start
        # if not we empty the backtracking procedure
        # assert qindices_left[-1] == start

        if self._debug > 1:
            print(qindices_right + qindices_left)

        # from relative to abslute indices
        qindices = defaultdict()
        for i in qindices_right + qindices_left:
            qindices[a[i].qindex] = 1

        elements = self.qindices2elements(list(qindices.keys()))
        score = self.score_lis(elements, self._cum_N_sum)

        return elements, score
    # end of dropoff_lis_score

    def gapscore(self, a_j, a_i):
        """
            Compute the gap score of chaining j to i, essentially
                qgap(j,i) + sgap(j,i)
            cum_N_sum : is used to correct fo query gaps (N's are omitted)
            stranded : i and j must be on the same strand
        """

        q_distance = q_real_nuc_distance(a_i, a_j, self._cum_N_sum)
        s_distance = sdistance(a_i, a_j)

        if self._debug > 2:
            print("    " + a_j.shrtstr())
            print("    " + a_i.shrtstr())
            print("    qdist=%2.0f sdist=%2.0f" % (q_distance, s_distance))
            print("    qgap_cutoff=%d sgap_cutoff=%d" % (self._qgap_cutoff,
                                                         self._sgap_cutoff))

        if (q_distance > self._qgap_cutoff or
                s_distance > self._sgap_cutoff):
            return math.inf

        if self._stranded:
            # if stranded hsps should be both increasing if strand = "+"
            # and increasing on x and decreasing on y if strand = "-"
            if a_i.strand == a_j.strand:
                if ((a_i.strand == "+" and
                     a_i.sindex > a_j.sindex)
                    or
                    (a_i.strand == "-" and
                        a_i.sindex < a_j.sindex)):
                    return self.calc_gapscore(q_distance, s_distance)
        else:
            # if not stranded chaining we just check for increasing
            # indices on both axes
            if a_i.sindex > a_j.sindex:
                return self.calc_gapscore(q_distance, s_distance)
        # previous clauses where false, we return -inf
        return math.inf
    # end of score function

    def score_lis(self, hsps, cum_N_sum=None):
        sorted_hsps = sorted(hsps, key=lambda x: x.qindex)
        prev = sorted_hsps[0]
        score = prev.score
        if self._debug > 2:
            print("%d" % (score))
        for i in range(1, len(sorted_hsps)):
            curr = sorted_hsps[i]
            q_distance = q_real_nuc_distance(curr, prev, cum_N_sum)
            s_distance = sdistance(curr, prev)
            score += curr.score - self.calc_gapscore(q_distance, s_distance)
            if self._debug > 2:
                print("%d\t%d\t%d" % (curr.score, q_distance, s_distance))
            prev = curr
        return score

    def gapscore_simple(self, dq, dt):
        return dq + dt

    def calc_gapscore(self, dq, dt):
        # return self.gapscore_simple(dq, dt)
        return self._gapCalc.gapCalcCost(dq, dt)


def valid_chromosome(chrom,regexp):
    p = re.compile(regexp)
    return p.match(chrom) is not None


def read_alignfile(alignfile,
                   min_scaffold_size=20000,
                   min_align_length=200,
                   autosomes=False,
                   validchromreg='((C|c)hr|Scaffold)?(\d+|Z|W)$'):
    """
       Reading the last tab file
    """
    i = 0
    hsps = defaultdict(list)
    scaffolds = defaultdict(list)
    # valid chromosomes
    chromosomes = [str(x) for x in range(1, 29)]
    if not autosomes:
        chromosomes += ['W', 'Z']
    with open(alignfile) as fin:
        for index, line in enumerate(fin):
            if not line.startswith('#'):
                fields = line.rstrip().split()
                # the alignment object
                align = Alignment(*fields[0:11])
                # only alignments passing the filters are kept
                if (align.qsize < min_scaffold_size or
                    align.qlen < min_align_length or
                    not valid_chromosome(align.sname, validchromreg)):
                    continue
                # keep a dictionnary of scaffolds and associated sizes
                if align.qname not in scaffolds:
                    scaffolds[align.qname] = align.qsize
                hsps[align.qname].append(align)
    return hsps, scaffolds


def read_psl_alignfile(alignfile,
                       min_scaffold_size=20000,
                       min_align_length=200,
                       autosomes=False):
    """
       Reading the last psl file
    """
    i = 0
    hsps = defaultdict(list)
    scaffolds = defaultdict(list)
    # valid chromosomes
    chromosomes = [str(x) for x in range(1, 29)]
    if not autosomes:
        chromosomes += ['W', 'Z']
    with open(alignfile) as fin:
        for index, line in enumerate(fin):
            if not line.startswith('#'):
                fields = line.rstrip().split()
                # the alignment object
                align = Alignment(*fields[0:11])
                # only alignments passing the filters are kept
                if (align.qsize < min_scaffold_size or
                    align.qlen < min_align_length or
                        align.sname not in chromosomes):
                    continue
                # keep a dictionnary of scaffolds and associated sizes
                if align.qname not in scaffolds:
                    scaffolds[align.qname] = align.qsize
                hsps[align.qname].append(align)
    return hsps, scaffolds


def chaining(elements, qgap_cutoff=5000,
             sgap_cutoff=5000,
             cum_N_sum=None,
             stranded=True,
             dropoff=-15000,
             score_cutoff=10000,
             debug=False):

    hsps = list(elements)

    chains = []
    num = 0
    while len(hsps) > 0:
        align_set = AlignSet(hsps)
        align_set.debug(debug)
        elements, score = align_set.dropoff_chain(qgap_cutoff,
                                                  sgap_cutoff,
                                                  cum_N_sum,
                                                  stranded,
                                                  dropoff)
        if score < score_cutoff:
            break
        chain = Chain(num, elements, score=score)
        chains.append(chain)
        # we remove the chained alignments from the list
        for elt in elements:
            hsps.remove(elt)
        num += 1

    return chains


class AlignSetOld(object):
    def strict_chaining_old(self, qgap_cutoff=5000,
                            sgap_cutoff=5000,
                            cum_N_sum=None,
                            stranded=True):
        """
           Chain the hsps using strict consecutive criteria
              - hsps are sorted along both genomes
              - two hsps are chained if
                  - they are consecutive in both genomes
                  - they are in the same strand in both genomes
                  - the distance separating them is less than
                     qgap_cutoff (resp sgap_cutoff)
        """
        prev = None
        chains = []
        num = 0
        for a in sorted(self._hsps, key=lambda x: x.qindex):
            if not prev:
                current_chain = [a]
            elif a.sname != prev.sname:
                score = self.score_lis(current_chain, cum_N_sum)
                chains.append(Chain(num, current_chain, score=score))
                current_chain = [a]
                num += 1
            else:
                if a.strand == "+":
                    if (a.sindex == prev.sindex + 1
                        and prev.strand == "+"
                        and a.qstart - prev.qend < qgap_cutoff
                            and a.sstart - prev.send < sgap_cutoff):
                        current_chain.append(a)
                    else:
                        score = self.score_lis(current_chain, cum_N_sum)
                        chains.append(Chain(num, current_chain, score=score))
                        current_chain = [a]
                        num += 1
                else:
                    if (a.sindex == prev.sindex - 1
                        and prev.strand == "-"
                        and a.qstart - prev.qend < qgap_cutoff
                            and prev.sstart - a.send < sgap_cutoff):
                        current_chain.append(a)
                    else:
                        score = self.score_lis(current_chain, cum_N_sum)
                        chains.append(Chain(num, current_chain, score=score))
                        current_chain = [a]
                        num += 1
            prev = a

        score = self.score_lis(current_chain, cum_N_sum)
        chains.append(Chain(num, current_chain, score=score))
        return chains
    # lis returns length of the longest increasing subsequence
    # in arr of size n

    def lischain(self):
        arr = sorted(self._hsps, key=lambda x: (chromNum(x.qname), x.qstart))
        n = len(arr)

        # Declare the list (array) for LIS and initialize LIS
        # values for all indexes
        lis = np.array([1] * n)
        back = [None] * n
        # Compute optimized LIS values in bottom up manner
        for i in range(1, n):
            for j in range(0, i):
                if (arr[i].strand == arr[j].strand and
                        arr[i].sname == arr[j].sname):
                    if arr[i].strand == "+":
                        if arr[i].sindex > arr[j].sindex and lis[i] < lis[j] + 1:
                            lis[i] = lis[j] + 1
                            back[i] = j
                    else:
                        if arr[i].sindex < arr[j].sindex and lis[i] < lis[j] + 1:
                            lis[i] = lis[j] + 1
                            back[i] = j
        chain_indices = self.backtrack(lis, back)
        chains = []
        num = 0
        for qindices in chain_indices:
            chains.append(Chain(num, self.qindices2elements(qindices)))
            num += 1
        return chains
    # end of lis function

    def lis_simple_score(self):

        a = sorted(self._hsps, key=lambda x: x.qindex)
        n = len(a)

        # Declare the list (array) for LIS and initialize LIS
        # values for all indexes
        lis = np.array([1] * n)
        # Declare the list (array) for backtrack and initialize to None
        back = [None] * n

        # Compute optimized LIS values in bottom up manner
        for i in range(1, n):
            for j in range(0, i):
                s_ij = self.simple_score(a[i], a[j])
                if lis[j] + s_ij > lis[i]:
                    lis[i] = lis[j] + s_ij
                    back[i] = j

        chain_indices = self.backtrack(lis, back)
        chains = []
        num = 0
        for qindices in chain_indices:
            # The first index is always the index with the chain score
            score = len(qindices)
            chains.append(Chain(num,
                                self.qindices2elements(qindices),
                                score=score))
            num += 1
        return chains

    # end of lisscore function
    def simple_score(self, a_i, a_j):

        q_distance = q_real_nuc_distance(a_i, a_j, self._cum_N_sum)
        s_distance = sdistance(a_i, a_j)

        if (a_i.strand != a_j.strand or
            q_distance > qgap_cutoff or
                s_distance > sgap_cutoff):
            return -math.inf

        if ((a_i.strand == "+" and
             a_i.sindex > a_j.sindex)
            or
            (a_i.strand == "-" and
             a_i.sindex < a_j.sindex)):
            return 1

        return -math.inf

    def lis_score(self):
        """
            Compute the highest scoring increasing sequencing stopping
            at i for all i.
        """
        a = self._qindex
        n = len(a)

        # Declare the list (array) for LIS and initialize LIS
        # values for all indexes
        # We initialize the lis array with the hsp score
        lis = np.array([0] * n)
        for i in range(0, n):
            lis[i] = a[i].score

        # Declare the list (array) for backtrack and initialize to None
        back = [None] * n

        # Compute optimized LIS values in bottom up manner
        for i in range(1, n):
            for j in range(0, i):
                gap_ij = self.gapscore(a[j], a[i])
                s_ij = a[i].score - gap_ij
                if lis[j] + s_ij > lis[i]:
                    lis[i] = lis[j] + s_ij
                    back[i] = j

        chain_indices = self.backtrack(lis, back)
        chains = []
        num = 0
        for qindices in chain_indices:
            elements = self.qindices2elements(qindices)
            score = self.score_lis(elements, cum_N_sum)
            chains.append(Chain(num, elements, score=int(score)))
            num += 1
        return chains

    def dropoff_lis_score_debug(self, hsps, qgap_cutoff=5000,
                                sgap_cutoff=5000,
                                cum_N_sum=None,
                                stranded=True,
                                dropoff=-15000):

        a = sorted(hsps, key=lambda x: x.qindex)
        n = len(a)

        # Declare the list (array) for LIS and initialize LIS
        # values for all indexes
        # We initialize the lis array with the segment score
        lis = np.array([0] * n)
        for i in range(0, n):
            lis[i] = a[i].score

        # Declare the list (array) for backtrack and initialize to None
        back = [None] * n

        # we start at the highest score
        candidates = [i for i in reversed(lis.argsort().tolist())]
        start = candidates[0]

        # we extend to the right
        # Compute optimized LIS values extending to the right
        best_score = lis[start]
        best_end = start
        for i in range(start + 1, n):
            max_sij = -math.inf
            if self.is_a_successor(a[i], a[start]):
                for j in range(start, i):
                    gap_ij = self.gapscore(a[j], a[i], qgap_cutoff,
                                           sgap_cutoff,
                                           cum_N_sum,
                                           stranded)
                    s_ij = a[i].score - gap_ij
                    if self._debug > 1:
                        self.debug_lis(lis, i, j, gap_ij, s_ij)
                    if lis[j] + s_ij > lis[i]:
                        lis[i] = lis[j] + s_ij
                        max_sij = max(s_ij, max_sij)
                        back[i] = j
                if self._debug > 1:
                    print("right l[%d]=%3.1f max_sij = %3.1f prev = %r"
                          % (i, lis[i], max_sij, back[i]))
                if max_sij < dropoff:
                    break
                best_end = i

        # we need a first backtracking here
        qindices_right = self.single_backtrack(best_end, back)
        # last element of the chain should be the start
        assert qindices_right[-1] == start

        # we extend to the left
        # Compute optimized LIS values extending to the left
        # we proceed here in decreasing i values (hence qindex values)
        best_score = lis[start]
        best_end = start
        for i in range(start - 1, -1, -1):
            if self.is_a_successor(a[start], a[i]):
                max_sij = -math.inf
                for j in range(start, i, -1):
                    gap_ij = self.gapscore(a[i], a[j], qgap_cutoff,
                                           sgap_cutoff,
                                           cum_N_sum,
                                           stranded)
                    s_ij = a[i].score - gap_ij
                    if self._debug > 0:
                        print("s_i=%d" % (a[i].score))
                        self.debug_lis(lis, i, j, gap_ij, s_ij)
                    if lis[j] + s_ij > lis[i]:
                        lis[i] = lis[j] + s_ij
                        max_sij = max(s_ij, max_sij)
                        back[i] = j
                if self._debug > 1:
                    print("left l[%d]=%3.1f max_sij = %3.1f prev = %r"
                          % (i, lis[i], max_sij, back[i]))
                if max_sij < dropoff:
                    break
                best_end = i

        # we need a second backtracking here
        qindices_left = self.single_backtrack(best_end, back)
        # last element of the chain should be the start
        assert qindices_left[-1] == start
        # remove this last element already in qindices_right
        qindices_left.pop()

        elements = self.qindices2elements(qindices_right + qindices_left)
        score = self.score_lis(elements, cum_N_sum)

        return elements, score
