import Bio.pairwise2
import ab1
import sys
import re

def as_key(s):
    return tuple(sorted(set(s)))

iupac = {('A','C'): 'M',
         ('A','G'): 'R',
         ('A','T'): 'W',
         ('C','G'): 'S',
         ('C','T'): 'Y',
         ('G','T'): 'K',
         ('A',): 'A',
         ('C',): 'C',
         ('T',): 'T',
         ('G',): 'G',
         ('A','C','G'): 'V',
         ('A','C','T'): 'H',
         ('A','G','T'): 'D',
         ('C','G','T'): 'B',
         ('A','C','G','T'): 'N'}

iupac_table = {}
for k1,v1 in iupac.iteritems():
    for k2,v2 in iupac.iteritems():
        new_key = as_key([v1,v2])
        new_base = iupac[as_key(set(k1).union(k2))]
        iupac_table[new_key] = new_base
for i in iupac.itervalues():
    iupac_table[(i,)] = i
    iupac_table[('-',i)] = i
    iupac_table[('.',i)] = i
iupac_table[tuple()] = 'N'

def mask_poor_bases(s, threshold=40):
    new_seq = ''
    for i in range(len(s)):
        if s.conf[i] > threshold:
            new_seq += s.sequence[i]
        else:
            new_seq += 'N'
    assert len(s) == len(new_seq)
    return new_seq
    

def merge(s1, s2, min_segment_len=20, quality_threshold=40):
    # s1 and s2 are SequenceWithConfidence objects containing the two 
    # reads. The preconditions are:
    #   * s1 is a noisy subsequence of some template sequence x,
    #     and s2 is a noisy subsequence of revcomp(x).
    #   * s1 and s2 both consist of two low quality segments flanking
    #     a (possibly empty) high quality segment.
    #   * s1 and s2 are words over the alphabet AGCTN.
    #   * s1 and s2 may contain incorrect bases and indels.

    # First we extract the high quality segment in the middle.
    # Extremely short high quality segments are not useful to us, so
    # we put a threshold on the lowest length of interest to us.  If
    # we cannot get a long enough segment out of either sequence, our
    # sequence estimation fails.  If only one sequence has a long
    # enough high quality segment, we use its high quality segment and
    # ignore the other sequence entirely.
    h1 = s1.high_quality_segment(threshold=quality_threshold)
    h2 = s2.high_quality_segment(threshold=quality_threshold).revcomp()
    if len(h1) < min_segment_len and len(h2) < min_segment_len:
        return None
    elif len(h1) < min_segment_len:
        # But we don't want to accept low quality bases.  We set them
        # back to N if the quality is below threshold.  If there are
        # more than 15% N's in the sequence, we fail.
        masked = mask_poor_bases(h2, quality_threshold)
        if len([x for x in masked if x == 'N'])/float(len(masked)) > 0.15:
            return None
        else:
            return (masked,
                    ''.join(['-' for i in range(len(h2))]), 
                    h2.sequence)
    elif len(h2) < min_segment_len:
        masked = mask_poor_bases(h1, quality_threshold)
        if len([x for x in masked if x == 'N'])/float(len(masked)) > 0.15:
            return None
        else:
            return (masked,
                    h1.sequence,
                    ''.join(['-' for i in range(len(h1))]))
    else:
        pass

    # At this point h1 and h2 are a pair of high quality segments of
    # adequate length.  We align them locally, imposing a gap penalty
    # to keep each sequence clumped together.  This produces a1 and
    # a2, which are aligned, have the same length, and are words over
    # the alphabet ACGTN-.
    alignment = Bio.pairwise2.align.localxs(h1.sequence, h2.sequence, -5, -10)[0]
    # Bio.pairwise2 doesn't understand SequenceWithConfidence objects,
    # so we have to reapply the confidences in order to proceed.
    assert(isinstance(alignment[0], str))
    assert(isinstance(alignment[1], str))
    a1 = ab1.reapply_confidences(alignment[0], h1)
    a2 = ab1.reapply_confidences(alignment[1], h2)

    # Finally we choose a base at each position by a two pass voting
    # algorithm.  In the first round, each sequence votes on whether
    # this base should be skipped.  '-' is a vote to skip, anything
    # else is a vote to insert.  A vote with low confidence is ignored
    # if there is a vote to skip.  Upon skipping, an 'N' is
    # accumulated, and upon the next vote not to skip, the number of
    # accumulated N's is counted.  If there are at least 2 N's, then
    # the N's are inserted in the sequence.  Otherwise the accumulated
    # N's are abandoned.

    # If there is a second round, that is, if the base is not skipped,
    # then there is a vote on what base it is.  A low confidence base
    # leads to an ignored vote.  High confidence bases are combined
    # into their IUPAC codes, ignoring any vote of '-'.

    # This all applies during the overlap of the sequences.  In end
    # regions where there is only one sequence, we need to behave
    # differently.  In this case, we replace the -'s at the beginning
    # and end with .'s, representing uncalled bases.  A dot does not
    # vote for skipping, but also does not vote for a base.
    offset = len(re.split(r'[^-]', a1.sequence, maxsplit=1)[0])
    if offset > 0:
        a1.sequence = ''.join(['.' for i in range(offset)]) + a1.sequence[offset:]
    offset = len(re.split(r'[^-]', a2.sequence, maxsplit=1)[0])
    if offset > 0:
        a2.sequence = ''.join(['.' for i in range(offset)]) + a2.sequence[offset:]
    offset = len(re.split(r'[^-]', a1.sequence[::-1], maxsplit=1)[0])
    if offset > 0:
        a1.sequence = a1.sequence[:(-offset)] + ''.join(['.' for i in range(offset)])
    offset = len(re.split(r'[^-]', a2.sequence[::-1], maxsplit=1)[0])
    if offset > 0:
        a2.sequence = a2.sequence[:(-offset)] + ''.join(['.' for i in range(offset)])

    y = ''
    accum = ''
    for i in range(len(a1)): # len(a1) == len(a2)
        base1 = a1.sequence[i]
        base2 = a2.sequence[i]
        conf1 = a1.conf[i]
        conf2 = a2.conf[i]
        
        # Voting round 1: skipping
        if (base1 == '-' or base2 == '-') and \
                max(conf1,conf2) <= quality_threshold:
            accum += 'N'
            continue
        else: # if not skipping, try to insert N's
            if len(accum) > 1:
                y += accum
            accum = ''

        # Voting round 2: base call
        base_key = set()
        if conf1 > quality_threshold:
            base_key.update([base1])
        if conf2 > quality_threshold:
            base_key.update([base2])
        if base1 == base2 and conf1 + conf2 > quality_threshold:
            base_key.update([base1])
        call = iupac_table[as_key(base_key)]
        y += call
    return (y, a1.sequence, a2.sequence)

