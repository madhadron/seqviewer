"""
fasta.py - Module to run and parse FASTA alignment

This module only binds the ssearch36 program, which does pairwise,
global alignment.
"""
import re
import os
import shutil
import tempfile
import subprocess
import contextlib
import Bio.SeqIO

@contextlib.contextmanager
def as_fasta(seq, tmpdir=None, label='sequence'):
    """Create a temporary FASTA files containing the sequence *seq*.

    Many external programs insist on a file, not being fed from stdin.
    as_fasta takes a string, writes it to a temporary file under the
    label *label*, returns that file's name, then deletes the file at
    the end of the with block.
    """
    (db_fd, db_name) = tempfile.mkstemp(text=True, dir=tmpdir)
    db_handle = os.fdopen(db_fd, 'w')
    seqrecord = Bio.SeqIO.SeqRecord(id=label, seq=Bio.SeqIO.Seq(seq))
    Bio.SeqIO.write([seqrecord], db_handle, 'fasta')
    db_handle.close()
    try:
        yield db_name
    finally:
        os.unlink(db_name)

def fasta(seq1, seq2, ssearch36_path="ssearch36", tmpdir='/tmp'):
    with as_fasta(seq1, tmpdir) as fasta1, as_fasta(seq2, tmpdir) as fasta2:
        command = ' '.join([ssearch36_path, '-d', '1', fasta1, fasta2])
        pipe = subprocess.Popen(str(command), shell=True, 
                                stdout=subprocess.PIPE)
        (alignment, _) = pipe.communicate()
        return parse_fasta(alignment, seq1, seq2)

def parse_hunk(h):
    r = re.search(r'^(?:(?:\s+[0-9]+\s*)+\n.{6} ( *)([A-Z-]*)\s*\n)?(?:.+\n)?(?:.{6} ( *)([A-Z-]+)\s*\n(?:\s+[0-9]+\s*)+)', h)

    if r:
        offset1, seq1, offset2, seq2 = r.groups()
        if offset1 != None:
            offset1 = len(offset1)
        if offset2 != None:
            offset2 = len(offset2)
        return (offset1,seq1,offset2,seq2)
    else:
        return None

def parse_fasta(alignment, seg1, seg2):
    lines = alignment.split('\n')
    while True:
        if re.match(r'^>>', lines[0]):
            break
        else:
            lines.pop(0)
    for i in range(4):
        lines.pop(0)
    alseq1, alseq2 = '', ''
    offset1, offset2 = 0, 0
    # Six line hunks: coords, seq1, match, seq2, coords, blank
    both_started = False
    quit_now = False
    while True:
        hunk = '\n'.join(lines[:6])
        lines = lines[6:]
        x = parse_hunk(hunk)
        if x == None:
            break
        maybeoffset1,seq1,maybeoffset2,seq2 = x
        if seq1 and seq2:
            both_started = True
        if maybeoffset1 != None and maybeoffset1 != 0:
            offset1 = maybeoffset1
        if maybeoffset2 != None and maybeoffset2 != 0:
            offset2 = maybeoffset2
        if seq1 != None:
            alseq1 += seq1
            if both_started:
                quit_now = True
        if seq2 != None:
            alseq2 += seq2
            if both_started:
                quit_now = True
        if quit_now:
            break
    # ssearch36 doesn't return the whole sequence if it would mean a
    # line of alignment with only one sequence on it.
    ralseq1 = alseq1.replace('-','')
    if len(ralseq1) < len(seg1):
        i = seg1.find(ralseq1)
        assert i > -1
        alseq1 = seg1[:i] + alseq1 + seg1[i+len(ralseq1):]
        offset1 -= len(ralseq1[:i])
    ralseq2 = alseq2.replace('-','')
    if len(ralseq2) < len(seg2):
        i = seg2.find(ralseq2)
        assert i > -1
        alseq2 = seg2[:i] + alseq2 + seg2[i+len(ralseq2):]
        offset2 -= len(seg2[:i])

    if min(offset1,offset2) < 0:
        offset1 -= min(offset1,offset2)
        offset2 -= min(offset1,offset2)
    assert len([c for c in alseq1 if c != '-']) == len(seg1)
    assert len([c for c in alseq2 if c != '-']) == len(seg2)

    return ((offset1, alseq1), (offset2, alseq2))
        

def test_ssearch():
    s1 = 'CTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCMAGTCGAGCGAACAGATAAGGAGCTTGCTCCTTTGACGTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTACCTATAAGACTGGGACAACTTCGGGAAACCGGAGCTAATACCGGATAATATGTTGAACCGCATGGTTCAATAGTGAAAGATGGTTTTGCTATCACTTATAGATGGACCCGCGCCGTATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATACGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTAGGATCGTAAAACTCTGTTATTAGGGAAGAACAAACGTGTAAGTAACTGTGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACG'
    s2 = 'GATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAACAGATAAGGAGCTTGCTCCTTTGACGTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTACCTATAAGACTGGGACAACTTCGGGAAACCGGAGCTAATACCGGATAATATGTTGAACCGCATGGTTCAATAGTGAAAGATGGTTTTGCTATCACTTATAGATGGACCCGCGCCGTATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATACGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTAGGATCGTAAAACTCTGTTATTAGGGAAGAACAAACGTGTAAGTAACTGTGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTA'
    fasta(s1,s2)

test_ssearch()
