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
        return parse_fasta(alignment)

def offset_and_fill(s):
    (spaces, fill) = re.match(r'^(\s*)([^\s]+)\s*$', s).groups()
    offset = len(spaces)
    return (offset, fill)

def test_offset_and_fill():
    assert offset_and_fill('ACTG') == (0, 'ACTG')
    assert offset_and_fill('    ACTG') == (4, 'ACTG')
    assert offset_and_fill('ACTG    ') == (0, 'ACTG')
    assert offset_and_fill('    ACTG   ') == (4, 'ACTG')

def parse_fasta(alignment):
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
    while lines[0] != '':
        hunk = lines[:6]
        lines = lines[6:]
        maybeoffset1, seq1 = offset_and_fill(hunk[1][7:])
        maybeoffset2, seq2 = offset_and_fill(hunk[3][7:])
        if maybeoffset1 != 0:
            offset1 = maybeoffset1
        if maybeoffset2 != 0:
            offset2 = maybeoffset2
        alseq1 += seq1
        alseq2 += seq2
    return ((offset1, alseq1), (offset2, alseq2))
        
