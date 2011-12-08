import sys; sys.path.insert(0, '../')
from seqviewer.assemble import assemble
from seqviewer.tracks import standalone, sequence
import os
import Bio.SeqIO
import cPickle

for n in [x[:-6] for x in os.listdir('.') if x.endswith('-1.ab1')]:
    if not(os.path.exists('%s.pickle' % n)):
        if os.path.exists('%s.fasta' % n):
            extra_seqs = [('lab assembly', 
                           sequence(Bio.SeqIO.read('%s.fasta'%n,
                                                   'fasta').seq.tostring()))]
        else:
            extra_seqs = []
        t, fate = assemble('%s-1.ab1' % n, '%s-2.ab1' % n, *extra_seqs)
        with open('%s.pickle' % n, 'wb') as h:
            cPickle.dump((t,fate), h)
          
