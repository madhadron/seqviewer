import sys; sys.path.insert(0, '../')
from seqviewer.assemble import assemble
from seqviewer.tracks import standalone, sequence
import os
import Bio.SeqIO
import cPickle

ts = []
for n in ['tmp_zrpuq', 'tmpzSwR7u', 'tmpzsxOCM',
          'tmpzubQwp', 'tmpzRpKiy', 'tmpzTyEvV', 'tmpzth38k']:
    if os.path.exists('%s.pickle' % n):
        with open('%s.pickle' % n) as h:
            ts,fate = cPickle.load(h)
    else:
        if os.path.exists('%s.fasta' % n):
            extra_seqs = [('lab assembly', 
                           sequence(Bio.SeqIO.read('%s.fasta'%n,
                                                   'fasta').seq.tostring()))]
        else:
            extra_seqs = []
        ts, fate = (n, assemble('%s-1.ab1' % n, '%s-2.ab1' % n, *extra_seqs))
        with open('%s.pickle' % n, 'wb') as h:
            cPickle.dump((ts,fate), h)
    with open('%s.html' % n, 'w') as h:
        print >>h, standalone([(n,ts)])

          
