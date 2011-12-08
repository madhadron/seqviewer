import sys; sys.path.insert(0, '../')
import os
import cPickle
import seqviewer.tracks
import collections

AssemblyStats = collections.namedtuple('AssemblyStats',
 ['ref_fate','lab_fate','percent_id','frac_overlap','lab_len_minus_ref_len'])

lab_fate = []
alg_fate = []
percent_ids = []
frac_overlap = []
lab_minus_alg_len = []
names = []

for n in [x[:-7] for x in os.listdir('.') if x.endswith('.pickle')]:
    with open('%s.pickle' % n) as h:
        ts,fate = cPickle.load(h)
    ref = filter(lambda x: x.name == 'reference', ts)
    lab = filter(lambda x: x.name == 'lab assembly', ts)

    names.append(n)
    lab_fate.append(lab != [] and 'lab assembled' or 'lab unassembled')
    alg_fate.append(fate)

    if lab != [] and fate != 'none':
        lab = lab[0]
        ref = ref[0]
        diff = min(ref.offset+len(ref.track), lab.offset+len(lab.track)) - \
            max(ref.offset, lab.offset)
        frac_overlap.append(diff / (0.5 * (len(ref.track) + len(lab.track))))
        lab_minus_alg_len.append(len(lab.track) - len(ref.track))

        refseg = ref.track[lab.offset - min(ref.offset,lab.offset):]
        labseg = lab.track[ref.offset - min(ref.offset,lab.offset):]
        numer = sum([1 for x,y in zip(refseg,labseg) if x==y])
        denom = sum([1 for x,y in zip(refseg,labseg)])
        percent_ids.append(float(numer)/float(denom))
    else:
        percent_ids.append(None)
        frac_overlap.append(None)
        lab_minus_alg_len.append(None)
    if (percent_ids[-1] != None and percent_ids[-1] < 0.99) or \
            (frac_overlap[-1] != None and frac_overlap[-1] < 0.90):
        with open('%s.html' % n, 'w') as h:
            print >>h, seqviewer.tracks.standalone([(n,ts)])

print 'd <- data.frame('
print '    names=c(' + ','.join(['"'+x+'"' for x in names]) + '),'
print '    lab_fate=c(' + ','.join(['"'+x+'"' for x in lab_fate]) + '),'
print '    alg_fate=c(' + ','.join(['"'+x+'"' for x in alg_fate]) + '),'
print '    percent_ids=c(' + ','.join([str(x) for x in percent_ids]) + '),'
print '    frac_overlap=c(' + ','.join([str(x) for x in frac_overlap]) + '),'
print '    lab_minus_alg_len=c(' + ','.join([str(x) for x in lab_minus_alg_len]) + '))'
    

