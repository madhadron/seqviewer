import tracks
import contig
import ab1
import fasta

def assemble(read1, read2, *extra_seqs):
    tracks1 = ab1.read(read1)
    tracks2 = ab1.read(read2)
    ref = contig.contig(tracks1['sequence'], tracks1['confidences'],
                        tracks.revcomp(tracks2['sequence']), 
                        tracks.revcomp(tracks2['confidences']),
                        threshold=20)
    t = tracks.TrackSet()

    read1_offset, read1_sequence = ref['read1']
    read2_offset, read2_sequence = ref['read2']
    read1_confs = tracks.regap(read1_sequence, tracks1['confidences'])
    read2_confs = tracks.regap(read2_sequence, tracks.revcomp(tracks2['confidences']))
    read1_traces = tracks.regap(read1_sequence, tracks1['traces'])
    read2_traces = tracks.regap(read2_sequence, tracks.revcomp(tracks2['traces']))

    t.extend([
              tracks.TrackEntry('read 1 traces', read1_offset, read1_traces),
              tracks.TrackEntry('read 1 confidences', read1_offset, read1_confs),
              tracks.TrackEntry('read 1 bases', read1_offset, read1_sequence),
              tracks.TrackEntry('read 2 traces', read2_offset, read2_traces),
              tracks.TrackEntry('read 2 confidences', read2_offset, read2_confs),
              tracks.TrackEntry('read 2 bases', read2_offset, read2_sequence)])


    if ref['reference'] != None:
        reference_offset, reference_sequence = ref['reference']
        t.append(tracks.TrackEntry('reference', reference_offset, reference_sequence))

    for (name,s) in extra_seqs:
        if ref['reference'] != None:
            (roffset, _), (soffset, saligned) = fasta.fasta(reference_sequence, s)
            t.append(tracks.TrackEntry(name, reference_offset + soffset - roffset,
                                       tracks.sequence(saligned)))
        else:
            t.append(tracks.TrackEntry(name, 0, s))
    return t

# def test_assemble():
#     import Bio.SeqIO
#     s = Bio.SeqIO.read('../test_data/traces/tmpZRPl7_.fasta', 'fasta').seq.tostring()
#     t = assemble('../test_data/traces/tmpZRPl7_-1.ab1',
#                  '../test_data/traces/tmpZRPl7_-2.ab1',
#                  s)
#     with open('test.html', 'w') as h:
#         print >>h, tracks.standalone(t)
        
