import re
import sys
import numpy
import contextlib

import ab1
import contig

def base_color(base):
    base_coloring = {'A': 'green', 'C': 'blue', 'T': 'red', 'G': 'black'}
    try:
        return base_coloring[base]
    except KeyError:
        return 'yellow'


class SequenceTrack(object):
    def __init__(self, sequence):
        self.sequence = sequence
    def render_row(self, limit=None, offset=0):
        return div(classes=["track","sequence"],
                   [div(classes="track-entry") for i in range(offset)] + \
                       [self.render(i) for i 
                        in range(limit == None and len(self) or min(len(self),limit))])
    def render(self, i):
        base = self.sequence[i]
        return """<div class="track-entry %d" style="color: %s">%s</div>""" % \
            (i, base_color(base), base)
    def __getitem__(self, i):
        return self.render(i)
    def __str__(self):
        return """SequenceTrack("%s")""" % (self.sequence)
    def __len__(self):
        return len(self.sequence)
    def reverse(self):
        return SequenceTrack(self.sequence[::-1])
    def complement(self):
        return SequenceTrack(self.sequence. \
                                 replace('A','t'). \
                                 replace('T','a'). \
                                 replace('C','g'). \
                                 replace('G','c').upper())
    def reverse_complement(self):
        return self.reverse().complement()



class IntegerTrack(object):
    def __init__(self, sequence):
        self.sequence = sequence
    def render_row(self, limit=None, offset=0):
        xml = """<div class="track integer">"""
        for i in range(offset):
            xml += """<div class="track-entry"></div>"""
        for i in range(limit==None and len(self) or limit):
            xml += self.render(i)
        xml += """</div>"""
        return xml
    def render(self, i):
        val = self.sequence[i]
        return """<div class="track-entry %d">%d</div>""" % (i, val)
    def __getitem__(self, i):
        return self.render(i)
    def __str__(self):
        return """IntegerTrack(%s)""" % (repr(self.sequence))
    def __len__(self):
        return len(self.sequence)
    def reverse(self):
        return IntegerTrack(self.sequence[::-1])
    def reverse_complement(self):
        return self.reverse()

def close_enough(Lx,Ly, Rx,Ry, px,py):
    """Is px,py close enough to the line given by L and R to be approximated by it?"""
    # Find the vertical distance of px,py from the line through Lx,Ly
    # and Rx,Ry.  px,py is defined to be "close enough" if it no more
    # than a fraction alpha of the average height of the line away
    # from it.  The value of alpha here was selected by looking at the
    # output by eye and taking the highest value that left the curves
    # still looking reasonably smooth.
    alpha = 0.005
    return abs(py - ((Ry-Ly)/float(Rx-Lx))*(px-Lx) - Ly) < alpha * (Ly + Ry)/2.0

def cutoff(a, n_hinges=6.1):
    m = numpy.median(numpy.log(a+1))
    h = sorted(numpy.log(a+1))[int(0.75*len(a))]
    d = h-m
    c = numpy.exp(m + n_hinges*d) - 1
    return c


class ChromatogramTrack(object):
    def __init__(self, A, C, T, G, centers):
        assert len(A) == len(C)
        assert len(A) == len(T)
        assert len(A) == len(G)
        assert len(centers) < len(A)
        assert centers[-1] < len(A)
        self.A = numpy.array(A).astype(numpy.float)
        self.C = numpy.array(C).astype(numpy.float)
        self.G = numpy.array(G).astype(numpy.float)
        self.T = numpy.array(T).astype(numpy.float)
        all_traces = numpy.concatenate([self.A,self.C,self.T,self.G])
        self.max_value = min(max(all_traces), cutoff(all_traces))

        self.centers = numpy.array(centers)

        self.boundaries = numpy.zeros(len(centers)+1)
        self.boundaries[1:-1] = (self.centers[1:] + self.centers[:-1])/2.0
        self.boundaries[-1] = len(A)-1

    def complement(self):
        return ChromatogramTrack(A=self.T, C=self.G, T=self.A, G=self.C, centers=self.centers)

    def reverse(self):
        return ChromatogramTrack(T=self.T[::-1], G=self.G[::-1], 
                                 A=self.A[::-1], C=self.C[::-1],
                                 centers=len(self.T) - 1 - self.centers[::-1])

    def reverse_complement(self):
        return self.reverse().complement()

    def render_row(self, limit=None, offset=0):
        xml = """<div class="track chromatogram">"""
        for i in range(offset):
            xml += """<div class="track-entry"></div>"""
        for i in range(limit==None and len(self) or limit):
            xml += self.render(i)
        xml += """</div>"""
        return xml

    def render(self, i):
        left = int(self.boundaries[i])
        right = int(self.boundaries[i+1])
        assert left < self.centers[i]
        assert right > self.centers[i]
        width = float(right-left)
        xml = """<div class="track-entry %d"><div class="svg-container">
                   <svg preserveAspectRatio="none" viewbox="0 -0.05 1 1.05" version="1.1">""" % i
        start = max(left-1, 0)
        end = min(right, len(self.trace('A'))-1)
        m = self.max_value
        for b in 'ACTG':
            path = "M%1.2f,%1.2f" % ((start-left)/width,
                                     1-self.trace(b)[start]/m)
            i = start
            while i < end:
                Lx, Ly, Rx, Ry = \
                    (i-left)/width, 1-self.trace(b)[i]/m, \
                    (i+1-left)/width, 1-self.trace(b)[i+1]/m

                # This code sparsifies the lines in the SVG: any
                # points that can be approximated adequately (as
                # defined by the function close_enough) by just
                # passing a line on through are omitted.  This makes a
                # huge difference: the HTML output of the example
                # trace file I'm working with drops from 14.4kb to
                # 5.1kb with this code in place.  The rendering time
                # in Firefox goes from 9s to under 7s.

                # It's also important that I generate a single path
                # per trace instead of a bunch of separate line
                # elements.  This drops the file size another factor
                # of 5, and decreases the rendering time to trivial (~0.7s).

                # I also tried shortening all the div classes and
                # names, but that gave a negligible improvement.

                skipped = []
                # The len(skipped)<10 check speeds the processing up
                # massively, since otherwise on some files we end up
                # checking longer and longer lists over and over.
                while i+1 < end and len(skipped) < 10:
                    nextx, nexty = (i+2-left)/width, 1-self.trace(b)[i+2]/m
                    if all([close_enough(Lx,Ly,nextx,nexty,px,py)
                            for px,py in skipped + [(Rx,Ry)]]):
                        skipped += [(Rx,Ry)]
                        Rx, Ry = nextx, nexty
                        i += 1
                    else:
                        break
                path += """L%1.2f,%1.2f""" % (Rx, Ry)
                i += 1
            xml += """<path stroke-width="0.01" stroke="%s" fill="none" d="%s" />""" % \
                (base_color(b), path)
        xml += "</svg></div></div>"
        return xml
    def trace(self, base):
        if base == 'A':
            return self.A
        elif base == 'C':
            return self.C
        elif base == 'T':
            return self.T
        elif base == 'G':
            return self.G
        else:
            raise ValueError('Invalid base: %s' % base)
    def __getitem__(self, i):
        return self.render(i)
    def __str__(self):
        return "ChromatogramTrack"
    def __len__(self):
        return len(self.centers)


@contextlib.contextmanager
def liftW(x):
    yield x

def _ab1_to_tracks(a, revcomp=False):
    bases = SequenceTrack(a.bases())
    confidences = IntegerTrack(a.base_confidences())
    chromatogram = ChromatogramTrack(A=a.trace('A'),
                                     C=a.trace('C'),
                                     T=a.trace('T'),
                                     G=a.trace('G'),
                                     centers=a.base_centers())
    if revcomp:
        return {'bases': bases.reverse_complement(),
                'confidences': confidences.reverse_complement(),
                'traces': chromatogram.reverse_complement()}
    else:
        return {'bases': bases,
                'confidences': confidences,
                'traces': chromatogram}

def ab1_to_tracks(handle_or_filename):
    with isinstance(handle_or_filename, basestring) and \
            open(handle_or_filename, 'rb') or \
            liftW(handle_or_filename) as handle:
        a = ab1.Ab1File(handle)
        return _ab1_to_tracks(a)


class TrackSet(object):
    def __init__(self):
        self.tracks = {}
        self.offsets = {}
        self.order = []
    def add_track(self, name, track, offset=0, position=None):
        if name in self.tracks.keys():
            raise ValueError("There is already a track named "
                             "%s in this Trackset" % name)
        self.tracks[name] = track
        self.offsets[name] = offset
        if position == None:
            self.order.append(name)
        else:
            if position < len(self.order):
                self.order.insert(position, name)
            else:
                raise ValueError(("Could not insert track %s at "
                                  "position %d in a TrackSet with "
                                  "%d rows.") % (name, position,
                                                 len(self.order)))
    def render(self):
        xml = """<div class="trackset">"""
        xml += """<div class="label-column"><div class="label"><span>Index</span></div>"""
        for n in self.order:
            if isinstance(self.tracks[n], ChromatogramTrack):
                xml += """<div class="chromatogram-label"><span>%s</span></div>""" % n
            else:
                xml += """<div class="label"><span>%s</span></div>""" % n
        xml += """</div>"""
        xml += """<div class="scrolling-container"><div class="trackset-rows">"""
        max_len = max(len(t) for t in self.tracks.itervalues())
        xml += """<div class="track integer">"""
        for i in range(max_len):
            xml += """<div class="track-entry">%d</div>""" % i
        xml += """</div>"""
        for n in self.order:
            xml += self.tracks[n].render_row(offset=self.offsets[n])
        xml += """</div></div>"""
        xml += """</div>"""
        return xml


class ContigDisplay(TrackSet):
    def __init__(self, read1, read2):
        TrackSet.__init__(self)
        read1 = ab1.read_ab1(read1)
        read2 = ab1.read_ab1(read2)

        contig_seq, a1_seq, a2_seq = contig.merge(a1.bases_with_confidence(),
                                                  a2.bases_with_confidence())
            # re.search(r'([^-\.]+)', a1_seq).groups(0)
            

            a1_contig_offset = len(re.split(r'[^\.]', a1_seq, 1)[0])
            a1_contig_end_no_gaps = a1_seq.strip('.').find('-')
            a1_orig_offset = read1_tracks['bases'].sequence.find(a1_seq[a1_contig_offset:a1_contig_end_no_gaps])
            a1_offset = a1_contig_offset - a1_orig_offset

            a2_contig_offset = len(re.split(r'[^\.]', a2_seq, 1)[0])
            a2_contig_end_no_gaps = a2_seq.strip('.').find('-')
            a2_orig_offset = read1_tracks['bases'].sequence.find(a2_seq[a2_contig_offset:a2_contig_end_no_gaps])
            a2_offset = a2_contig_offset - a2_orig_offset

            self.add_track('Read 1 bases', read1_tracks['bases'], offset=a1_offset)
            self.add_track('Read 1 confidences', read1_tracks['confidences'], 
                           offset=a1_offset)
            self.add_track('Read 1 trace', read1_tracks['traces'],
                           offset=a1_offset)
            self.add_track('Read 2 bases', read2_tracks['bases'],
                           offset=a2_offset)
            self.add_track('Read 2 confidences', read2_tracks['confidences'],
                           offset=a2_offset)
            self.add_track('Read 2 trace', read2_tracks['traces'],
                           offset=a2_offset)

            contig_offset = len(re.split(r'[^\.]', contig_seq, 1)[0])
            self.contig_seq = contig_seq.strip('.')
            

            self.add_track('Contig', SequenceTrack(self.contig_seq), 
                           offset=contig_offset, position=0)
        

