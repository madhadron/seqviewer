import sys
import ab1
import numpy
import contextlib

base_coloring = {'A': 'green', 'C': 'blue', 'T': 'red', 'G': 'black'}

class SequenceTrack(object):
    def __init__(self, sequence):
        self.sequence = sequence
    def render_row(self, limit=None):
        xml = """<div class="track sequence">"""
        for i in range(limit==None and len(self) or limit):
            xml += self.render(i)
        xml += """</div>"""
        return xml
    def render(self, i):
        base = self.sequence[i]
        return """<div class="track-entry %d" style="color: %s">%s</div>""" % \
            (i, base_coloring[base], base)
    def __getitem__(self, i):
        return self.render(i)
    def __str__(self):
        return """SequenceTrack("%s")""" % (self.sequence)
    def __len__(self):
        return len(self.sequence)


class IntegerTrack(object):
    def __init__(self, sequence):
        self.sequence = sequence
    def render_row(self, limit=None):
        xml = """<div class="track integer">"""
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
        self.max_value = float(max(numpy.max(self.trace(b)) for b in 'ACTG'))

        self.centers = numpy.array(centers)

        self.boundaries = numpy.zeros(len(centers)+1)
        self.boundaries[1:-1] = (self.centers[1:] + self.centers[:-1])/2.0
        self.boundaries[-1] = len(A)-1

        # self.Qx = numpy.array([i+1/3.0 for i in range(len(A)-1)])
        # self.Px = numpy.array([i+2/3.0 for i in range(len(A)-1)])

        # self.Py = dict([(b, qs(self.trace(b)))
        #                 for b in 'ACTG'])
        # self.Qy = dict([(b, ps_to_qs(self.Py[b], self.trace(b)))
        #                 for b in 'ACTG'])

    def render_row(self, limit=None):
        xml = """<div class="track chromatogram">"""
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
        m = numpy.sqrt(float(self.max_value))
        for b in 'ACTG':
            i = start
            while i < end:
                Lx, Ly, Rx, Ry = \
                    (i-left)/width, 1-numpy.sqrt(self.trace(b)[i])/m, \
                    (i+1-left)/width, 1-numpy.sqrt(self.trace(b)[i+1])/m
                # This code sparsifies the lines in the SVG: any
                # points that can be approximated adequately (as
                # defined by the function close_enough) by just
                # passing a line on through are omitted.  This makes a
                # huge difference: the HTML output of the example
                # trace file I'm working with drops from 14.4kb to
                # 5.1kb with this code in place.  The rendering time
                # in Firefox goes from 9s to under 7s.

                skipped = []
                while i+1 < end:
                    nextx, nexty = (i+2-left)/width, 1-numpy.sqrt(self.trace(b)[i+2])/m
                    if all([close_enough(Lx,Ly,nextx,nexty,px,py)
                            for px,py in skipped + [(Rx,Ry)]]):
                        skipped += [(Rx,Ry)]
                        Rx, Ry = nextx, nexty
                        i += 1
                    else:
                        break
                xml += """<line x1="%f" y1="%f" x2="%f" y2="%f" stroke-width="0.01"
                              stroke="%s" fill="none" />""" % \
                    (Lx, Ly, Rx, Ry, base_coloring[b])
                i += 1
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


def ab1_to_tracks(handle_or_filename):
    with isinstance(handle_or_filename, basestring) and \
            open(handle_or_filename, 'rb') or \
            liftW(handle_or_filename) as handle:
        a = ab1.Ab1File(handle)
        bases = SequenceTrack(a.bases())
        confidences = IntegerTrack(a.base_confidences())
        chromatogram = ChromatogramTrack(A=a.trace('A'),
                                        C=a.trace('C'),
                                        T=a.trace('T'),
                                        G=a.trace('G'),
                                        centers=a.base_centers())
        return [bases,confidences,chromatogram]



def htmlize(tracks, spacing):
    pass

def fasta_to_track(handle_or_filename):
    pass

