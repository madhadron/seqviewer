import ab1
import numpy
import contextlib

base_coloring = {'A': 'green', 'C': 'blue', 'T': 'red', 'G': 'black'}

class SequenceTrack(object):
    def __init__(self, sequence):
        self.sequence = sequence
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
    def render(self, i):
        val = self.sequence[i]
        return """<div class="track-entry %d">%d</span>""" % (i, val)
    def __getitem__(self, i):
        return self.render(i)
    def __str__(self):
        return """IntegerTrack(%s)""" % (repr(self.sequence))
    def __len__(self):
        return len(self.sequence)


def several(x,n):
    return [x for i in range(n)]

def solve_tridiag(a, b, c, d):
    assert len(a) == len(b)-1
    assert len(c) == len(b)-1
    assert len(b) == len(d)
    N = len(b)
    bprime = numpy.zeros(len(b), dtype=numpy.float)
    bprime[0] = b[0]
    bprime[1:] = b[1:] - a*c/b[:-1].astype(numpy.float)
    dprime = numpy.zeros(len(d), dtype=d.dtype)
    dprime[0] = d[0]
    dprime[1:] = d[1:] - d[:-1]*a/b[:-1].astype(numpy.float)
    xs = numpy.zeros(len(dprime), dtype=d.dtype)
    xs[N-1] = dprime[N-1] / bprime[N-1]
    for i in range(N-2, -1, -1):
        xs[i] = (dprime[i] - xs[i+1]*c[i-1])/bprime[i]
    return xs

class Point(object):
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)
    def __repr__(self):
        return "Point(%s,%s)" % (self.x,self.y)
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)
    def __mul__(self, other):
        return Point(self.x*other, self.y*other)
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    def __div__(self, other):
        return Point(self.x/other, self.y/other)
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

def qs(points):
    N = len(points) - 1
    right_a = numpy.array(several(1.0, N-1))
    right_b = numpy.array([7.0] + several(4.0, N-2) + [2.0])
    right_c = numpy.array([2.0] + several(1.0, N-2))
    right_d = numpy.array([points[0] + points[1]*8] + \
                             [points[i]*2 + points[i+1]*4 for i in range(1,N-1)] + \
                             [points[N-1]*2 + points[N]], dtype=numpy.float)
    assert len(right_a) == N-1
    assert len(right_b) == N
    assert len(right_c) == N-1
    assert len(right_d) == N
    right_control_points = solve_tridiag(right_a, right_b, right_c, right_d)
    return right_control_points

def ps_to_qs(ps, points):
    assert len(ps)+1 == len(points)
    N = len(ps)
    qs = numpy.zeros(N)
    qs[0] = (ps[0] - points[0])/2.0
    qs[1:] = 2*points[1:-1] - ps[:-1]
    return qs



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



        self.Qx = numpy.array([i+1/3.0 for i in range(len(A)-1)])
        self.Px = numpy.array([i+2/3.0 for i in range(len(A)-1)])

        self.Py = dict([(b, qs(self.trace(b)))
                        for b in 'ACTG'])
        self.Qy = dict([(b, ps_to_qs(self.Py[b], self.trace(b)))
                        for b in 'ACTG'])

    def render(self, i):
        left = int(self.boundaries[i])
        right = int(self.boundaries[i+1])
        assert left < self.centers[i]
        assert right > self.centers[i]
        width = float(right-left)
        xml = """<div class="track-entry %d">
                   <svg preserveAspectRatio="none" viewbox="0 0 1 1" version="1.1">""" % i
        start = max(left-1, 0)
        end = min(right, len(self.trace('A'))-1)
        m = numpy.sqrt(float(self.max_value))
        for b in 'ACTG':
            path = "M %f,%f " % (start, self.trace(b)[start]/m)
            for i in range(start, end):
                Lx, Ly, Rx, Ry = \
                    (i-left)/width, 1-numpy.sqrt(self.trace(b)[i])/m, \
                    (i+1-left)/width, 1-numpy.sqrt(self.trace(b)[i+1])/m
                xml += """<line x1="%f" y1="%f" x2="%f" y2="%f" stroke-width="0.01"
                              stroke="%s" fill="none" />""" % \
                    (Lx, Ly, Rx, Ry, base_coloring[b])
        xml += "</svg></div>"
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

