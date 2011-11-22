
import numpy
from collections import namedtuple

def base_color(base):
    base_coloring = {'A': 'green', 'C': 'blue', 'T': 'red', 
                     'G': 'black', 'U': 'red'}
    try:
        return base_coloring[base]
    except KeyError:
        return 'yellow'

class Traces(list):
    def __comp__(self):
        return Traces([{'A': d['T'], 'C': d['G'], 'T': d['A'], 'G': d['C']}
                       for d in self])
    def __rev__(self):
        return [{'A': rev_entry(d['A']), 
                 'T': rev_entry(d['T']),
                 'C': rev_entry(d['C']),
                 'G': rev_entry(d['G'])}
                for d in self[::-1]]
    def __render__(self, offset=0):
        entries = [div(classes=['track-entry','global'+str(i)], body='')
                   for i in range(offset)]
        for i in range(len(self)):
            entry = self[i]
            paths = ''
            for b in 'ACTG':
                paths += path(M(entry[b][0]) + \
                                  ''.join([L(x) for x in entry[b][1:]]),
                              stroke=base_color(b))
            entries.append(div(classes=['track-entry',
                                        'global'+str(i+offset),
                                        'local'+str(i)],
                               body=div(classes='svg-container',
                                        body=unit_svg(paths))))
        return div(classes=['track','chromatogram'],
                   body=''.join(entries))
    gap = {'A': numpy.array([(0.5,0)]), 'C': numpy.array([(0.5,0)]),
           'T': numpy.array([(0.5,0)]), 'G': numpy.array([(0.5,0)])}


def rev_entry(d):
    return [(1-x,y) for x,y in d[::-1]]


def cutoff(a, n_hinges=6.1):
    m = numpy.median(numpy.log(a+1))
    h = sorted(numpy.log(a+1))[int(0.75*len(a))]
    d = h-m
    c = numpy.exp(m + n_hinges*d) - 1
    return c

def sparsify(xs, ys):
    new_points = [(xs[0], ys[0])]
    i = 0
    while i < len(ys)-1:
        Lx, Ly, Rx, Ry = xs[i], ys[i], xs[i+1], ys[i+1]
        skipped = []
        while i < len(ys)-2 and len(skipped) < 10:
            nextx, nexty = xs[i+2], ys[i+2]
            if all([close_enough(Lx,Ly, nextx,nexty, px,py)
                    for px,py in skipped + [(Rx,Ry)]]):
                skipped += [(Rx,Ry)]
                Rx,Ry = nextx,nexty
                i += 1
            else:
                break
        new_points.append((Rx,Ry))
        i += 1
    return new_points

def test_sparsify():
    xs = range(14)
    ys = numpy.arange(14)/5.0
    assert sparsify(xs,ys) == [(0,0), (11,11/5.0), (13,13/5.0)]
    xs = range(4)
    ys = numpy.arange(4)
    assert sparsify(xs,ys) == [(0,0), (3,3)]


def traces(A, C, T, G, centers):
    A = numpy.array(A).astype(numpy.float)
    C = numpy.array(C).astype(numpy.float)
    G = numpy.array(G).astype(numpy.float)
    T = numpy.array(T).astype(numpy.float)
    centers = numpy.array(centers).astype(numpy.integer)
    assert len(A) == len(C) == len(T) == len(G)
    assert all(centers >= 0) and all(centers < len(A))
    assert all(sorted(centers) == centers)
    _limits = [int(numpy.ceil((centers[i]+centers[i-1])/2.0)) 
               for i in range(1, len(centers))]
    limits = zip([0] + _limits, [x+1 for x in _limits] + [len(A)])
    print limits
    all_traces = numpy.concatenate([A,C,T,G])
    m = min(max(all_traces), cutoff(all_traces))
    t = Traces()
    for l,r in limits:
        xs = numpy.arange(0,r-l) / float(r-l-1)
        assert len(xs) == len(A[l:r]) == len(C[l:r]) \
            == len(T[l:r]) == len(G[l:r])
        t.append({'A': sparsify(xs, 1-A[l:r]/m), 
                  'C': sparsify(xs, 1-C[l:r]/m),
                  'T': sparsify(xs, 1-T[l:r]/m),
                  'G': sparsify(xs, 1-G[l:r]/m)}), 
    return t

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

    
def test_traces():
    centers = [2,4,7]
    ys = range(14)
    t = traces(ys,ys,ys,ys,centers)
    for x in t:
        assert sorted(list(x.keys())) == ['A','C','G','T']
    assert len(t[0]['A']) == 2
    assert len(t[0]['C']) == 2
    assert len(t[1]['A']) == 2



class Sequence(str):
    def __comp__(self):
        return Sequence(self. \
                            replace('A','t'). \
                            replace('T','a'). \
                            replace('C','g'). \
                            replace('G','c').upper())
    def __rev__(self):
        return Sequence(self[::-1])
    def __render__(self, offset=0):
        entries = \
            [div(classes=['track-entry','global'+str(i)], body='')
             for i in range(offset)] + \
            [div(classes=['track-entry','global'+str(i+offset),'local'+str(i)],
                 style="color: %s" % base_color(self[i]),
                 body=self[i])
             for i in range(len(self))]
        return div(classes=['track','sequence'],
                   body=''.join(entries))
    gap = '-'
    def __isgap__(self, other):
        return other == '-' or other == '.'

def sequence(s):
    return Sequence(s)

class Numeric(numpy.ndarray):
    def __rev__(self):
        return self[::-1]
    def __comp__(self):
        return self
    def __render__(self, offset=0):
        entries = \
            [div(classes=['track-entry', 'global'+str(i)],
                 body='') for i in range(offset)] + \
            [div(classes=['track-entry','global'+str(i+offset),'local'+str(i)],
                 body=str(self[i]))
             for i in range(len(self))]
        return div(classes=['track','integer'],
                   body=''.join(entries))
    gap = numpy.NaN
    def __isgap__(self, other):
        numpy.isnan(other)

def numeric(vals):
    a = Numeric(len(vals))
    a[:] = numpy.array(vals)
    return a

TrackEntry = namedtuple('TrackEntry', ['name', 'offset', 'track'])

class TrackSet(list):
    def __rev__(self):
        new = TrackSet()
        for t in self:
            pass # FIXME
    def __comp__(self):
        new = TrackSet()
        for t in self:
            new.append(TrackEntry(name=t.name, offset=t.offset,
                                  track=comp(t.track)))
        return new
    def __render__(self):
        labels = div(classes='label', body=span('Index'))
        for t in self:
            if isinstance(t.track,Traces):
                labels += div(classes='chromatogram-label',
                              body=span(t.name))
            else:
                labels += div(classes='label',
                              body=span(t.name))
        maxlen = max([len(x.track)+x.offset for x in self])
        tracks = div(classes=['track','integer'],
                     body=''.join([div(classes=['track-entry'], body=str(i))
                                   for i in range(maxlen)]))
        for t in self:
            tracks += render(t.track, offset=t.offset)

        return div(classes='trackset',
                   body=div(classes='label-column', body=labels) + \
                       div(classes='scrolling-container',
                           body=div(classes='trackset-rows',
                                    body=tracks)))


def rev(s):
    return s.__rev__()

def comp(s):
    return s.__comp__()

def revcomp(s):
    return rev(comp(s))

def render(s, *args, **kwargs):
    return s.__render__(*args, **kwargs)


def regap(template, target):
    gaps = [i for i,t in enumerate(template)
            if t == t.gap]
    


### XML
def tag(name):
    def f(body="", classes=[], style=""):
        if not(isinstance(classes, list)):
            classes = [classes]
        return ("""<%s class="%s" style="%s">""" % (name,
                                                    ' '.join(classes),
                                                    style)) + \
            body + \
            ("""</%s>""" % (name,))
    return f

div = tag('div')
span = tag('span')

def unit_svg(body=""):
    return """<svg preserveAspectRatio="none" viewbox="0 -0.05 1 1.05" version="1.1">""" + \
        body + """</svg>"""

def M((x,y)):
    return "M%0.3f,%0.3f" % (x,y)

def L((x,y)):
    return "L%0.3f,%0.3f" % (x,y)

def path(d, stroke="black", strokeWidth="0.01", fill="none"):
    dstr = ''.join(d)
    return """<path stroke="%s" stroke-width="%s" fill="%s" d="%s" />""" % \
        (stroke, strokeWidth, fill, dstr)
        
def standalone(tracksets):
    if isinstance(tracksets, TrackSet):
        tracksets = [tracksets]
    xml = """<html><head>\n"""
    xml += """<style>\n""" + stylesheet + "</style>\n"
    xml += """</head><body>"""
    for t in tracksets:
        xml += render(t)
    xml += """</body></html>"""
    return xml


stylesheet = """
* {
    margin: 0;
    padding: 0;
}

.scrolling-container {
    overflow: scroll;
    white-space: nowrap;
}

.trackset-rows {
    display: table;
}


.label-column {
    float: left;
    max-width: 10em;
    overflow: hidden;
    white-space: nowrap;
    border-right: 0.2em solid black;
}


.label-column div {
    display: block;
    font-family: Optima, Myriad, sans-serif;
    vertical-align: middle;
    color: white;
    border-top: 1px solid #eee;
    background-color: #111;
    background-image: url('concrete_wall.png');
    padding-right: 0.1em;
    padding-left: 0.2em;
    padding-top: 0.2em;
    padding-bottom: 0.2em;
    height: 0.6em;
    text-align: right;
}

.label-column div span {
    font-size: 0.6em;
}

.label-column div:first-child {
    border-top: none;
}

.chromatogram-label {
    height: 3.6em !important;
}

.chromatogram-label span {
    height: 6em !important;
    line-height: 6em;
}

.trackset {
    width: 100%;
    max-width: 100%;
}

div.track {
    display: table-row;
}

.track-entry {
    padding-left: 0.3em;
    padding-right: 0.3em;
    display: table-cell;
    text-align: center;
}


.track .track-entry:nth-of-type(odd) {
    background-color: #eee;
}

.chromatogram > .track-entry {
    height: 4em;
    padding: 0;
}

.track-entry > .svg-container {
    width: 100%;
    height: 100%;
    font-size: 0;
}

.integer {
    font-size: 70%;
    color: #666;
    line-height: 143%;
}

.track-name {
    position: fixed;
    vertical-align: middle;
    left: 0;
    text-align: right;
    white-space: wrap;
    overflow: wrap;
}
"""
