import sys
sys.path.insert(1, '..')

import numpy
import seqviewer

import numpy

def test_render_chromatogram_track():
    print """<html><head>
             <link rel="STYLESHEET" href="viz.css">
             <script src="jquery-1.7.js" type="text/javascript">
             <script type="text/javascript">
             $(document).ready(function() {
                 $("div.track-entry:even").css("background-color", "#ccc");
             });
             </script>
             </head><body>"""
    bases, confs, chrom = seqviewer.ab1_to_tracks('visualization_data/tmpZRPl7_-1.ab1')
    print """<div class="trackset">"""
    print bases.render_row()
    print confs.render_row()
    print chrom.render_row()
    print """</div>"""
    print "</body></html>"

    trace = numpy.array([-0.15, 0.275, 0.25, 0.9, 0.2])
    centers = numpy.array([2])
    track = seqviewer.ChromatogramTrack(trace,trace,trace,trace,centers)
    #print track.render(0)

    trace = numpy.array([1,3,5,12,5,3,4,6,4])
    centers = numpy.array([3,7])
    track = seqviewer.ChromatogramTrack(trace, trace, trace, trace, centers)
    #print track.render(0)
    #print track.render(1)



def test_solve_tridiag():
    a = numpy.array([1,1,1,2])
    b = numpy.array([2,4,4,4,7])
    c = numpy.array([1,1,1,1])
    d = numpy.array([0.4, 1.6, 2.8, 4.0, 7.4])
    xs = seqviewer.solve_tridiag(a,b,c,d)
    assert abs(max(xs - numpy.array([0.066344322, 0.267311333,
                                     0.464410256, 0.658461537,
                                     0.830769231]))) < 1e-7


if __name__=='__main__':
    test_solve_tridiag()
    test_render_chromatogram_track()
