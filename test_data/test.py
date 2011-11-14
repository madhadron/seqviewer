import sys
import os
sys.path.insert(1, '..')

import numpy
import seqviewer

import numpy

def test_render_chromatogram_track():
    t = seqviewer.TrackSet()
    for f in ['tmpzU_UJw-1.ab1']:
        # ['tmpZRPl7_-1.ab1']: #, 'tmpzU_UJw-1.ab1']:
        # 'tmpzUtYnm-1.ab1', 'tmpzV9RJg-1.ab1',
        # 'tmpzVcdSL-1.ab1']:
        bases, confs, chrom = seqviewer.ab1_to_tracks(os.path.join('traces', f))
        t.add_track(bases, name=f+'-bases')
        t.add_track(confs, name=f+'-confs')
        t.add_track(chrom, name=f+'-chrom')

    print """<html><head>
             <link rel="STYLESHEET" href="viz.css">
             </head><body>"""
    print t.render()
    print "</body></html>"


if __name__=='__main__':
    test_render_chromatogram_track()
