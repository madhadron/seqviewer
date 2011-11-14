import sys
import os
sys.path.insert(1, '..')

import numpy
import seqviewer

import numpy

def test_render_chromatogram_track():
    print """<html><head>
             <link rel="STYLESHEET" href="viz.css">
             </head><body>"""
    print """<div class="trackset">"""
    for f in ['tmpZRPl7_-1.ab1', 'tmpzU_UJw-1.ab1',
              'tmpzUtYnm-1.ab1', 'tmpzV9RJg-1.ab1',
              'tmpzVcdSL-1.ab1']:
        bases, confs, chrom = seqviewer.ab1_to_tracks(os.path.join('traces', f))
        print bases.render_row()
        print confs.render_row()
        print chrom.render_row()
    print """</div>"""
    print "</body></html>"


if __name__=='__main__':
    test_render_chromatogram_track()
