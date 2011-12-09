import os
import shutil
import tempfile

basedir = '/Volumes/clinical_workups'
targetdir = 'scraped_files'

def scrape(path):
    ab1_files = [x for x in os.listdir(path) if x.endswith('.ab1')]
    if ab1_files != []:
        print '%s:' % path
        filename = tempfile.mktemp(prefix='', dir=targetdir)
        for i,f in enumerate(ab1_files):
            src_file = os.path.join(path, f)
            tgt_file = '%s-%d.ab1' % (filename,i)
            print '  %s --> %s' % (src_file, tgt_file)
            shutil.copyfile(src_file, tgt_file)
        fasta_files = sorted([x for x in os.listdir(path) if x.endswith('.fasta')],
                             key=len)
        if fasta_files != []:
            src_file = os.path.join(path, fasta_files[0])
            tgt_file = '%s.fasta' % filename
            print '  %s --> %s' % (src_file, tgt_file)
            shutil.copyfile(src_file, tgt_file)
        
    else:
        for d in os.listdir(path):
            fullpath = os.path.join(path, d)
            if os.path.isdir(fullpath):
                scrape(fullpath)
            
scrape(basedir)
