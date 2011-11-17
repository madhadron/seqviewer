from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("seqviewer.ab1", ["seqviewer/ab1.pyx"])],
    include_dirs = [numpy.get_include(),],
)
