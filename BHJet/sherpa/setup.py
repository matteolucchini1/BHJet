"""Build the Sherpa model to the BHJet model"""

import glob
import os
import setuptools
import subprocess as sbp

# Needs to be updated to support Python 3.10
#
from distutils.core import Extension
from distutils.command.build_ext import build_ext

import numpy
from numpy.distutils import fcompiler

import sherpa

def get(tool, args):
    ans = sbp.run([tool] + args, stdout=sbp.PIPE)
    return ans.stdout.decode('ascii').strip()


def startswith(token, start):
    if token.startswith(start):
        return

    raise OSError(f"Expected {start}.. but sent '{token}'")


gsl_cflags = get("gsl-config", ["--cflags"])
gsl_iflags = get("gsl-config", ["--libs"])

# decode the gsl output; seems a bit excessive
#
if len(gsl_cflags.split()) != 1:
    raise OSError(f"Expected no spaces in '{gsl_cflags}'")

startswith(gsl_cflags, "-I")

gsl_incpath = gsl_cflags[2:]

gsl_iflags = gsl_iflags.split()

startswith(gsl_iflags[0], "-L")
for token in gsl_iflags[1:]:
    startswith(token, "-l")

gsl_libpath = gsl_iflags[0][2:]
gsl_libnames = [token[2:] for token in gsl_iflags[1:]]

# Build the interface to the models and the model code.
#
# Assume the code in Kariba/ does not need to be built in a
# specific order.
#
search = "../../Kariba/*.cpp"
karibas = glob.glob(search)
if len(karibas) == 0:
    raise OSError(f"No match for: {search}")

bhs = ["../bhjet.cpp", "../utils.cpp", "../jetpars.cpp", "../wrap_xspec.cpp"]

srcnames = ["src/bhjet/_models.cxx"] + bhs + karibas

this_incpath = os.getcwd() + '../'
numpy_incpath = numpy.get_include()
sherpa_incpath = sherpa.get_include()

includes = [this_incpath, numpy_incpath, sherpa_incpath, gsl_incpath]

libs = [gsl_libpath]
libnames = gsl_libnames

# Do we need this?
if os.uname().sysname == 'Darwin':
    cargs = ['-Wl,-no_compact_unwind']
else:
    cargs = []

mod = Extension('bhjet._models',
                include_dirs=includes,
                library_dirs=libs,
                libraries=libnames,
                sources=srcnames,
                extra_link_args=cargs,
                depends=[]
                )

setuptools.setup(ext_modules=[mod])
