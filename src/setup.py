from distutils.core import setup, Extension

module = Extension("myspkmeanssp", sources = ["spkmeansmodule.c", "spkmeans.c"])

setup(name="myspkmeanssp",
      version = 1.0,
      description = "spkmeans functionality for the final project",
      ext_modules = [module])