from distutils.core import setup, Extension

module = Extension("mykmeanssp", sources = ["kmeans.c"])

setup(name="mykmeanssp",
      version = 1.0,
      description = "Kmeans algo for swproj 2",
      ext_modules = [module])