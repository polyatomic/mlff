from setuptools import setup, Extension

sfc_sources = [
'module.cpp',
'../chemlib/Atom.cpp',
'../chemlib/chemutil.cpp',
'../chemlib/Descriptors.cpp',
'../chemlib/Mol.cpp',
'../chemlib/parsers.cpp',
'../mathlib/mathutil.cpp',
'../util/ThreadPool.cpp',
'../util/util.cpp'
]

sfc_include_dirs = [
'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include',
'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/fftw',
'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/include/intel64',
'../chemlib',
'../util',
'../mathlib'
]

sfc_libraries = [
'mkl_intel_lp64',
'mkl_sequential',
'mkl_core'
]

sfc_library_dirs = [
'C:/Program Files (x86)/Intel/oneAPI/mkl/2021.3.0/lib/intel64'
]

sfc_module = Extension('chemutil', sources = sfc_sources, include_dirs = sfc_include_dirs, libraries = sfc_libraries, library_dirs = sfc_library_dirs)

setup(
   name='chemutil',
   version='1.0',
   description='Python Package with chemutil C++ extension',
   ext_modules=[sfc_module]
)
