from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
import setuptools
import pybind11

__version__ = '0.0.1'

# Get the root directory of the project
root_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Try to import numpy and pybind11, but don't fail immediately if they're not available
try:
    import numpy as np
    numpy_includes = [np.get_include()]
except ImportError:
    numpy_includes = []

# Don't try to import pybind11 during the initial setup
include_dirs = numpy_includes + [
    os.path.join(root_dir, "bonsai/hll/include"),  # Path to sketch headers
    os.path.join(root_dir, "src"),                 # Path to dashing2 headers
    os.path.join(root_dir, "bonsai"),             # Path to bonsai headers
]

# Adapted from https://github.com/pybind/python_example
class get_pybind_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        return pybind11.get_include(self.user)

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': [],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.9']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        # Ensure numpy and pybind11 are available during the actual build
        import numpy as np
        import pybind11
        for ext in self.extensions:
            ext.include_dirs.extend([
                np.get_include(),
                pybind11.get_include(),
                pybind11.get_include(user=True)
            ])
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])

        if ct == 'unix':
            opts.append('-std=c++17')
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        super().build_extensions()

setup(
    name='dashing2',
    version=__version__,
    author='Jessica Bonnie',
    author_email='your.email@example.com',
    description='Python bindings for Dashing2 - fast sequence sketching',
    long_description='''
    Dashing2 is a fast and memory-efficient library for computing similarity
    between biological sequences using MinHash sketching.
    
    Features:
    - Fast sketching of FASTA/FASTQ files
    - Memory-efficient MinHash implementation
    - Computation of Jaccard similarities and distances
    - Support for k-mer based sequence comparison
    ''',
    long_description_content_type='text/markdown',
    ext_modules=[
        Extension(
            "dashing2.core",
            ["dashing2_bindings.cpp"],
            include_dirs=include_dirs,
            libraries=['z'],
            extra_compile_args=['-std=c++17'],
            language='c++'
        )
    ],
    packages=['dashing2'],
    package_data={
        'dashing2': ['*.so', 'core.*.so'],
    },
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    python_requires=">=3.6",
    setup_requires=[
        'numpy>=1.15.0',
        'pybind11>=2.6.0',
    ],
    install_requires=[
        'numpy>=1.15.0',
        'pybind11>=2.6.0',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
