from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import subprocess

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',  # FIX: sys.executable not self.python_executable
            '-DCMAKE_BUILD_TYPE=Release'
        ]

        build_temp = self.build_temp
        os.makedirs(build_temp, exist_ok=True)

        # Configure step
        subprocess.check_call(['cmake', os.path.abspath('.')] + cmake_args, cwd=build_temp)

        # Build step
        subprocess.check_call(['cmake', '--build', '.', '--target', 'pybind_wrapper', '--', '-j4'], cwd=build_temp)

setup(
    name='pyrat',
    version='0.1',
    author='Your Name',
    description='Nuclear resonance cross-section calculations',
    packages=find_packages(include=['scripts', 'scripts.*']),  # FIX: if you want to expose scripts
    ext_modules=[Extension('pyrat.bindings', sources=[])],
    cmdclass={'build_ext': CMakeBuild},
    zip_safe=False,
)
