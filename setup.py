from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import os

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={self.python_executable}',
            '-DCMAKE_BUILD_TYPE=Release'
        ]

        build_temp = self.build_temp or 'build/temp'
        os.makedirs(build_temp, exist_ok=True)
        os.makedirs(self.build_lib, exist_ok=True)

        os.system(f'cmake {" ".join(cmake_args)} {os.path.abspath(".")}')
        os.system(f'cmake --build . --target pybind_wrapper -- -j4')

setup(
    name='pyrat',
    version='0.1',
    author='Your Name',
    description='Nuclear resonance cross-section calculations',
    packages=['pyrat'],
    ext_modules=[Extension('pyrat.bindings', [])],
    cmdclass={'build_ext': CMakeBuild},
)
