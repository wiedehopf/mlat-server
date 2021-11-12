from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension("mlat.clocksync", ["mlat/clocksync.pyx"]),
    Extension("mlat.clocknorm", ["mlat/clocknorm.pyx"]),
    Extension("mlat.clocktrack", ["mlat/clocktrack.pyx"])
]
setup(
        ext_modules = cythonize(
            extensions,
            compiler_directives={'language_level' : "3"}
            )
        )
