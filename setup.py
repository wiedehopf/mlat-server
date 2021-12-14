from setuptools import Extension, setup
from Cython.Build import cythonize

import os

# remove potential bogus so files
pyFiles = [s[:-3] for s in os.listdir('mlat') if s[-2:] == 'py']
removeFiles = [s for s in os.listdir('mlat') if s[-2:] != 'py' and s.split('.')[0] in pyFiles]
for f in removeFiles:
    os.remove('mlat/' + f)

extensions = [
    Extension("modes_cython.message", ["modes_cython/message.pyx"]),
    Extension("mlat.geodesy", ["mlat/geodesy.pyx"]),
    Extension("mlat.clocknorm", ["mlat/clocknorm.pyx"]),
    Extension("mlat.clocktrack", ["mlat/clocktrack.pyx"])
]
setup(
        ext_modules = cythonize(
            extensions,
            compiler_directives={'language_level' : "3"}
            )
        )
