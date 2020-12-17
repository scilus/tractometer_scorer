import os

from Cython.Build import cythonize
import glob
import numpy
from setuptools import setup, find_packages
from setuptools.extension import Extension

opts = dict(name='ismrm_2015_tractography_challenge_scoring',
            maintainer='Jean-Christophe Houde',
            maintainer_email='jean-christophe.houde@usherbrooke.ca',
            description='Scoring system used for the ISMRM 2015 Tractography Challenge',
            url='https://github.com/scilus/ismrm_2015_tractography_challenge_scoring',
            author='The challenge team',
            version='1.0.1',
            packages=find_packages(),
            install_requires=['dipy', 'nibabel'],
            requires=['numpy', 'scipy'],
            scripts=glob.glob('scripts/*.py'))


extensions = [Extension(
    'challenge_scoring.tractanalysis.robust_streamlines_metrics',
    ['challenge_scoring/tractanalysis/robust_streamlines_metrics.pyx'],
    include_dirs=[numpy.get_include(),
                  os.path.join(
                  os.path.dirname(
                      os.path.realpath(__file__)),
                      'challenge_scoring/c_src')])]


opts['ext_modules'] = cythonize(extensions, compiler_directives={'language_level': "3"})


if __name__ == '__main__':
    setup(**opts)
