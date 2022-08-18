import codecs
import os

from setuptools import setup,find_packages


classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Operating System :: POSIX',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
]

def read(fname):
    return codecs.open(os.path.join(os.path.dirname(__file__), fname)).read()

long_des = read("README.rst")
platforms = ['linux']

install_requires = [
    'biopython==1.78',
    'matplotlib==3.3.3',
    'numpy==1.19.5',
    'pandas==1.2.0',
    'reportlab==3.6.8',
]

setup(name='pdifFiner',
      version='2.1',
      description='This program is designed for annotation of antimicrobal resistance(AMR), pdif site and pdif-ARGs module in bacteria',
      long_description=long_des,
      packages=['PdifFinder'],
      author = "Shao Mengjie",
      author_email = "1437819081@qq.com" ,
      url="https://github.com/mjshao06/PdifFinder.git",
      platforms=platforms,
      classifiers=classifiers,
      install_requires=install_requires,
    )