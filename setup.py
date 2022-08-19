import codecs
import os

from setuptools import setup,find_packages

def expand_data(data_to_expand):
    """
    From data structure like setup.py data_files (see http://)
     [(directory/where/to/copy/the/file, [path/to/file/in/archive/file1, ...]), ...]
    but instead of the original struct this one accept to specify a directory in elements to copy.
    This function will generate one entry for all *content* of the directory and subdirectory
    recursively, to in fine copy the tree in archive in dest on the host
    the first level of directory itself is not include (which allow to rename it)
    :param data_to_expand:
    :type  data_to_expand: list of tuple
    :return: list of tuple
    """
    def remove_prefix(prefix, path):
        prefix = os.path.normpath(prefix)
        path = os.path.normpath(path)
        to_remove = len([i for i in prefix.split(os.path.sep) if i])
        truncated = [i for i in path.split(os.path.sep) if i][to_remove:]
        truncated = os.path.sep.join(truncated)
        return truncated

    data_struct = []
    for base_dest_dir, src in data_to_expand:
        base_dest_dir = os.path.normpath(base_dest_dir)
        for one_src in src:
            if os.path.isdir(one_src):
                for path, _, files in os.walk(one_src):
                    if not files:
                        continue
                    path_2_create = remove_prefix(one_src, path)
                    data_struct.append((os.path.join(base_dest_dir, path_2_create), [os.path.join(path, f) for f in files]))
            if os.path.isfile(one_src):
                data_struct.append((base_dest_dir, [one_src]))
    return data_struct

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
]

setup(name='PdifFinder',
      version='1.1.1',
      description='This program is designed for annotation of antimicrobal resistance(AMR), pdif site and pdif-ARGs module in bacteria',
      long_description=long_des,
      packages=['PdifFinder'],
      author = "Shao Mengjie",
      author_email = "1437819081@qq.com" ,
      url="https://github.com/mjshao06/pdifFinder.git",
      platforms=platforms,
      classifiers=classifiers,
      install_requires=install_requires,
      data_files=expand_data([('PdifFinder/AMRDB',['PdifFinder/AMRDB']),('PdifFinder/data',['PdifFinder/data'])]),
      entry_points = {
          'console_scripts': [
              'pdifFinder = PdifFinder.pdifFinder:main',
          ]
      }
    )
