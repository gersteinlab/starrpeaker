from setuptools import setup

setup(name='starrpeaker',
      version='0.1',
      description='STARR-seq peak caller',
      keywords = "STARR-seq, peaks, bioinformatics",
      url='http://github.com/hoondy/starrpeaker',
      author='Donghoon Lee',
      author_email='donghoon.lee@yale.edu',
      license='GPL',
      packages=['starrpeaker'],
      install_requires=[
          'numpy',
          'scipy',
          'pandas',
          'statsmodels',
          'pysam',
          'pybedtools',
          'pyBigWig'
      ]
      )