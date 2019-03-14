from setuptools import setup

setup(name='starrpeaker',
      version='0.1',
      description='STARR-seq peak caller',
      keywords = "STARR-seq, peaks, bioinformatics",
      url='http://github.com/hoondy/starrpeaker',
      author='Donghoon Lee',
      author_email='flyingcircus@example.com',
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
      ],
      zip_safe=False
      )