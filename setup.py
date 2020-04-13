from setuptools import setup

setup(name='starrpeaker',
      version='1.0',
      description='STARRPeaker: STARR-seq peak caller',
      keywords = "STARRPeaker, STARR-seq, peaks, bioinformatics",
      url='http://github.com/gersteinlab/starrpeaker',
      author='Donghoon Lee',
      author_email='donghoon.lee@yale.edu',
      license='GPL',
      packages=['starrpeaker'],
      install_requires=[
            'numpy',
            'scipy',
            'pandas',
            'statsmodels<=0.10.2',
            'pysam',
            'pybedtools',
            'pyBigWig',
            'scikit-learn'
      ],
      entry_points = {
          'console_scripts': [
              'starrpeaker = starrpeaker.starrpeaker:main'
          ]
          }
      )
