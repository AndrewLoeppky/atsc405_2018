from setuptools import setup, find_packages

setup(
    name = "a405",
    packages=find_packages(),
    entry_points={
          'console_scripts': [
              'killjobs = a405.utils.killjobs:main',
              'pyncdump = a405.utils.ncdump:main'
          ]
    },

)
