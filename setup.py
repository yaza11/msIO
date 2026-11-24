import os
import re

from setuptools import setup, find_packages

# the version is kept in msIO/_version.py (also shown in the GUI), read without importing msIO and its dependencies
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'msIO', '_version.py')) as f:
    version = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", f.read()).group(1)

setup(
    name='msIO',
    version=version,
    url='',
    author='Yannick Zander',
    description='Reading and writing functions for different mass spectrometry '
                'related file types',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'pandas', 'tqdm', 'sqlalchemy', 'networkx', 'matchms', 'sqlalchemy'],
)

# pip install git+https://github.com/yaza11/LipidCalculator.git
