import os
import re

from setuptools import setup, find_packages

# the version is kept in msIO/_version.py (also shown in the GUI), read without importing msIO and its dependencies
HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, 'msIO', '_version.py')) as f:
    version = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", f.read()).group(1)

# the README is the project page on PyPI
with open(os.path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='msIO',
    version=version,
    url='https://github.com/yaza11/msIO',
    author='Yannick Zander',
    description='Reading and writing functions for different mass spectrometry '
                'related file types',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    python_requires='>=3.10',
    # only msIO, not the local folders next to it (testing, ...)
    packages=find_packages(include=['msIO', 'msIO.*']),
    # non-Python files used at runtime: the guide shown in the GUI, the SQL schema
    package_data={'msIO.gui': ['about'], 'msIO.sql': ['*.sql']},
    install_requires=['numpy', 'matplotlib', 'pandas', 'tqdm', 'sqlalchemy', 'networkx', 'matchms', 'rdkit'],
)

# pip install git+https://github.com/yaza11/LipidCalculator.git
