from setuptools import setup, find_packages

setup(
    name='msIO',
    version='0.0.1',
    url='',
    author='Yannick Zander',
    author_email='yannick.zander@gmail.com',
    description='Reading and writing functions for different mass spectrometry '
                'related file types',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'pandas', 'tqdm'],
)
