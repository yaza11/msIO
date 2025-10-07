from setuptools import setup, find_packages

setup(
    name='msIO',
    version='0.0.14',
    url='',
    author='Yannick Zander',
    author_email='yannick.zander@gmail.com',
    description='Reading and writing functions for different mass spectrometry '
                'related file types',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'pandas', 'tqdm', 'sqlalchemy', 'networkx'],
)

# pip install git+https://github.com/yaza11/LipidCalculator.git
