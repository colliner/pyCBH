from setuptools import setup, find_packages

setup(
    name = 'pyCBHt',
    version = '0.0.1',
    url = '',
    description = '',
    packages = find_packages(),
    install_requires = [
        # Githubrepo
        'xyz2mol @ git+ssh://git@github.com/colliner/xyz2mol.git'
    ]
)

