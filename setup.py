
"""
TCPG phonon tools: some tools for phonon using using ASE and dftb+ and xtb
"""

from __future__ import absolute_import
from setuptools import setup, find_packages
# from os.path import abspath, dirname
import versioneer
def main():
    """Install the package using setuptools"""
    # project_dir = abspath(dirname(__file__))

    setup(
        name='tcpg_phonon_tools',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        description='Tools for phonon simulation',
        long_description="""
    Just a few handy scripts using ASE dftb+ and xtb.
    """,
        url="https://github.com/ajjackson/mctools",
        author="Donald Chung",
        author_email="donald.chung@stfc.ac.uk",
        license='MIT',

        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            ],
        keywords='chemistry ase dftb+ xTB mantid',
        packages=find_packages(),
        install_requires=['ase', 'spglib', 'matplotlib'],
        entry_points={
            'console_scripts': [
                'ase-convert = mctools.ase_convert:main',
                ]
            },
        include_package_data=True,
        )
        

if __name__ == "__main__":
    main()