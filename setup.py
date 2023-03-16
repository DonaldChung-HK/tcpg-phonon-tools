
"""
TCPG phonon tools: some tools for phonon using using ASE and dftb+ and xtb
"""

from __future__ import absolute_import
from setuptools import setup, find_packages
# from os.path import abspath, dirname
def main():
    """Install the package using setuptools"""
    # project_dir = abspath(dirname(__file__))

    setup(
        name='tcpg_phonon_tools',
        version="0.0.1",
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
        install_requires=['ase', 'spglib', 'matplotlib','pandas', 'scipy','matplotlib', 'phonopy'],
        entry_points={
            'console_scripts': [
                'dftbp-socket-setup = tcpg_phonon_tools.dftbp.dftbp_socket_hsd_setup:main',
                'ase-get-primitive = tcpg_phonon_tools.ase_get_primitive:main',
                'ase-dftbp-optimise = tcpg_phonon_tools.dftbp.ase_dftbp_optimise:main',
                'abins-mantid-runner = tcpg_phonon_tools.Analysis.MantidHelper:cli_mantid_abins',
                'emd-chart = tcpg_phonon_tools.Analysis.Visualise:emd_and_chart_cli'
                ]
            },
        include_package_data=True,
        )
        

if __name__ == "__main__":
    main()