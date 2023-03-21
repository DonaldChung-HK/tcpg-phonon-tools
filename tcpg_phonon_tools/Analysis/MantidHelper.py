from mantid.simpleapi import Abins
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mantid.api import AnalysisDataService as ADS
from pathlib import Path
import argparse

def run_abins(
    label = "FOO",
    input_file = "phonopy.yaml",
    qpt_cutoff = 15.0,
    AbInitioProgram='FORCECONSTANTS',
    SumContributions=True,
    ScaleByCrossSection='Total',
    QuantumOrderEventsNumber = '2',
    Autoconvolution = True,
    Setting = 'All detectors (TOSCA)',
    BinWidthInWavenumber = 1,
    retrieve_items = ["total"],
    out_path = Path("."),
    **kwarg
    ):
    if not isinstance(out_path, Path):
        out_path = Path(out_path)
    
    import abins.parameters

    abins.parameters.sampling['force_constants']['qpt_cutoff'] = qpt_cutoff
    
    Abins(
        VibrationalOrPhononFile=input_file, 
        AbInitioProgram=AbInitioProgram, 
        OutputWorkspace=label, 
        SumContributions=SumContributions,
        ScaleByCrossSection=ScaleByCrossSection, 
        QuantumOrderEventsNumber=QuantumOrderEventsNumber, 
        Autoconvolution=Autoconvolution, 
        Setting=Setting, 
        BinWidthInWavenumber=BinWidthInWavenumber,
        **kwarg
    )

    for item in retrieve_items:
        retrieve_workspace = f'{label}_{item}'
        current = ADS.retrieve(retrieve_workspace)
        abin_y = current.extractY()[0]
        abin_x = current.extractX()[0]
        abin_x = (abin_x[1:] + abin_x[:-1]) / 2 
        abins_data = pd.DataFrame(np.column_stack((abin_x, abin_y)), columns=['E cm-1', 'count'])
        out_file_name = f"{retrieve_workspace}_abins.csv"
        abins_data.to_csv(out_path / out_file_name, index=False)


def cli_mantid_abins():
    parser = argparse.ArgumentParser(
        description="This script will run abins through python api and save the result to a csv file")
    parser.add_argument(
        '-f',
        '--input_file',
        type=str,
        default='phonopy.yaml',
        help="Input cell file")
    parser.add_argument(
        '-l',
        '--label',
        type=str,
        default='FOO',
        help="seed name of for file name")
    parser.add_argument(
        '-p',
        '--AbInitioProgram',
        type=str,
        default='FORCECONSTANTS',
        help="DFT program which was used for a phonon calculation. Allowed values: [FORCECONSTANT, CASTEP, CRYSTAL]")
    parser.add_argument(
        '-c',
        '--NotSumContributions',
        action="store_false",
        help="Not Sum the partial dynamical structure factors into a single workspace.")
    parser.add_argument(
        '-r',
        '--ScaleByCrossSection',
        type=str,
        default='Total',
        help="Scale the partial dynamical structure factors by the scattering cross section. Allowed values: [‘Total’, ‘Incoherent’, ‘Coherent’]")
    parser.add_argument(
        '-n',
        '--QuantumOrderEventsNumber',
        type=str,
        default="2",
        help="Number of quantum order effects included in the calculation (1 -> FUNDAMENTALS, 2-> first overtone + FUNDAMENTALS + 2nd order combinations, 3-> FUNDAMENTALS + first overtone + second overtone + 2nd order combinations + 3rd order combinations etc...). Allowed values: [‘1’, ‘2’, ‘3’, ‘4’]")
    parser.add_argument(
        '-a',
        '--NotAutoconvolution',
        action="store_false",
        help="Not Autoconvolution in Abins")
    parser.add_argument(
        '-s',
        '--Setting',
        type=str,
        default='All detectors (TOSCA)',
        help="Setting variable for Abins e.g 'All detectors (TOSCA)'")
    parser.add_argument(
        '-w',
        '--BinWidthInWavenumber',
        type=float,
        default=0.2,
        help="Bin width number you mush patch the Abins to remove resutiction for values to be 1.0 - 10.0")
    parser.add_argument(
        '-i',
        '--retrieve_items',
        type=str,
        nargs='+',
        default=["Total"],
        help="the spectra name you want to retrieve like 'total', 'H1' etc.")
    
    parser.add_argument(
        '-o',
        '--out_path',
        type=str,
        default='.',
        help="path to output the files")
    parser.add_argument(
        '-q',
        '--qpt_cutoff',
        type=float,
        default=15.0,
        help="qpoint cutoff")

    
    args = parser.parse_args()
    
    run_abins(
        label = args.label,
        input_file = args.input_file,
        AbInitioProgram=args.AbInitioProgram,
        SumContributions=args.NotSumContributions,
        ScaleByCrossSection=args.ScaleByCrossSection,
        QuantumOrderEventsNumber = str(args.QuantumOrderEventsNumber),
        Autoconvolution = args.NotAutoconvolution,
        Setting = args.Setting,
        BinWidthInWavenumber = args.BinWidthInWavenumber,
        retrieve_items = args.retrieve_items,
        out_path = args.out_path,
    )