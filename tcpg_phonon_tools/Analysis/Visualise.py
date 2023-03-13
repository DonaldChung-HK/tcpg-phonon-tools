import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.stats import wasserstein_distance

import argparse

def emd_and_chart(
        calculated_path = "foo_abins.csv",
        ref_path = "foo_ref.csv", 
        cut_off = 4000.0,
        out_path = "fig.png"
    ):
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref_x = ref.iloc[:,0].to_numpy()
    ref_y = ref.iloc[:,1].to_numpy()
    auc_ref = np.trapz(ref_y, ref_x)
    ref_y = ref_y / auc_ref

    calculated = pd.read_csv(calculated_path)
    calculated = calculated[calculated.iloc[:,0] >= ref_x.min()]
    calculated = calculated[calculated.iloc[:,0] <= ref_x.max()]
    calculated_x = calculated.iloc[:,0].to_numpy()
    calculated_y = calculated.iloc[:,1].to_numpy()
    auc_calculated = np.trapz(calculated_y, calculated_x)
    calculated_y = calculated_y / auc_calculated

    f = interp1d(ref_x, ref_y, kind='cubic')

    ref_y_inter = f(calculated_x)



    emd_calc = wasserstein_distance(ref_y_inter, calculated_y)
    emd_string = f"EMD: {emd_calc}"
    print(emd_string)

    print("Drawing Fig")
    
    fig, axes = plt.subplots(edgecolor='#ffffff')
    axes.plot(calculated_x, ref_y_inter, color='#1f77b4', label='Ref')
    axes.plot(calculated_x, calculated_y,color='#ff7f0e', label='Calc')
    axes.set_title(f'Ref VS Calculated | {emd_string}')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([0, 4000])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend

    if out_path != None:
        plt.savefig(out_path)
    else:
        plt.show()

def emd_and_chart_cli():
    parser = argparse.ArgumentParser(
        description="""
        conduct 1-d interpolation of reference y using f(calculated_x) scale the arbitry y using area under curve, calculate the earth moving distance and visulalise the spectrum
        """
    )
    parser.add_argument(
        '-c',
        '--calculated_input_file', 
        type=str,
        help="Path to csv of the abins calculation this should be one with smaller but even bin Width for interpolation"
    )
    parser.add_argument(
        '-r',
        '--ref_input_file',
        type=str,
        help="reference spectrum usually from ins_db in csv"
    )
    parser.add_argument(
        '-x',
        '--cut_off',
        type=float,
        default=4000.0,
        help="upper bound of the spectrum"
    )
    parser.add_argument(
        '-o',
        '--out_path',
        type=str,
        default=None,
        help="out_put the figure to file. If None it will just show it"
    )
    args = parser.parse_args()

    emd_and_chart(
        calculated_path = args.calculated_input_file,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        out_path = args.out_path        
    )
