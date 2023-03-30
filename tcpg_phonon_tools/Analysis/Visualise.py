import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.stats import wasserstein_distance

import argparse

from pathlib import Path

def emd_and_chart(
        calculated_path = "foo_abins.csv",
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "fig.png"
    ):
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref = ref[ref[ref.columns[0]] >= min_energy]
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

    y_top = max(max(ref_y_inter), max(calculated_y)) * 1.05

    emd_calc = wasserstein_distance(ref_y_inter, calculated_y) * (cut_off - min_energy) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
    emd_string = f"EMD: {emd_calc}"
    print(emd_string)

    print("Drawing Fig")
    
    fig, axes = plt.subplots(edgecolor='#ffffff')
    axes.plot(calculated_x, ref_y_inter, color='#1f77b4', label='Ref')
    axes.plot(calculated_x, calculated_y,color='#ff7f0e', label='Calc')
    axes.set_title(f'Ref VS Calculated | {emd_string}')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([0, cut_off])
    axes.set_ylim([0, y_top])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend

    if out_path != None:
        plt.savefig(out_path)
    else:
        plt.show()

    return emd_calc

def emd_and_chart_multi(
        calculated_paths = ["foo_abins.csv"],
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "fig.png"
    ):
    emd_result = []
    fig, axes = plt.subplots(edgecolor='#ffffff')
    
    print(ref_path)
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref = ref[ref[ref.columns[0]] >= min_energy]
    ref_x = ref.iloc[:,0].to_numpy()
    ref_y = ref.iloc[:,1].to_numpy()
    auc_ref = np.trapz(ref_y, ref_x)
    ref_y = ref_y / auc_ref

    f = interp1d(ref_x, ref_y, kind='cubic')
    y_top = 0
    for calculated_path in calculated_paths:
        calculated_path = Path(calculated_path) if not isinstance(calculated_path, Path) else calculated_path
        calculated = pd.read_csv(calculated_path)
        calculated = calculated[calculated.iloc[:,0] >= ref_x.min()]
        calculated = calculated[calculated.iloc[:,0] <= ref_x.max()]
        calculated_x = calculated.iloc[:,0].to_numpy()
        calculated_y = calculated.iloc[:,1].to_numpy()
        auc_calculated = np.trapz(calculated_y, calculated_x)
        calculated_y = calculated_y / auc_calculated
        if max(calculated_y) * 1.05 > y_top:
            y_top = max(calculated_y) * 1.05

        ref_y_inter = f(calculated_x) #assuming you put stuff with the same BinWidth

        emd_calc = wasserstein_distance(ref_y_inter, calculated_y) * (cut_off - min_energy) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
        emd_string = f"{calculated_path.name} EMD: {emd_calc}"
        axes.plot(calculated_x, calculated_y, label=calculated_path.name)
        print(emd_string)
        emd_result.append([calculated_path.name, emd_calc])
    
    y_top = max(y_top, max(ref_y_inter)* 1.05)

    axes.plot(calculated_x, ref_y_inter, label='Ref')
    axes.set_title('Ref VS Calculated')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([min_energy, cut_off])
    axes.set_ylim([0, y_top])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend

    emd_result_df = pd.DataFrame(emd_result, columns=['Label', 'EMD'])

    if out_path != None:
        plt.savefig(out_path)
        emd_result_df.to_csv(f'{out_path}.csv')
    else:
        print(emd_result_df)
        plt.show()

    return emd_result_df

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
        '-m',
        '--min_energy',
        type=float,
        default=0.0,
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
        min_energy=args.minmin_energy,
        out_path = args.out_path        
    )

def emd_and_chart_multi_cli():
    parser = argparse.ArgumentParser(
        description="""
        conduct 1-d interpolation of reference y using f(calculated_x) scale the arbitry y using area under curve, calculate the earth moving distance and visulalise the spectrum
        """
    )
    parser.add_argument(
        '-c',
        '--calculated_input_files', 
        type=str,
        nargs='+',
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
        '-m',
        '--min_energy',
        type=float,
        default=0.0,
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

    emd_and_chart_multi(
        calculated_paths = args.calculated_input_files,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path        
    )
