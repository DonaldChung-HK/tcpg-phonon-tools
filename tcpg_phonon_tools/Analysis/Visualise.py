import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.metrics import mean_squared_error,r2_score
import math

from scipy.interpolate import interp1d
from scipy.stats import wasserstein_distance, pearsonr
from scipy.signal import butter, filtfilt

import argparse

from pathlib import Path

import sys

import re

# reserving '#1f77b4' for reference
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler('color', ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
# reference: https://gist.github.com/thriveth/8560036
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

def color_picker(name):
    name = str(name)
    castep_pattern = re.compile("CASTEP")
    crystal_pattern = re.compile("CRYSTAL")
    gfn1_pattern = re.compile("GFN1")
    gfn2_pattern = re.compile("GFN2")
    if castep_pattern.search(name):
        return '#ff7f00'
    elif crystal_pattern.search(name):
        return '#e41a1c'
    elif gfn1_pattern.search(name):
        return '#984ea3'
    elif gfn2_pattern.search(name):
        return '#a65628'
    else:
        return '#999999'

def label_picker(name):
    castep_pattern = re.compile("CASTEP")
    crystal_pattern = re.compile("CRYSTAL")
    gfn1_pattern = re.compile("GFN1")
    gfn2_pattern = re.compile("GFN2")
    name = str(name)
    if castep_pattern.search(name):
        return 'CASTEP'
    elif crystal_pattern.search(name):
        return 'CRYSTAL'
    elif gfn1_pattern.search(name):
        return 'GFN1'
    elif gfn2_pattern.search(name):
        return 'GFN2'
    else:
        return 'Unknown'

def emd_and_chart(
        calculated_path = "foo_abins.csv",
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "fig.png",
        do_cumulative_emd = True,
        cumulative_emd_data = None,
        do_error_range = True,
    ):
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref = ref[ref[ref.columns[0]] >= min_energy]
    ref_x = ref.iloc[:,0].to_numpy()
    ref_y = ref.iloc[:,1].to_numpy()
    try:
        ref_error = ref.iloc[:,2].to_numpy()
    except:
        print("No error data in reference")
        do_error_range = False
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

    emd_calc = wasserstein_distance(ref_y_inter, calculated_y) * (max(calculated_x) - min(calculated_x)) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
    emd_string = f"EMD: {emd_calc}"
    print(emd_string)
    cumulative_emd = []
    if do_cumulative_emd:
        print("Doing cumulative EMD:")
        for i in range(1, len(calculated_x) + 1):
            sys.stdout.write(f"{i/(len(calculated_x) + 1)*100:.2f}%")
            sys.stdout.write('\r')
            ref_y_inter_current = ref_y_inter[0:i]
            calculated_y_current = calculated_y[0:i]
            calculated_x_current = calculated_x[0:i]
            emd_current = wasserstein_distance(ref_y_inter_current, calculated_y_current) * (max(calculated_x_current) - min(calculated_x_current))
            cumulative_emd.append([calculated_x[i-1], emd_current])
        cumulative_emd_df = pd.DataFrame(cumulative_emd, columns=['x', 'CumEMD'])
    elif cumulative_emd_data != None:
        cumulative_emd_df = pd.read_csv(cumulative_emd_data)

    print("Drawing Fig")



    if cumulative_emd_data != None or do_cumulative_emd:
        fig, axes = plt.subplots(edgecolor='#ffffff')
        axes2 = axes.twinx()
        axes2.plot(cumulative_emd_df.iloc[:,0], cumulative_emd_df.iloc[:,1],color="brown", label='Cumulative_emd', linestyle='-', alpha=0.6, linewidth=0.6)
        axes2.set_ylabel('Cumulative EMD')
    else:
        fig, axes = plt.subplots(edgecolor='#ffffff')
    axes.plot(calculated_x, ref_y_inter,color='#1f77b4', label='Ref',  linestyle='-', alpha=0.7, linewidth=0.7)
    axes.plot(calculated_x, calculated_y,color='#ff7f0e', label='Calc',  linestyle='-', alpha=0.7, linewidth=0.7)
    
    if do_error_range:
        ref_error = ref_error / auc_ref
        # plotting the error range is basically a linear interpolation
        # need to thin about this
        f_error = interp1d(ref_x, ref_error, kind='cubic')
        ref_error_inter = f_error(calculated_x)
        
        axes.fill_between(calculated_x, ref_y_inter-ref_error_inter, ref_y_inter+ref_error_inter, color='lightpink',alpha=0.6)

        # axes.fill_between(ref_x, ref_y-ref_error, ref_y+ref_error, color='#1f77b4',alpha=0.4, interpolate="cubic")


    axes.set_title(f'Ref VS Calculated | {emd_string}')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([0, cut_off])
    axes.set_ylim([0, y_top])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend

    if out_path != None:
        plt.tight_layout()
        plt.savefig(out_path, dpi=400)
        if do_cumulative_emd:
            cumulative_emd_df.to_csv(f"{out_path}_cumulative_emd.csv", index=False)
            command = f"emd-chart --cumulative_emd_data {out_path}_cumulative_emd.csv -c {calculated_path} -r {ref_path} -m {min_energy} -x {cut_off}"
        elif cumulative_emd_data != None:
            command = f"emd-chart --cumulative_emd_data {cumulative_emd_data} -c {calculated_path} -r {ref_path} -m {min_energy} -x {cut_off}"
        else:
            command = f"emd-chart -c {calculated_path} -r {ref_path} -m {min_energy} -x {cut_off}"


        if do_error_range:
            command += " --do_error_range"
        
        with open(f"{out_path}_plot.sh", "w") as f:
            f.write(command)
        
            
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

        emd_calc = wasserstein_distance(ref_y_inter, calculated_y) * (max(calculated_x) - min(calculated_x)) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
        emd_string = f"{calculated_path.name} EMD: {emd_calc}"
        axes.plot(calculated_x, calculated_y, label=calculated_path.name, linestyle='-', alpha=0.6, linewidth=0.7)
        print(emd_string)
        emd_result.append([calculated_path.name, emd_calc])
    
    y_top = max(y_top, max(ref_y_inter)* 1.05)

    axes.plot(calculated_x, ref_y_inter, label='Ref', alpha=0.7, color='#1f77b4', linewidth=0.7)
    axes.set_title('Ref VS Calculated')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([min_energy, cut_off])
    axes.set_ylim([0, y_top])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend

    emd_result_df = pd.DataFrame(emd_result, columns=['Label', 'EMD'])

    if out_path != None:
        plt.tight_layout()
        plt.savefig(out_path, dpi=400)
        emd_result_df.to_csv(f'{out_path}.csv')
        with open(f"{out_path}_plot.sh", "w") as f:
            f.write(f"emd-chart-multi -c {' '.join(calculated_paths)} -r {ref_path} -m {min_energy} -x {cut_off}")
    else:
        print(emd_result_df)
        plt.show()

    return emd_result_df

def chart_multi(
        calculated_paths = ["foo_abins.csv"],
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "fig.png"
    ):
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
    zeros = [0]
    y_labels = ["Ref 0"]
    y_top = 0
    first_time = True
    for calculated_path in calculated_paths:
        calculated_path = Path(calculated_path) if not isinstance(calculated_path, Path) else calculated_path
        calculated = pd.read_csv(calculated_path)
        calculated = calculated[calculated.iloc[:,0] >= ref_x.min()]
        calculated = calculated[calculated.iloc[:,0] <= ref_x.max()]
        calculated_x = calculated.iloc[:,0].to_numpy()
        if first_time:
            ref_y_inter = f(calculated_x)#assuming you put stuff with the same BinWidth
            next_zero = 0 + max(ref_y_inter)
            first_time = False
        calculated_y = calculated.iloc[:,1].to_numpy() 
        auc_calculated = np.trapz(calculated_y, calculated_x)
        calculated_y = calculated_y / auc_calculated
        calculated_y += next_zero
        if max(calculated_y) * 1.05 > y_top:
            y_top = max(calculated_y) * 1.05
        axes.plot(calculated_x, calculated_y, label=label_picker(calculated_path.parent.name), linestyle='-', alpha=0.6, linewidth=0.7, color=color_picker(calculated_path.parent.name))
        zeros.append(next_zero)
        y_labels.append(f"{label_picker(calculated_path.parent.name)} 0")
        next_zero = max(calculated_y)
    

    axes.plot(calculated_x, ref_y_inter, label='Ref', alpha=0.7, color='#1f77b4', linewidth=0.7)
    axes.set_title('Ref VS Calculated')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units')
    axes.set_yticks(zeros)
    axes.set_yticklabels(y_labels)
    axes.set_xlim([min_energy, cut_off])
    #legend = axes.legend(fontsize=8.0).set_draggable(True).legend
    axes.legend('', frameon=False)

    if out_path != None:
        plt.tight_layout()
        plt.savefig(out_path, dpi=400)
        with open(f"{out_path}_plot.sh", "w") as f:
            f.write(f"emd-chart-multi -c {' '.join(calculated_paths)} -r {ref_path} -m {min_energy} -x {cut_off}")
    else:
        plt.show()

def chart_multi_filter(
        calculated_paths = ["foo_abins.csv"],
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "fig.png",
        high_pass_cutoff = 0.001
    ):
    fig, axes = plt.subplots(edgecolor='#ffffff')
    
    print(ref_path)
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref = ref[ref[ref.columns[0]] >= min_energy]
    ref_x = ref.iloc[:,0].to_numpy()
    ref_y = ref.iloc[:,1].to_numpy()

    f = interp1d(ref_x, ref_y, kind='cubic')
    zeros = [0]
    y_labels = ["Ref 0"]
    y_top = 0
    first_time = True

    for calculated_path in calculated_paths:
        calculated_path = Path(calculated_path) if not isinstance(calculated_path, Path) else calculated_path
        calculated = pd.read_csv(calculated_path)
        calculated = calculated[calculated.iloc[:,0] >= ref_x.min()]
        calculated = calculated[calculated.iloc[:,0] <= ref_x.max()]
        calculated_x = calculated.iloc[:,0].to_numpy()
        b, a = butter(2, high_pass_cutoff, btype='high')
        if first_time:
            ref_y_inter = f(calculated_x)#assuming you put stuff with the same BinWidth
            high_passed_ref = filtfilt(b, a, ref_y_inter)
            auc_ref = np.trapz(np.abs(high_passed_ref), calculated_x)
            high_passed_ref = high_passed_ref / auc_ref
            next_zero = 0 + max(high_passed_ref)
            first_time = False
        calculated_y = calculated.iloc[:,1].to_numpy() 
        high_passed_calc = filtfilt(b, a, calculated_y)
        auc_calculated = np.trapz(np.abs(high_passed_calc), calculated_x)
        high_passed_calc = high_passed_calc / auc_calculated 
        next_zero += abs(min(min(high_passed_calc),0))
        original_max = max(high_passed_calc)  
         
        high_passed_calc += next_zero
        zeros.append(next_zero)
        axes.plot(calculated_x, high_passed_calc, label=label_picker(calculated_path.parent.name), linestyle='-', alpha=0.6, linewidth=0.7, color=color_picker(calculated_path.parent.name))
        y_labels.append(f"{label_picker(calculated_path.parent.name)} 0")
        next_zero += original_max 
    

    axes.plot(calculated_x, high_passed_ref, label='Ref', alpha=0.7, color='#1f77b4', linewidth=0.7)
    #axes.set_title('Ref VS Calculated')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units')
    axes.set_yticks(zeros)
    axes.set_yticklabels(y_labels)
    axes.set_xlim([min_energy, cut_off])
    #legend = axes.legend(fontsize=8.0).set_draggable(True).legend
    axes.legend('', frameon=False)

    if out_path != None:
        plt.tight_layout()
        plt.savefig(out_path, dpi=400)
        with open(f"{out_path}_plot.sh", "w") as f:
            f.write(f"emd-chart-multi -c {' '.join(calculated_paths)} -r {ref_path} -m {min_energy} -x {cut_off}")
    else:
        plt.show()
def chart_multi_filter_cli():
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
    parser.add_argument(
        '--high_pass_cutoff',
        type=float,
        default=0.001,
        help="cutoff for high pass filter"
    )    
    args = parser.parse_args()

    chart_multi_filter(
        calculated_paths = args.calculated_input_files,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,
        high_pass_cutoff=args.high_pass_cutoff,     
    )


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
    parser.add_argument(
        '--do_cumulative_emd',
        action="store_true",
        help="calculate and show cumulative emd in seondary y axis"
    )
    parser.add_argument(
        '--do_error_range',
        action="store_true",
        help="calculate and show error range"
    )
    parser.add_argument(
        '--cumulative_emd_data',
        type=str,
        default=None,
        help="use calculated cumulative_emd data"
    )
    args = parser.parse_args()

    emd_and_chart(
        calculated_path = args.calculated_input_file,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,
        do_cumulative_emd=args.do_cumulative_emd,
        cumulative_emd_data=args.cumulative_emd_data,
        do_error_range =args.do_error_range,      
    )

def chart_multi_cli():
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

    chart_multi(
        calculated_paths = args.calculated_input_files,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,     
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

    chart_multi(
        calculated_paths = args.calculated_input_files,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,     
    )


def quant_metrics(
        calculated_path = "foo_abins.csv",
        ref_path = "foo_ref.csv", 
        min_energy = 0.0, 
        cut_off = 4000.0,
        out_path = "metrics.csv",
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

    #y_top = max(max(ref_y_inter), max(calculated_y)) * 1.05

    emd_calc = wasserstein_distance(ref_y_inter, calculated_y) * (max(calculated_x) - min(calculated_x)) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
    emd_string = f"EMD: {emd_calc}"
    print(emd_string)
    mse = mean_squared_error(ref_y_inter, calculated_y) * (max(calculated_x) - min(calculated_x)) # normalising the weight to the range for easier number
    rmse = math.sqrt(mse)
    print(f"Root Mean Square Error:{rmse}")
    pcc = pearsonr(ref_y_inter, calculated_y)
    pcc_score = pcc[0]
    print(f"Pearson correlation coefficient: {pcc_score}")
    r2 = r2_score(ref_y_inter, calculated_y)
    print(f"r2 score: {r2}")

    result_dict = {
        "EMD":emd_calc,
        "RMSE":rmse,
        "PCC":pcc_score,
        "R2":r2,
    }
    result = pd.DataFrame.from_records(list(result_dict.items()),columns=["Metrics", "Quantity"])
    if out_path != None:
        result.to_csv(out_path, index=False)
    
    return result_dict

def metrics_cli():
    parser = argparse.ArgumentParser(
        description="""
        conduct 1-d interpolation of reference y using f(calculated_x) scale the arbitry y using area under curve, calculate the quantitative metrics
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

    quant_metrics(
        calculated_path = args.calculated_input_file,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,
    )



def emd_and_chart_filter(
    calculated_path = "foo_abins.csv",
    ref_path = "foo_ref.csv", 
    min_energy = 0.0, 
    cut_off = 4000.0,
    out_path = "fig.png",
    high_pass_cutoff = 0.001,
    supress_output = False
    ):
    ref = pd.read_csv(ref_path)
    ref = ref[ref[ref.columns[0]] <= cut_off]
    ref = ref[ref[ref.columns[0]] >= min_energy]
    ref_x = ref.iloc[:,0].to_numpy()
    ref_y = ref.iloc[:,1].to_numpy()

    calculated = pd.read_csv(calculated_path)
    calculated = calculated[calculated.iloc[:,0] >= ref_x.min()]
    calculated = calculated[calculated.iloc[:,0] <= ref_x.max()]
    calculated_x = calculated.iloc[:,0].to_numpy()
    calculated_y = calculated.iloc[:,1].to_numpy()

    f = interp1d(ref_x, ref_y, kind='cubic')

    ref_y_inter = f(calculated_x)

   
    b, a = butter(2, high_pass_cutoff, btype='high')
    high_passed_ref = filtfilt(b, a, ref_y_inter)
    high_passed_calc = filtfilt(b, a, calculated_y)
    auc_calculated = np.trapz(np.abs(high_passed_calc), calculated_x)
    high_passed_calc = high_passed_calc / auc_calculated
    aus_ref = np.trapz(np.abs(high_passed_ref), calculated_x)
    high_passed_ref = high_passed_ref / aus_ref

    y_top = max(max(high_passed_ref), max(high_passed_calc)) * 1.05
    y_min = min(min(high_passed_ref), min(high_passed_calc)) * 0.95
    emd_calc = wasserstein_distance(high_passed_ref, high_passed_calc) * (max(calculated_x) - min(calculated_x)) #the default weight is normalised to 1 so multiply it by the range will have a more readable number
    emd_string = f"EMD: {emd_calc}"
    print(emd_string)
    mse = mean_squared_error(high_passed_ref, high_passed_calc) * (max(calculated_x) - min(calculated_x)) # normalising the weight to the range for easier number
    rmse = math.sqrt(mse)
    print(f"Root Mean Square Error:{rmse}")
    pcc = pearsonr(high_passed_ref, high_passed_calc)
    pcc_score = pcc[0]
    print(f"Pearson correlation coefficient: {pcc_score}")
    r2 = r2_score(high_passed_ref, high_passed_calc)
    print(f"r2 score: {r2}")

    print("Drawing Fig")

    fig, axes = plt.subplots(edgecolor='#ffffff')
    axes.plot(calculated_x, high_passed_ref,color='#1f77b4', label='Ref',  linestyle='-', alpha=0.7, linewidth=0.7)
    axes.plot(calculated_x, high_passed_calc,color='#ff7f0e', label='Calc',  linestyle='-', alpha=0.7, linewidth=0.7)

    axes.set_title(f'Ref VS Calculated | {emd_string}')
    axes.set_xlabel('Energy transfer ($cm^{-1}$)')
    axes.set_ylabel('S / Arbitrary Units ($cm^{-1}$)$^{-1}$')
    axes.set_xlim([0, cut_off])
    axes.set_ylim([y_min, y_top])
    legend = axes.legend(fontsize=8.0).set_draggable(True).legend
    if not supress_output:
        if out_path != None:
            plt.tight_layout()
            plt.savefig(out_path, dpi=400)
            command = f"emd-chart-filter -c {calculated_path} -r {ref_path} -m {min_energy} -x {cut_off}"
            
            with open(f"{out_path}_plot.sh", "w") as f:
                f.write(command)
        else:
            plt.show()
    else:
        plt.cla()
    return (high_pass_cutoff, emd_calc, rmse, pcc_score, r2)

def emd_and_chart_filter_cli():
    parser = argparse.ArgumentParser(
        description="""
        conduct 1-d interpolation of reference y using f(calculated_x), perform a high pass filter using cut off scale the arbitry y using area under curve, calculate the earth moving distance and visulalise the spectrum
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

    parser.add_argument(
        '--high_pass_cutoff',
        type=float,
        default=0.0,
        help="cutoff for high pass filter"
    )    

    args = parser.parse_args()

    emd_and_chart_filter(
        calculated_path = args.calculated_input_file,
        ref_path = args.ref_input_file, 
        cut_off = args.cut_off,
        min_energy=args.min_energy,
        out_path = args.out_path,
        high_pass_cutoff = args.high_pass_cutoff, 
    )