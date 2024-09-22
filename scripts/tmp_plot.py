from collections import defaultdict
from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats

import matplotlib.pyplot as plt
from scipy.stats import siegelslopes

DIR = Path('/tmp/dpd')


RE_XY = re.compile(r'deltapd_xy_(.+)_(\d+)')
RE_INF = re.compile(r'deltapd_relative_influence_(.+)_(\d+)')

def find_files():
    out_xy = dict()
    out_inf = dict()
    for file in DIR.glob('*.tsv'):
        if file.stem.startswith('deltapd_xy'):
            hits = RE_XY.match(file.stem)
            taxon_id, replicate = hits.groups()
            replicate = int(replicate)
            if taxon_id not in out_xy:
                out_xy[taxon_id] = dict()
            out_xy[taxon_id][replicate] = file
        elif file.stem.startswith('deltapd_relative_influence'):
            hits = RE_INF.match(file.stem)
            taxon_id, replicate = hits.groups()
            replicate = int(replicate)
            if taxon_id not in out_inf:
                out_inf[taxon_id] = dict()
            out_inf[taxon_id][replicate] = file
    return out_xy, out_inf

def read_xy_file(path: Path):
    out = dict()
    with path.open() as f:
        f.readline()
        for line in f.readlines():
            qry_taxon_i, qry_taxon_j, ref_taxon_i, ref_taxon_j, qry_dist, ref_dist = line.strip().split('\t')
            qry_dist, ref_dist = float(qry_dist), float(ref_dist)
            out[(qry_taxon_i, qry_taxon_j)] = (qry_dist, ref_dist)
    return out

def aggregate_xy_files(d_xy):
    out = dict()
    for taxon_id, d_replicates in d_xy.items():
        out[taxon_id] = dict()
        for replicate, path in d_replicates.items():
            cur_data = read_xy_file(path)
            out[taxon_id][replicate] = cur_data
    return out

def read_inf_file(path: Path):
    out = dict()
    with path.open() as f:
        f.readline()
        for line in f.readlines():
            taxon, r_inf, std_err, grad, intercept = line.strip().split('\t')
            r_inf = float(r_inf)
            std_err = float(std_err)
            grad = float(grad)
            intercept = float(intercept)
            out[taxon] = (r_inf, std_err, grad, intercept)
    return out

def aggregate_inf_files(d_inf):
    out = dict()
    for taxon_id, d_replicates in d_inf.items():
        out[taxon_id] = dict()
        for replicate, path in d_replicates.items():
            cur_data = read_inf_file(path)
            out[taxon_id][replicate] = cur_data
    return out

def plot_data_2(xy_files, inf_files):

    # Iterate over the X/Y points collected for each taxon
    for taxon, d_xy_replicate in xy_files.items():

        # Iterate over each replicate
        for rep_id, d_xy in d_xy_replicate.items():

            # Collect all X/Y points
            all_x_values = list()
            all_y_values = list()
            all_label_values = list()
            for (qi, qj), (xi, xj) in d_xy.items():
                all_x_values.append(xj)
                all_y_values.append(xi)
                all_label_values.append(f'{qi}\n{qj}')

            d_inf_values = inf_files[taxon][rep_id]
            n_plots = min(len(d_inf_values.keys()), 25)
            n_plots = math.ceil(math.sqrt(n_plots))
            fig, axes = plt.subplots(ncols=n_plots, nrows=n_plots, figsize=(15, 15), sharex=True, sharey=True)
            axes = axes.flatten()

            inf_value_subset = sorted(d_inf_values.items(), key=lambda x: -x[1][1])[:len(axes)]

            for (unq_taxon, (cur_inf, cur_std_err, cur_grad, cur_intercept)), cur_ax in zip(inf_value_subset, axes):

                cur_x_values = list()
                cur_y_values = list()
                cur_label_values = list()

                for (qi, qj), (xi, xj) in d_xy.items():
                    if qi == unq_taxon or qj == unq_taxon:
                        cur_x_values.append(xj)
                        cur_y_values.append(xi)
                        cur_label_values.append(f'{qi}\n{qj}')

                cur_ax.scatter(all_x_values, all_y_values, alpha=1, s=3)
                cur_ax.scatter(cur_x_values, cur_y_values, label=f'{unq_taxon} {cur_std_err:.4f}', alpha=1, s=3)

                # Also plot the line of best fit
                x = np.linspace(0, max(cur_x_values), 100)
                y = cur_grad * x + cur_intercept
                cur_ax.plot(x, y, label=f'y={cur_grad:.4f}x+{cur_intercept:.6f}', alpha=0.2)

                # Put a title above each plot
                cur_ax.set_title(f'{unq_taxon} ({cur_std_err:.2%})')


                print()

            # ax.legend()
            plt.show()
            print()


        print()


    return


def plot_data(xy_files, inf_files):

    # The first axis will contain the X/Y points seried by the replicate
    # The second axis will contain a heatmap of the influence values
    # 3rd will be heatmap of std error

    REPLICATE_LIMIT = None

    for taxon, d_xy_replicate in xy_files.items():

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 10))

        # Collect the XY data for each replicate and plot it as a different series
        text_values_to_add = list()
        max_x = 0
        max_y = 0
        for replicate, d_xy in sorted(d_xy_replicate.items()):
            x = [v[1] for v in d_xy.values()]
            y = [v[0] for v in d_xy.values()]
            max_x = max(max_x, max(x))
            max_y = max(max_y, max(y))

            # Re-compute the repeated median for these values
            slope, intercept = siegelslopes(y, x)
            print(taxon, replicate, slope, intercept)

            for (qi, qj), (xi, xj) in d_xy.items():
                if qi == taxon or qj == taxon:
                    text_values_to_add.append((xi, xj, f'{qi}\n{qj}'))
            ax1.scatter(x, y, label=f'Replicate {replicate}', alpha=0.2, s=1)

            if REPLICATE_LIMIT and replicate > REPLICATE_LIMIT:
                break

        # Add the text values wherever the taxon of interest appears
        for xi, xj, text in text_values_to_add:
            ax1.text(xj, xi, 'x', size=5)

        d_inf_replicate = inf_files[taxon]
        inf_values_lst = list()
        std_values_lst = list()

        for replicate, d_inf in sorted(d_inf_replicate.items()):
            cur_inf, cur_std, cur_grad, cur_intercept = d_inf[taxon]
            inf_values_lst.append(cur_inf)
            std_values_lst.append(cur_std)

            # Plot the line of best fit
            x = np.linspace(0, max_x, 100)
            y = cur_grad * x + cur_intercept
            ax1.plot(x, y, label=f'r{replicate} {cur_grad:.4f}+{cur_intercept:.6f}', alpha=0.2)
            # ax1.legend()

            if REPLICATE_LIMIT and replicate > REPLICATE_LIMIT:
                break

        ax1.set_ylim(0, max_y)


        ax2.hist(inf_values_lst)
        ax3.hist(std_values_lst)

        print(f'Mean std: {np.mean(std_values_lst)}')
        print(f'Median stf: {np.median(std_values_lst)}')
        print(f'Total std: {np.sum(std_values_lst)}')

        ax2.set_title('Relative Influence')
        ax3.set_title('Standard Error')

        # Add a title above the two subplots
        fig.suptitle(f'{taxon}')
        plt.show()
        # plt.savefig(f'/tmp/dpd_figs/{taxon}.png')
        plt.close()
        print('Saved', taxon)


    return


def calculate_method_using_python(xy_files):

    out = dict()

    for taxon, d_xy_replicate in xy_files.items():
        out[taxon] = dict()
        for rep_id, d_xy in d_xy_replicate.items():
            out[taxon][rep_id] = dict()
            unq_taxa = sorted({xi for xi, xj in d_xy.keys()}.union({xj for xi, xj in d_xy.keys()}))

            # Create the jackknife sets
            d_taxon_to_xy = dict()
            for unq_taxon in unq_taxa:

                # Get all data points without this taxon in it
                cur_x_values, cur_y_values = list(), list()
                for (qi, qj), (xi, xj) in d_xy.items():
                    if qi == unq_taxon or qj == unq_taxon:
                        continue
                    cur_x_values.append(xj)
                    cur_y_values.append(xi)
                d_taxon_to_xy[unq_taxon] = (cur_x_values, cur_y_values)

            # Compute the model parameters using the repeated median method
            d_taxon_to_params = dict()
            for unq_taxon, (x_values, y_values) in d_taxon_to_xy.items():
                slope, intercept = siegelslopes(y_values, x_values)
                y_hat = slope * np.array(x_values) + intercept
                rmse = calc_rmse(y_values, y_hat)
                r2 = calc_r2(y_values, y_hat)
                d_taxon_to_params[unq_taxon] = (slope, intercept, rmse, r2)

            # Compute the influence values
            error_values = [d_taxon_to_params[x][2] for x in unq_taxa]
            influence_vals = calc_relative_jackknife_influence(error_values)
            std_error = calc_prop_std_error_of_jackknife(influence_vals)
            std_error_total = sum(std_error)

            for i, unq_taxon in enumerate(unq_taxa):
                slope, intercept, rmse, r2 = d_taxon_to_params[unq_taxon]
                out[taxon][rep_id][unq_taxon] = (influence_vals[i], std_error[i], slope, intercept)
    return out

def calc_rmse(y, y_i):
    return float(math.sqrt(np.mean((y - y_i) ** 2)))

def calc_r2(y, y_hat):
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1 - (ss_res / ss_tot))

def calc_relative_jackknife_influence(x):
    x_bar = np.mean(x)
    n = len(x)
    sum_squares = sum([((n - 1) * (x_bar - x_i)) ** 2 for x_i in x])
    denominator = math.sqrt(sum_squares / (n - 1))
    return [float((n - 1) * (x_bar - x_i) / denominator) for x_i in x]

def calc_prop_std_error_of_jackknife(x):
    n_inv = 1.0 / (len(x) - 1.0)
    return [float(x_i ** 2 * n_inv) for x_i in x]


def run_tests():
    law_jackknifed = [0.89294714566676303,0.76370684340942141,0.75499836753560889,0.77609676665137606,0.73131966574578589,0.77996866647699103,0.78453597565726862,0.73616182518142603,0.75173907951007512,0.7761230992718563,0.81810070810147284,0.7857184432079275,0.74035089460061443,0.76704134075534491,0.77987252287711384]
    law_jackknifed_influence = [-2.97, 0.31, 0.53, 0.0, 1.13, -0.1, -0.22, 1.01, 0.61, -0.01, -1.07, -0.25, 0.9, 0.22, -0.1]
    law_jackknifed_test = calc_relative_jackknife_influence(law_jackknifed)
    law_jackknifed_test_rounded = [round(x, 2) for x in law_jackknifed_test]

    if law_jackknifed_test_rounded != law_jackknifed_influence:
        print('Jackknife test failed')

    prop_std = calc_prop_std_error_of_jackknife(law_jackknifed_influence)
    prop_std_rounded = [round(x, 2) for x in prop_std]

    prop_std_total = round(sum(prop_std_rounded), 2)

    if prop_std_rounded[0] != 0.63:
        print('Prop std failed')

    if prop_std_total != 1.0:
        print(f'Prop std total failed: {prop_std_total} != 1.0 ')
    return




def main():

    run_tests()



    print('Locating files')
    d_xy, d_inf = find_files()

    print('Getting XY files')
    xy_files = aggregate_xy_files(d_xy)

    # inf_files_py = calculate_method_using_python(xy_files)

    print('Getting INF files')
    inf_files = aggregate_inf_files(d_inf)

    print('Plotting data')
    # plot_data_2(xy_files, inf_files)
    plot_data(xy_files, inf_files)

    return


if __name__ == '__main__':
    main()