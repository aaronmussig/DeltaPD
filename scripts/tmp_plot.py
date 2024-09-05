from collections import defaultdict
from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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
            taxon, r_inf, std_err = line.strip().split('\t')
            r_inf = float(r_inf)
            std_err = float(std_err)
            out[taxon] = (r_inf, std_err)
    return out

def aggregate_inf_files(d_inf):
    out = dict()
    for taxon_id, d_replicates in d_inf.items():
        out[taxon_id] = dict()
        for replicate, path in d_replicates.items():
            cur_data = read_inf_file(path)
            out[taxon_id][replicate] = cur_data
    return out

def plot_data(xy_files, inf_files):

    # The first axis will contain the X/Y points seried by the replicate
    # The second axis will contain a heatmap of the influence values
    # 3rd will be heatmap of std error

    REPLICATE_LIMIT = None

    for taxon, d_xy_replicate in xy_files.items():

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 10))

        # Collect the XY data for each replicate and plot it as a different series
        text_values_to_add = list()
        for replicate, d_xy in sorted(d_xy_replicate.items()):
            x = [v[0] for v in d_xy.values()]
            y = [v[1] for v in d_xy.values()]


            for (qi, qj), (xi, xj) in d_xy.items():
                if qi == taxon or qj == taxon:
                    text_values_to_add.append((xi, xj, f'{qi}\n{qj}'))
            ax1.scatter(x, y, label=f'Replicate {replicate}', alpha=0.2, s=1)

            if REPLICATE_LIMIT and replicate > REPLICATE_LIMIT:
                break

        # Add the text values wherever the taxon of interest appears
        for xi, xj, text in text_values_to_add:
            ax1.text(xi, xj, '.', size=4)

        d_inf_replicate = inf_files[taxon]
        inf_values_lst = list()
        std_values_lst = list()
        for replicate, d_inf in sorted(d_inf_replicate.items()):
            cur_inf, cur_std = d_inf[taxon]
            inf_values_lst.append(cur_inf)
            std_values_lst.append(cur_std)

            if REPLICATE_LIMIT and replicate > REPLICATE_LIMIT:
                break


        ax2.hist(inf_values_lst)
        ax3.hist(std_values_lst)

        print(f'Mean std: {np.mean(std_values_lst)}')
        print(f'Median stf: {np.median(std_values_lst)}')

        ax2.set_title('Relative Influence')
        ax3.set_title('Standard Error')

        # Add a title above the two subplots
        fig.suptitle(f'{taxon}')
        plt.show()
        # plt.savefig(f'/tmp/dpd_figs/{taxon}.png')
        plt.close()
        print('Saved', taxon)


    return


def main():

    print('Locating files')
    d_xy, d_inf = find_files()

    print('Getting XY files')
    xy_files = aggregate_xy_files(d_xy)

    print('Getting INF files')
    inf_files = aggregate_inf_files(d_inf)

    print('Plotting data')
    plot_data(xy_files, inf_files)

    return


if __name__ == '__main__':
    main()