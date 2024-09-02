from collections import defaultdict
from pathlib import Path

import typer
from deltapd import file_md5, PyDeltaPD, PyDistMatrix, PyParams, PyLinearModelCorr, PyLinearModelError, PyLinearModelType
from deltapd.tree import Tree as DPDTree
import matplotlib.pyplot as plt
import time

import numpy as np
from ete3 import Tree, TreeStyle

app = typer.Typer()

def log(msg: str):
    print(f'{time.ctime()}: {msg}')

def read_meta(path: Path):
    print('Update the read meta to be consistent')

    out = dict()
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f.readlines():
            columns = line.strip().split('\t')
            out[columns[0]] = {k: v for k, v in zip(header, columns)}
            # sid, gid, contig_len = line.strip().split('\t')
            # out[gid] = (sid, int(contig_len))
    return out

@app.command()
def run(name: str):

    # path_ref = Path('/Users/aaron/phd/DeltaPD/examples/reference.tree')
    # path_qry = Path('/Users/aaron/phd/DeltaPD/examples/query.tree')
    # path_meta = Path('/Users/aaron/phd/DeltaPDNew/example/metadata.tsv')


    # Temp
    path_ref = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.tree')
    path_qry = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.tree')
    path_meta = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv')

    ref_dm_path = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.dm')
    qry_dm_path = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.dm')


    # ROOT_DIR = Path('/Users/aaron/phd/DeltaPDNew/data')
    # # ROOT_DIR = Path('/srv/home/uqamussi/projects/deltapd/code/data')
    #
    # path_ref = ROOT_DIR / 'bac120_ref_tree.tree'
    # path_qry = ROOT_DIR / 'bac120_ssu_tree.tree'
    # path_meta = ROOT_DIR / 'bac120_ssu_reps_metadata.tsv'

    log('Reading metadata file')
    d_meta = read_meta(path_meta)

    # Output directory
    out_dir = Path('/tmp/deltapd')
    out_dir.mkdir(parents=True, exist_ok=True)

    log('REF: Getting nodes and edges for reference')
    ref = DPDTree(path_ref)
    log('REF: Found reference nodes')
    ref_nodes = ref.get_nodes_and_edges_for_deltapd()
    log('REF: Creating reference DM')
    ref_dm = PyDistMatrix(ref_nodes[0], ref_nodes[1])
    ref_dm.to_file(str(ref_dm_path.absolute()))

    log('QRY: Getting nodes and edges for query')
    qry = DPDTree(path_qry)
    qry_nodes = qry.get_nodes_and_edges_for_deltapd()
    log('QRY: Found query nodes')
    qry_dm = PyDistMatrix(qry_nodes[0], qry_nodes[1])
    log('QRY: Creating query DM')
    qry_dm.to_file(str(qry_dm_path.absolute()))

    log('Creating the DeltaPD object')
    # path_meta = ROOT_DIR / 'bac120_ssu_reps_metadata_rs.tsv'
    dpd = PyDeltaPD(qry_dm, ref_dm, str(path_meta.absolute()), '\t')

    log('Running DeltaPD')
    params = PyParams(50, 1, PyLinearModelType.RepeatedMedian, PyLinearModelError.RMSE, PyLinearModelCorr.Pearson)

    log("getting results")
    results = dpd.run(params)

    log('creating plots')
    # load_as_ete3(path_ref, path_qry)
    create_embedded_plots(results)

    log('writing results')
    write_results(results, d_meta)


    """
    The goal for today is to get the output data from DeltaPD correctly.
    Ideally this would be in a format that is compatible with snakemake.
    
    - For each query taxon we find the 100 nearest neighbours
    - We then get the corresponding reference distances
    - We create a base model
    - We then test the model by removing each point in the dataset to calculate the relative influence
    
    I kind of feel like it should be doing the 100 nearest neighbours for the reference tree 
    as it would then sort of ensure that poor parts of the query tree (i.e. neighbours are all on long 
    branches) are not influencing the model.
    
    Output data:
    
        - 
        - X/Y points for KNN search + labels
        - Model params (grad, intercept, corr, err)
        
        - Relative influence of each point
        - Standard error of each point
        ^ Both of these contain:
             - Model params (grad, intercept, corr, err)
             - Taxon removed
             - Std error, relative influence
             
    
    """

    return

def create_embedded_plots(results):

    fig, ax = plt.subplots()

    return



def load_tree_as_ete3(path_tree):
    with path_tree.open() as f:
        content = f.read()
    t = Tree(content, format=1, quoted_node_names=True)
    return t

def load_as_ete3(path_ref_tree, path_qry_tree):
    # ref_tree = load_tree_as_ete3(path_ref_tree)
    t = load_tree_as_ete3(path_qry_tree)

    ts = TreeStyle()
    ts.show_leaf_name = True

    t.show(tree_style=ts)

    return


def write_results(results, d_meta):

    # fig, ax = plt.subplots()

    # From the results we want to get the relative influence of each point
    d_taxon_relative_inf = defaultdict(list)


    # Create plots of the points and the line of best fit
    for result in results:
        # taxon_removed = result.taxon_removed

        # qry_data = result.knn_as_linalg.qry_data
        qry_labels = result.knn_as_linalg.qry_labels
        # ref_data = result.knn_as_linalg.ref_data
        # ref_labels = result.knn_as_linalg.ref_labels
        #
        # base_corr = result.linear_model_base.eval.corr
        # base_err = result.linear_model_base.eval.error
        #
        # base_grad = result.linear_model_base.params.gradient
        # base_int = result.linear_model_base.params.intercept

        # for cur_ref, cur_qry, cur_jk in zip(ref_labels, qry_labels, result.linear_model_jk):
        #     print()

        relative_inf = result.relative_influence
        std_err = result.std_error

        for cur_qry, cur_inf in zip(qry_labels, std_err):
            d_taxon_relative_inf[cur_qry].append(cur_inf)

    d_results = dict()

    for cur_taxon, cur_vec in sorted(d_taxon_relative_inf.items()):
        cur_mean = np.mean(cur_vec)
        cur_median = np.median(cur_vec)
        cur_std = np.std(cur_vec)
        d_results[cur_taxon] = (cur_mean, cur_median, cur_std)
        # print(f'{cur_taxon} = {cur_mean:.4f}+/-{cur_std:.4f} (median={cur_median:.4f})')

    with open('/tmp/deltapd.tsv', 'w') as f:
        f.write('taxon\tmean\tmedian\tstd\tcontig_len\n')
        for cur_taxon, cur_vec in d_results.items():
            contig_len = str(d_meta[cur_taxon])
            f.write(f'{cur_taxon}\t{cur_vec[0]:.4f}\t{cur_vec[1]:.4f}\t{cur_vec[2]:.4f}\t{contig_len}\n')

    # plt.show()

    return



@app.command()
def goodbye(name: str, formal: bool = False):
    if formal:
        print(f"Goodbye Ms. {name}. Have a good day.")
    else:
        print(f"Bye {name}!")


if __name__ == '__main__':
    app()
