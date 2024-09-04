import time
from collections import defaultdict
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import typer
from ete3 import Tree, faces, TreeStyle

from deltapd import PyDeltaPD, PyDistMatrix, PyParams, PyLinearModelCorr, PyLinearModelError, \
    PyLinearModelType
from deltapd.tree import Tree as DPDTree

app = typer.Typer()

GIDS_TO_PLOT = {
    "GB_GCA_016219485.1",
    "RS_GCF_901686465.1"
}


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


def read_taxonomy(path: Path):
    out = dict()
    with path.open() as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax.split(';')
    return out


def eg_ar53_reps():
    path_ref = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.tree')
    path_qry = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.tree')
    path_meta = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv')
    path_tax = Path('/Users/aaron/phd/DeltaPDNew/example/ar53_taxonomy.tsv')
    return path_ref, path_qry, path_meta, path_tax


def eg_ar53_all():
    path_ref = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53.tree')
    path_qry = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200.tree')
    path_meta = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200_metadata.tsv')
    path_tax = Path('/Users/aaron/phd/DeltaPDNew/example/ar53_taxonomy.tsv')
    return path_ref, path_qry, path_meta, path_tax


def eg_basic():
    path_ref = Path('/Users/aaron/phd/DeltaPD/examples/reference.tree')
    path_qry = Path('/Users/aaron/phd/DeltaPD/examples/query.tree')
    path_meta = Path('/Users/aaron/phd/DeltaPD/examples/metadata.tsv')
    path_tax = None
    return path_ref, path_qry, path_meta, path_tax


@app.command()
def run(name: str):
    # Load the example data

    # path_ref, path_qry, path_meta, path_tax = eg_basic()
    # path_ref, path_qry, path_meta, path_tax = eg_ar53_reps()
    path_ref, path_qry, path_meta, path_tax = eg_ar53_all()

    create_dm = True

    log('Reading metadata file')
    d_meta = read_meta(path_meta)

    # Output directory
    out_dir = Path('/tmp/deltapd')
    out_dir.mkdir(parents=True, exist_ok=True)

    log('REF: Getting nodes and edges for reference')
    ref = DPDTree(path_ref)
    log('REF: Found reference nodes')
    ref_nodes = ref.get_nodes_and_edges_for_deltapd(save=create_dm)
    log('REF: Creating reference DM')
    ref_dm = PyDistMatrix(ref_nodes[0], ref_nodes[1])

    log('QRY: Getting nodes and edges for query')
    qry = DPDTree(path_qry)
    qry_nodes = qry.get_nodes_and_edges_for_deltapd(save=create_dm)
    log('QRY: Found query nodes')
    qry_dm = PyDistMatrix(qry_nodes[0], qry_nodes[1])

    log('Creating the DeltaPD object')
    # path_meta = ROOT_DIR / 'bac120_ssu_reps_metadata_rs.tsv'
    dpd = PyDeltaPD(qry_dm, ref_dm, str(path_meta.absolute()), '\t')

    """
    Instead of taking the KNN, take the x closest points until the total branch distance is > some value
    maybe 10% of the total sum of branch lengths?
    """

    log('Running DeltaPD')  # there's an error where 0 can be all of the kNN distances # sample more??
    params = PyParams(100, 10, PyLinearModelType.RepeatedMedian, PyLinearModelError.RMSE, PyLinearModelCorr.Pearson)

    log("getting results")
    results = dpd.run(params)

    log('writing results')
    agg_results, agg_xy = write_results(results, d_meta)
    create_agg_plots(agg_results, agg_xy)

    log('creating plots')
    # load_as_ete3(qry, results, agg_results, d_taxonomy)
    create_embedded_plots(results)

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
    
    Not quite, this doesn't work for the case of multiple genes in the query tree.
    
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
             
    There is an idea of neighbourhood stability that I would like to explore. I.e. for a knn set
    what is the overall outlier / relative proportion.
    
    """

    return


def create_agg_plots(agg_results, agg_xy):
    """
    These plots are good as they show the points that only contain those for the query taxon.
    It does fall apart a bit when you look at ones that are on long branches, as they are very rarely ever
    the nearest neighbour.

    Poor linearity for these models does seem to indicate contamination.
    """

    for qry_taxon, qry_inf in sorted(agg_results.items()):
        #
        # if qry_taxon not in {'GB_GCA_030751995.1', 'GB_GCA_003164295.1', 'GB_GCA_021852135.1', 'GB_GCA_014874415.1', 'GB_GCA_014729945.1', 'GB_GCA_018659375.1', 'GB_GCA_014728835.1', 'GB_GCA_026014615.1', 'GB_GCA_002508315.1', 'GB_GCA_016839385.1', 'GB_GCA_003096255.1', 'GB_GCA_018815265.1', 'GB_GCA_018609935.1', 'GB_GCA_900318035.1', 'GB_GCA_027016915.1', 'GB_GCA_018304285.1', 'GB_GCA_016202745.1', 'GB_GCA_023254445.1', 'RS_GCF_023256325.1'}:
        #     continue

        if qry_taxon not in GIDS_TO_PLOT:
            continue

        cur_xy = agg_xy[qry_taxon]
        cur_x = [x[0] for x in cur_xy]
        cur_y = [x[1] for x in cur_xy]

        cmap = cm.get_cmap('viridis')
        norm = plt.Normalize(vmin=0, vmax=1)
        colors = cmap(norm(qry_inf))

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
        ax1.scatter(cur_x, cur_y, c=colors)
        ax1.set_title(f'KNN for: {qry_taxon}')

        ax2.hist(qry_inf)

        plt.show()
        # plt.savefig('/tmp/am/embedded.png')
        # plt.close()
        print()

    return


def create_embedded_plots(results):
    for result in sorted(results, key=lambda x: x.query_taxon):

        if result.query_taxon not in GIDS_TO_PLOT:
            continue

        fig, ax = plt.subplots()

        x_data = result.knn_as_linalg.ref_data
        y_data = result.knn_as_linalg.qry_data

        # Colour each point by the proportion of standard error
        std_error = result.std_error
        cmap = cm.get_cmap('viridis')
        norm = plt.Normalize(vmin=0, vmax=1)
        colors = cmap(norm(std_error))

        ax.scatter(x_data, y_data, c=colors)

        # title
        ax.set_title(f'KNN for: {result.query_taxon}')

        # Add the line of best fit for the base model
        base_grad = result.linear_model_base.params.gradient
        base_int = result.linear_model_base.params.intercept
        lobf_y = [base_grad * x + base_int for x in x_data]
        ax.plot(x_data, lobf_y, color='red', label='Base Model')

        plt.show()
        # plt.savefig('/tmp/am/embedded.png')
        # plt.close()

        print('next.')

    return


def load_tree_as_ete3(path_tree):
    with path_tree.open() as f:
        content = f.read()
    t = Tree(content, format=1, quoted_node_names=True)
    return t


def annotate_with_taxonomy(tree, d_taxonomy):
    # Add the taxonomy to the tree
    for leaf_node in tree.iter_leaves():
        cur_node = leaf_node.up
        while cur_node:
            descendants = [x.name for x in cur_node.get_descendants()]

            desc_tax = [d_taxonomy[x] for x in descendants]

            cur_node = cur_node.up
            print()

    return


def generate_data_for_export(qry: DPDTree, results, agg_results, d_taxonomy):
    # For each query taxon we want to find the maximal set of other query nodes that it will appear in

    return


def load_as_ete3(qry: DPDTree, results, agg_results, d_taxonomy):
    # Simple test, I want to know how many times each query taxon appears in a different knn set
    d_taxon_to_set = defaultdict(set)
    for result in results:
        for qry_gid in result.knn_as_linalg.qry_labels:
            d_taxon_to_set[qry_gid].update(result.knn_as_linalg.qry_labels)
    print()

    # So we can expect at minimum a given taxon to appear in knn sets with itself (e.g. 50).
    # The maximum varies a lot, for example ar53 gives 682 for the maximum. but drops sharply to 200

    def ete3_layout_fn(node):
        if node.is_leaf():
            data = agg_results[node.name]
            data_avg = np.mean(data)
            data_std = np.std(data)
            data_txt = f'{data_avg:.4f}+/-{data_std:.4f} {d_taxonomy[node.name]}'
            anno = faces.TextFace(data_txt, fsize=10)
            faces.add_face_to_node(anno, node, column=0)

    for result in sorted(results, key=lambda x: x.query_taxon):
        print(f'Query taxon: {result.query_taxon}')

        # Load the knn for the query taxon
        qry_labels = result.knn_as_linalg.qry_labels

        # Subset the dendropy tree to these nodes
        subtree = qry.subset_to_taxa(qry_labels)

        # Load it as an ete3 tree

        t = Tree(subtree.as_string(schema='newick'), format=1, quoted_node_names=True)
        ts = TreeStyle()
        ts.show_leaf_name = True

        # x = annotate_with_taxonomy(t, d_taxonomy)

        ts.layout_fn = ete3_layout_fn

        t.show(tree_style=ts)

        print()

    # t = load_tree_as_ete3(path_qry_tree)
    #
    # ts = TreeStyle()
    # ts.show_leaf_name = True
    #
    # t.show(tree_style=ts)

    return


def write_results(results, d_meta):
    # fig, ax = plt.subplots()

    # From the results we want to get the relative influence of each point
    d_taxon_relative_inf = defaultdict(list)
    d_taxon_data_points = defaultdict(list)
    # Get the relative influence of each point in the models that were created

    # Create plots of the points and the line of best fit
    for result in results:
        # taxon_removed = result.taxon_removed
        qry_labels = result.knn_as_linalg.qry_labels
        qry_data = result.knn_as_linalg.qry_data
        ref_labels = result.knn_as_linalg.ref_labels
        ref_data = result.knn_as_linalg.ref_data
        # relative_inf = result.relative_influence
        std_err = result.std_error

        for cur_qry, cur_ref, cur_qry_data, cur_ref_data, cur_inf in zip(qry_labels, ref_labels, qry_data, ref_data,
                                                                         std_err):
            d_taxon_relative_inf[cur_qry].append(cur_inf)
            d_taxon_data_points[cur_qry].append((cur_ref_data, cur_qry_data))
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

    return d_taxon_relative_inf, d_taxon_data_points


@app.command()
def goodbye(name: str, formal: bool = False):
    if formal:
        print(f"Goodbye Ms. {name}. Have a good day.")
    else:
        print(f"Bye {name}!")


if __name__ == '__main__':
    app()
