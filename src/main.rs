use std::path::PathBuf;
use rs_deltapd::method::common_taxa::find_common_ids_in_pdms;
use rs_deltapd::method::deltapd::run_deltapd;
use rs_deltapd::method::knn::{create_vecs_from_knn2, get_pdm_k_nearest_neighbours_from_matrix};
use rs_deltapd::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType};
use rs_deltapd::model::metadata::MetadataFile;
use rs_deltapd::model::params::Params;
use rs_deltapd::model::pdm::{DistMatrix, QryDistMatrix, RefDistMatrix};
use rs_deltapd::python::deltapd::DeltaPD;

fn main() {



    let knn = 3;

    let ref_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.dm");
    let qry_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.dm");
    let meta_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv");

    let qry_matrix = DistMatrix::from_path(&qry_path).unwrap();
    let ref_matrix = DistMatrix::from_path(&ref_path).unwrap();

    let meta = MetadataFile::read(&meta_path, b'\t').unwrap();

    let deltapd = DeltaPD::new(
        QryDistMatrix { dm: qry_matrix},
        RefDistMatrix { dm: ref_matrix},
        meta
    );

    // // Identify the common taxa between both trees
    // let common_ids = find_common_ids_in_pdms(&deltapd.ref_dm.dm.taxa, &deltapd.qry_dm.dm.taxa, &deltapd.metadata).unwrap();
    //
    // // Get the knn
    // let knn_qry = get_pdm_k_nearest_neighbours_from_matrix(
    //     &deltapd.qry_dm.dm.taxa,
    //     &common_ids.qry_taxa,
    //     &deltapd.qry_dm.dm.matrix,
    //     knn,
    // ).unwrap();
    //
    // let knn_ref = get_pdm_k_nearest_neighbours_from_matrix(
    //     &deltapd.ref_dm.dm.taxa,
    //     &common_ids.ref_taxa,
    //     &deltapd.ref_dm.dm.matrix,
    //     knn,
    // ).unwrap();

    let query_taxon = &deltapd.qry_dm.dm.taxa[0];

    let params = Params::new(4, 1, LinearModelType::RepeatedMedian, LinearModelError::RMSE, LinearModelCorr::R2);

    let x = run_deltapd(&deltapd.qry_dm, &deltapd.ref_dm, &deltapd.metadata, &params).unwrap();

    // let x = create_vecs_from_knn2(&deltapd.qry_dm, &deltapd.ref_dm, query_taxon, &deltapd.metadata, knn_qry, knn_ref).unwrap();


    println!("Done.");
}