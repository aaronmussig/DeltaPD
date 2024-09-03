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

    // Ref total sum: 1235.005060000002
    // Qry total sum: 132.6799049183003
    let ref_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53_ref.dm");
    let qry_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53_qry.dm");
    let meta_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv");

    // let ref_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/example/reference.dm");
    // let qry_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/example/query.dm");
    // let meta_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/example/metadata.tsv");

    // let ref_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53.dm");
    // let qry_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200.dm");
    // let meta_path = PathBuf::from("/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200_metadata.tsv");

    println!("Loading query matrix...");
    let qry_matrix = DistMatrix::from_path(&qry_path).unwrap();

    println!("Loading reference matrix...");
    let ref_matrix = DistMatrix::from_path(&ref_path).unwrap();

    println!("Loading metadata...");
    let meta = MetadataFile::read(&meta_path, b'\t').unwrap();

    let deltapd = DeltaPD::new(
        QryDistMatrix { dm: qry_matrix},
        RefDistMatrix { dm: ref_matrix},
        meta
    );


    let params = Params::new(100, 1, LinearModelType::TheilSen, LinearModelError::RMSE, LinearModelCorr::R2);

    let _x = run_deltapd(&deltapd.qry_dm, &deltapd.ref_dm, &deltapd.metadata, &params).unwrap();


    println!("Done.");
}