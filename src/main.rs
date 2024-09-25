use rs_deltapd::method::deltapd::run_deltapd;
use rs_deltapd::model::linalg::{LinearModelCorr, LinearModelError, LinearModelType};
use rs_deltapd::model::metadata::MetadataFile;
use rs_deltapd::model::params::Params;
use rs_deltapd::model::pdm::{DistMatrix, QryDistMatrix, RefDistMatrix};
use rs_deltapd::python::deltapd::DeltaPD;
use std::collections::HashSet;
use std::path::PathBuf;


fn main() {
    // Ref total sum: 1235.005060000002
    // Qry total sum: 132.6799049183003
    let ref_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.dm");
    let qry_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.dm");
    let meta_path = PathBuf::from("/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv");


    let out_dir = PathBuf::from("/tmp/out");
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
        QryDistMatrix { dm: qry_matrix },
        RefDistMatrix { dm: ref_matrix },
        meta,
    );

    /*
    include the distances to doubel cpoy genes
     */

    let taxa_subset: Vec<String> = vec![
        "GB_GCA_016185435.1".to_string(), // This one is a good long-branch false positive
        "RS_GCF_001307315.1".to_string(), // Within the correct genus
        "GB_GCA_012329085.1".to_string(),
        "GB_GCA_016219485.1".to_string(), // possible a true positive?
        "GB_GCA_018693315.1".to_string(),
    ];
    let taxa_subset: HashSet<String> = taxa_subset.into_iter().collect();

    let replicates = 3;
    let sample_size = 20.0; // if > 1 then this is the number of samples to use

    println!("Running with sample size: {}", sample_size);
    println!("Running with replicates: {}", replicates);

    let params = Params::new(1, sample_size, replicates, taxa_subset, LinearModelType::RepeatedMedian, LinearModelError::RMSE, LinearModelCorr::R2, true, out_dir);

    let _results = run_deltapd(&deltapd.qry_dm, &deltapd.ref_dm, &deltapd.metadata, &params).unwrap();


    println!("Done.");
}