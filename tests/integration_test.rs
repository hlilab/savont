use assert_cmd::Command;
use minimap2;
use std::fs;
use tempfile::TempDir;

const REF_FA: &str = "tests/data/zymo_ref_asvs.fa.gz";
const READS_FQ: &str = "tests/data/ont_zymo_1000.fq.gz";

/// Run `savont asv` on the bundled 1000-read zymo dataset, then verify that
/// every output ASV aligns to the zymo reference ASVs with zero mismatches,
/// using the minimap2 Rust bindings (no external binary required).
#[test]
fn test_asv_generation_and_perfect_alignment() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().to_str().unwrap();

    // --- 1. run savont asv ---
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", READS_FQ,
            "-o", out_dir,
            "-t", "4",
            // 1000 reads → depth per cluster can be modest; loosen size gate
            "--min-cluster-size", "5",
        ])
        .assert()
        .success();

    let asv_fasta = tmp.path().join("final_asvs.fasta");
    assert!(asv_fasta.exists(), "final_asvs.fasta was not created");

    // --- 2. load ASV sequences ---
    let asvs = savont::taxonomy::load_fasta_with_needletail(&asv_fasta)
        .expect("failed to read final_asvs.fasta");
    assert!(!asvs.is_empty(), "savont produced zero ASVs");

    // --- 3. build minimap2 index from reference, align each ASV ---
    let aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar()
        .with_index(REF_FA, None)
        .expect("failed to build minimap2 index from zymo reference");

    let mut mapped = 0usize;
    let mut imperfect: Vec<(String, i32)> = Vec::new();

    for (header, seq) in &asvs {
        let hits = aligner
            .map(seq, true, false, None, None, None)
            .expect("minimap2 mapping failed");

        // Keep only the primary hit (lowest NM / highest score comes first)
        let primary = hits
            .iter()
            .find(|h| h.is_primary);

        match primary {
            None => { /* unmapped – will be caught by the mapped == asvs.len() assert */ }
            Some(hit) => {
                mapped += 1;
                let nm = hit.alignment.as_ref().map(|a| a.nm).unwrap_or(i32::MAX);
                if nm > 0 {
                    imperfect.push((header.clone(), nm));
                }
            }
        }
    }

    assert!(mapped > 0, "no ASVs aligned to the zymo reference");
    assert_eq!(
        mapped,
        asvs.len(),
        "only {}/{} ASVs mapped to the zymo reference",
        mapped,
        asvs.len()
    );
    assert!(
        imperfect.is_empty(),
        "ASVs with NM > 0 (not perfect): {:?}",
        imperfect
    );
}

// ── database download + load tests ──────────────────────────────────────────
// These are skipped by default (`cargo test`) because they download large
// files from the internet.  Run them explicitly with:
//   cargo test -- --ignored

/// Download the EMU database, check the expected files are present, and load
/// the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_emu_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "emu"])
        .assert()
        .success();

    let db_dir = tmp.path().join("emu");
    assert!(db_dir.join("species_taxid.fasta").exists(), "species_taxid.fasta missing");
    assert!(db_dir.join("taxonomy.tsv").exists(), "taxonomy.tsv missing");
    assert!(db_dir.join(".savont_db").exists(), "marker file missing");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load EMU database via registry");
    assert!(!db.taxonomy.is_empty(), "EMU taxonomy map is empty");
    assert!(db.fasta_path.exists(), "FASTA path does not exist on disk");
}

/// Download the SILVA database, check the expected files are present, and load
/// the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_silva_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "silva-138.2"])
        .assert()
        .success();

    let db_dir = tmp.path().join("silva-138.2");
    let has_fasta = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.ends_with(".fasta.gz") || n.ends_with(".fasta")));
    assert!(has_fasta, "no FASTA file found after download");

    let has_taxmap = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.starts_with("taxmap_") && n.ends_with(".txt")));
    assert!(has_taxmap, "no taxmap file found after download");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load SILVA database via registry");
    assert!(!db.taxonomy.is_empty(), "SILVA taxonomy map is empty");
}

/// Download the GTDB r232 SSU database, check the expected file is present,
/// and load the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_gtdb_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "gtdb-r232"])
        .assert()
        .success();

    let db_dir = tmp.path().join("gtdb-r232");
    let has_fna = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str().map_or(false, |n| n.ends_with(".fna.gz")));
    assert!(has_fna, "no .fna.gz file found after download");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load GTDB database via registry");
    assert!(!db.taxonomy.is_empty(), "GTDB taxonomy map is empty");

    let bad: Vec<_> = db.taxonomy.values()
        .filter(|e| e.superkingdom.is_empty())
        .take(5).map(|e| e.tax_id.clone()).collect();
    assert!(bad.is_empty(), "GTDB entries missing superkingdom: {:?}", bad);
}

/// Parse a small hand-crafted GTDB FASTA header without hitting the network,
/// to confirm the taxonomy parser is correct.
#[test]
fn test_gtdb_taxonomy_parser_unit() {
    use std::io::Write;

    let tmp = TempDir::new().unwrap();

    // Write a two-entry mock GTDB FASTA into a temp directory
    let fna_path = tmp.path().join("mock_gtdb.fna");
    {
        let mut f = fs::File::create(&fna_path).unwrap();
        writeln!(
            f,
            ">RS_GCF_000001405.40~NC_000001.11 d__Bacteria;p__Pseudomonadota;\
c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;\
g__Escherichia;s__Escherichia coli [location=1..1500] [ssu_len=1500]"
        )
        .unwrap();
        writeln!(f, "ACGT").unwrap();
        writeln!(
            f,
            ">GB_GCA_000007185.1~AE017221.1 d__Archaea;p__Thermoproteota;\
c__Thermoprotei;o__Thermoproteales;f__Thermoproteaceae;\
g__Thermoproteus;s__Thermoproteus tenax [location=1..1200] [ssu_len=1200]"
        )
        .unwrap();
        writeln!(f, "TTTT").unwrap();
    }

    let db = savont::taxonomy::Database::load_gtdb(tmp.path())
        .expect("load_gtdb failed on mock file");

    assert_eq!(db.taxonomy.len(), 2);

    let ecoli = db
        .taxonomy
        .get("RS_GCF_000001405.40~NC_000001.11")
        .expect("E. coli entry missing");
    assert_eq!(ecoli.superkingdom, "Bacteria");
    assert_eq!(ecoli.phylum, "Pseudomonadota");
    assert_eq!(ecoli.class, "Gammaproteobacteria");
    assert_eq!(ecoli.order, "Enterobacterales");
    assert_eq!(ecoli.family, "Enterobacteriaceae");
    assert_eq!(ecoli.genus, "Escherichia");
    assert_eq!(ecoli.species, "Escherichia coli");

    let archaea = db
        .taxonomy
        .get("GB_GCA_000007185.1~AE017221.1")
        .expect("Archaea entry missing");
    assert_eq!(archaea.superkingdom, "Archaea");
    assert_eq!(archaea.genus, "Thermoproteus");
    assert_eq!(archaea.species, "Thermoproteus tenax");
}
