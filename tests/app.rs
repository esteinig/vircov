use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec!["--alignment", "file/doesnt/exist.paf"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("does not exist"));

    Ok(())
}

#[test]
fn valid_inputs_raise_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "-v",
    ]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn valid_output_string_default_filters() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(
        "21172389_LCMV_L-segment_final\t3\t2\t3\t321\t7194\t0.0446\t-\t-\n21172389_LCMV_S-segment_final\t2\t2\t2\t137\t3407\t0.0402\t-\t-\n"
    ));

    Ok(())
}

#[test]
fn valid_output_string_default_filters_verbose_one_no_descr(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "-v",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(
        "21172389_LCMV_L-segment_final\t3\t2\t3\t321\t7194\t0.0446\t-\t-",
    ));

    Ok(())
}

#[test]
fn valid_output_string_default_filters_region_threshold_pass(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "--regions",
        "2",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(
        "21172389_LCMV_L-segment_final\t3\t2\t3\t321\t7194\t0.0446\t-\t-",
    ));

    Ok(())
}

#[test]
fn valid_output_string_default_filters_region_threshold_none(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "--regions",
        "3",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(""));

    Ok(())
}

#[test]
fn valid_output_string_default_filters_refseq_length_pass() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "--length",
        "5000",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(
        "21172389_LCMV_L-segment_final\t3\t2\t3\t321\t7194\t0.0446\t-\t-",
    ));

    Ok(())
}

#[test]
fn valid_output_string_default_filters_refseq_length_none() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "--alignment",
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_paf_ok.fasta",
        "--length",
        "8000",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(""));

    Ok(())
}
