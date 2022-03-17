use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec!["file/doesnt/exist.paf"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("does not exist"));

    Ok(())
}

#[test]
fn valid_inputs_raise_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
    cmd.args(vec![
        "tests/cases/test_ok.paf",
        "--fasta",
        "tests/cases/test_ok.fasta",
        "-v",
    ]);

    cmd.assert().success();

    Ok(())
}
