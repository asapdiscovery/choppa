from click.testing import CliRunner
import traceback
import pytest
import json

from choppa.data.metadata.resources import NEXTSTRAIN_METADATA


def click_success(result):
    if result.exit_code != 0:  # -no-cov-  (only occurs on test error)
        print(result.output)
        traceback.print_tb(result.exc_info[2])
        print(result.exc_info[0], result.exc_info[1])
    return result.exit_code == 0


from choppa.cli.cli import cli
from choppa.data.toy_data.resources import (
    TOY_COMPLEX,
    TOY_FITNESS_DATA_COMPLETE,
    TOY_FITNESS_DATA_TRUNCATED,
)
from choppa.data.toy_data.resources import (
    TOY_FITNESS_DATA_SECTIONED,
    TOY_FITNESS_DATA_COMPLETE_NOCONF,
)


def test_choppa_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert click_success(result)


def test_choppa_cli_render(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "render",
            "--pdb-file",
            TOY_COMPLEX,
            "--fitness-file",
            TOY_FITNESS_DATA_SECTIONED,
            "--fitness-threshold",
            0.7,
            "--outfile-publication",
            tmp_path / "out.pse",
            "--outfile-interactive",
            tmp_path / "out.html",
        ],
    )
    assert click_success(result)
    # check the files exist
    assert (tmp_path / "out.pse").exists()
    assert (tmp_path / "out.html").exists()


def test_nextstrain_cli(tmp_path):
    """Runs through all available NextStrain viruses and pulls down mutation data for each. Will fail on e.g. a 404"""

    runner = CliRunner()

    with open(NEXTSTRAIN_METADATA) as f:
        d = json.load(f)
    for virus, query_data in d.items():
        gene = query_data["genes"][
            0
        ]  # just take the first gene for this virus, we're testing the connection
        result = runner.invoke(
            cli,
            [
                "nextstrain",
                "-v",
                virus,
                "-g",
                gene,
                "-o",
                f"{tmp_path}/nextstrain_output.csv",
            ],
        )
        assert result.exit_code == 0
