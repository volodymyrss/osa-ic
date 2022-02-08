from click.testing import CliRunner
from osaic.integralicindex import cli

def test_create():
    cli_runner = CliRunner()
    cli_runner.invoke(cli, ["create"])
    

def test_list():
    cli_runner = CliRunner()
    cli_runner.invoke(cli, ["list"])
