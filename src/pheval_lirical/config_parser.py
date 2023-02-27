from dataclasses import dataclass
from pathlib import Path

import yaml
from serde import serde
from serde.yaml import from_yaml


@serde
@dataclass
class ExomiserLiricalExomiserConfigs:
    path_to_exomiser_data_directory: Path


@serde
@dataclass
class LiricalConfigRun:
    environment: str
    phenotype_only: bool
    version: str
    path_to_lirical_software_directory: Path
    exomiser_configurations: ExomiserLiricalExomiserConfigs


@serde
@dataclass
class LiricalPostProcess:
    sort_order: str


@serde
@dataclass
class LiricalConfig:
    run: LiricalConfigRun
    post_process: LiricalPostProcess


def parse_lirical_config(config_path: Path) -> LiricalConfig:
    """Parse the config file."""
    with open(config_path, "r") as config_file:
        config = yaml.safe_load(config_file)
    config_file.close()
    return from_yaml(LiricalConfig, yaml.dump(config))
