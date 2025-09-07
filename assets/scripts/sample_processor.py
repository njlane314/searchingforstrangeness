from __future__ import annotations

import concurrent.futures
import json
import logging
import re
import subprocess
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import uproot


logger = logging.getLogger(__name__)


def setup_logger() -> None:
    """Configure basic logging for the module."""
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


def get_xml_entities(xml_path: Path) -> dict:
    content = xml_path.read_text()
    entity_regex = re.compile(r"<!ENTITY\s+([^\s]+)\s+\"([^\"]+)\">")
    return {match.group(1): match.group(2) for match in entity_regex.finditer(content)}


def run_command(command: list[str], execute: bool) -> bool:
    logger.info("COMMAND %s", " ".join(command))
    if not execute:
        logger.info("Dry run mode. HADD command not executed.")
        return True

    try:
        subprocess.run(command, check=True)
        logger.info("HADD Execution successful.")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        logger.error("HADD Execution failed: %s", exc)
        return False


def get_total_pot_from_single_file(file_path: Path) -> float:
    if not file_path or not file_path.is_file():
        return 0.0
    try:
        with uproot.open(file_path) as root_file:
            if "nuselection/SubRun" in root_file:
                subrun_tree = root_file["nuselection/SubRun"]
                if "pot" in subrun_tree:
                    return float(subrun_tree["pot"].array(library="np").sum())
    except Exception as exc:  # pylint: disable=broad-except
        logger.warning("Could not read POT from %s: %s", file_path, exc)
    return 0.0


def get_total_triggers_from_single_file(file_path: Path) -> int:
    if not file_path or not file_path.is_file():
        return 0
    try:
        with uproot.open(file_path) as root_file:
            tree_path = "nuselection/EventSelectionFilter"
            if tree_path in root_file:
                tree = root_file[tree_path]
                if "software_trigger" in tree.keys():
                    data = tree["software_trigger"].array(library="np")
                    return int(data.sum())

                algo_branches = [
                    "software_trigger_post",
                    "software_trigger_pre",
                    "software_trigger_post_ext",
                    "software_trigger_pre_ext",
                ]

                arrays = [
                    tree[branch].array(library="np") != 0
                    for branch in algo_branches
                    if branch in tree.keys()
                ]
                if arrays:
                    combined = arrays[0]
                    for arr in arrays[1:]:
                        combined |= arr
                    return int(combined.sum())
    except Exception as exc:  # pylint: disable=broad-except
        logger.warning("Could not read triggers from %s: %s", file_path, exc)
    return 0


def get_total_pot_from_files_parallel(file_paths: list[str]) -> float:
    with concurrent.futures.ThreadPoolExecutor() as executor:
        return sum(executor.map(get_total_pot_from_single_file, map(Path, file_paths)))


def get_total_triggers_from_files_parallel(file_paths: list[str]) -> int:
    with concurrent.futures.ThreadPoolExecutor() as executor:
        return sum(
            executor.map(get_total_triggers_from_single_file, map(Path, file_paths))
        )


def resolve_input_dir(
    stage_name: str | None, stage_outdirs: dict, entities: dict
) -> str | None:
    if not stage_name or stage_name not in stage_outdirs:
        return None
    input_dir = stage_outdirs[stage_name]
    for name, value in entities.items():
        input_dir = input_dir.replace(f"&{name};", value)
    return input_dir


def process_sample_entry(
    entry: dict,
    processed_analysis_path: Path,
    stage_outdirs: dict,
    entities: dict,
    nominal_pot: float,
    is_detvar: bool = False,
) -> bool:
    execute_hadd = entry.pop("do_hadd", False)

    if not entry.get("active", True):
        sample_key = entry.get("sample_key", "UNKNOWN")
        logger.info(
            "Skipping %s: %s (marked as inactive)",
            "detector variation" if is_detvar else "sample",
            sample_key,
        )
        return False

    stage_name = entry.get("stage_name")
    sample_key = entry.get("sample_key")
    sample_type = entry.get("sample_type", "mc")

    logger.info(
        "Processing %s: %s (from stage: %s)",
        "detector variation" if is_detvar else "sample",
        sample_key,
        stage_name,
    )
    logger.info(
        "HADD execution for this %s: %s",
        "sample" if not is_detvar else "detector variation",
        "Enabled" if execute_hadd else "Disabled",
    )

    input_dir = resolve_input_dir(stage_name, stage_outdirs, entities)
    if not input_dir:
        logger.warning(
            "Stage '%s' not found in XML outdirs. Skipping %s '%s'.",
            stage_name,
            "detector variation" if is_detvar else "sample",
            sample_key,
        )
        return False

    output_file = processed_analysis_path / f"{sample_key}.root"
    entry["relative_path"] = output_file.name

    root_files = sorted(str(f) for f in Path(input_dir).rglob("*.root"))
    if not root_files:
        logger.warning("No ROOT files found in %s. HADD will be skipped.", input_dir)

    if root_files and execute_hadd:
        if not run_command(["hadd", "-f", str(output_file), *root_files], True):
            logger.error(
                "HADD failed for %s. Skipping further processing for this entry.",
                sample_key,
            )
            return False
    elif not root_files and execute_hadd:
        logger.info(
            "No ROOT files found for '%s'. Skipping HADD but proceeding to record metadata (if applicable).",
            sample_key,
        )
    elif root_files and not execute_hadd:
        logger.info("HADD not requested for '%s'. Skipping HADD command.", sample_key)

    source_files = [str(output_file)] if execute_hadd else root_files

    if sample_type == "mc" or is_detvar:
        entry["pot"] = get_total_pot_from_files_parallel(source_files)
        entry["triggers"] = get_total_triggers_from_files_parallel(source_files)
        if entry["pot"] == 0.0:
            entry["pot"] = nominal_pot
    elif sample_type == "data":
        entry["pot"] = nominal_pot
        entry["triggers"] = get_total_triggers_from_files_parallel(source_files)

    entry.pop("stage_name", None)
    return True


def main() -> None:
    setup_logger()
    definitions_path = Path("config/data.json")
    xml_paths = [
        "/exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs/ubana/ubana/searchingforstrangeness/xml/numi_fhc_workflow.xml"
    ]
    config_path = "config/samples.json"
    runs_process = ["run1"]

    logger.info("Loading configurations")
    with open(definitions_path) as f:
        config = json.load(f)

    entities: dict[str, str] = {}
    stage_outdirs: dict[str, str] = {}
    for xml in xml_paths:
        xml_path = Path(xml)
        entities.update(get_xml_entities(xml_path))
        tree = ET.parse(xml_path)
        project_node = tree.find("project")
        if project_node is None:
            logger.error("Could not find <project> tag in XML file '%s'.", xml)
            continue
        stage_outdirs.update(
            {
                s.get("name"): s.find("outdir").text
                for s in project_node.findall("stage")
                if s.find("outdir") is not None
            }
        )

    logger.info("Hadding files and calculating metadata")
    processed_analysis_path = Path(config["ntuple_base_directory"])
    processed_analysis_path.mkdir(parents=True, exist_ok=True)

    for beam, run_configs in config.get("run_configurations", {}).items():
        for run, run_details in run_configs.items():
            if runs_process and run not in runs_process:
                logger.info("Skipping run: %s (not in specified RUNS_PROCESS list)", run)
                continue

            logger.info("Processing run: %s", run)
            nominal_pot = run_details.get("nominal_pot", 0.0)

            if nominal_pot == 0.0:
                logger.warning(
                    "No nominal_pot specified for run '%s'. MC scaling might be incorrect.",
                    run,
                )
            else:
                logger.info("Using nominal POT for this run: %.4e", nominal_pot)

            for sample in run_details.get("samples", []):
                if process_sample_entry(
                    sample,
                    processed_analysis_path,
                    stage_outdirs,
                    entities,
                    nominal_pot,
                ):
                    if "detector_variations" in sample:
                        for detvar_sample in sample["detector_variations"]:
                            process_sample_entry(
                                detvar_sample,
                                processed_analysis_path,
                                stage_outdirs,
                                entities,
                                nominal_pot,
                                is_detvar=True,
                            )

    output_path = Path(config_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    samples_cfg = {
        "samples": {
            "ntuple_directory": config["ntuple_base_directory"],
            "beamlines": config["run_configurations"],
        }
    }
    with open(output_path, "w") as f:
        json.dump(samples_cfg, f, indent=4)

    logger.info("Workflow complete")
    logger.info("Successfully generated configuration at '%s'", output_path)


if __name__ == "__main__":
    main()

