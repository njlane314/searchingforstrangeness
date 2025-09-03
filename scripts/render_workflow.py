import json
import re
from string import Template
from pathlib import Path


def load_entity_mapping(template_path: Path) -> dict:
    content = template_path.read_text()
    # extract simple entity definitions
    pairs = re.findall(r'<!ENTITY\s+(\w+)\s+"([^"\n]*)">', content)
    mapping = {value: f"&{name};" for name, value in pairs}
    # jobsub_config spans multiple lines; handle explicitly
    mapping["--tar_file_name=file:///pnfs/uboone/scratch/users/nlane/tarballs/strangeness_assets.tar.gz"] = "&jobsub_config;"
    # apply longer replacements first
    return dict(sorted(mapping.items(), key=lambda kv: len(kv[0]), reverse=True))


def apply_entities(text: str, mapping: dict) -> str:
    for value, entity in mapping.items():
        text = text.replace(value, entity)
    return text


def render_stage(stage: dict, mapping: dict, template: Template) -> str:
    lines = []
    if 'inputdef' in stage:
        lines.append(f"    <inputdef>{apply_entities(stage['inputdef'], mapping)}</inputdef>")
    if 'inputstage' in stage:
        lines.append(f"    <inputstage>{apply_entities(stage['inputstage'], mapping)}</inputstage>")
    for fcl in stage.get('fcl', []):
        lines.append(f"    <fcl>{apply_entities(fcl, mapping)}</fcl>")
    uses_selection = any('run_neutrinoselection.fcl' in fcl for fcl in stage.get('fcl', []))
    if uses_selection:
        lines.append("    <initsource>&initsrc;</initsource>")
    elif 'initsource' in stage:
        lines.append(f"    <initsource>{apply_entities(stage['initsource'], mapping)}</initsource>")
    lines.append(f"    <outdir>{apply_entities(stage['outdir'], mapping)}</outdir>")
    lines.append(f"    <logdir>{apply_entities(stage['logdir'], mapping)}</logdir>")
    lines.append(f"    <workdir>{apply_entities(stage['workdir'], mapping)}</workdir>")
    lines.append(f"    <datatier>{apply_entities(stage['datatier'], mapping)}</datatier>")
    lines.append(f"    <numjobs>{apply_entities(str(stage['numjobs']), mapping)}</numjobs>")
    if uses_selection:
        lines.append("    <jobsub>&jobsub_config;</jobsub>")
    elif 'jobsub' in stage:
        jobsub_lines = '\n'.join('        ' + apply_entities(line, mapping) for line in stage['jobsub'])
        lines.append('    <jobsub>\n' + jobsub_lines + '\n    </jobsub>')
    content = '\n'.join(lines)
    return template.substitute(name=stage['name'], content=content)


def main() -> None:
    base = Path(__file__).resolve().parents[1]
    template_path = base / 'xml' / 'numi_fhc_workflow_template.xml'
    stage_template_path = base / 'xml' / 'stage_template.xml'
    stages_path = base / 'xml' / 'stages.json'

    mapping = load_entity_mapping(template_path)
    stage_template = Template(stage_template_path.read_text())
    stages = json.loads(stages_path.read_text())
    rendered_stages = [render_stage(stage, mapping, stage_template) for stage in stages]
    stages_str = '\n\n'.join(rendered_stages)

    workflow_template = Template(template_path.read_text())
    output = workflow_template.substitute(stages=stages_str)

    (base / 'xml' / 'numi_fhc_workflow.xml').write_text(output)


if __name__ == '__main__':
    main()
