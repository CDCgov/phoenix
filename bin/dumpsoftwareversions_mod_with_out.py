#!/usr/bin/env python

import argparse
from pathlib import Path
import yaml
import platform
from textwrap import dedent

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--versions", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--workflow_name", required=False, default="Workflow")
parser.add_argument("--workflow_version", required=False, default="dev")
parser.add_argument("--nextflow_version", required=False, default="unknown")
parser.add_argument("--source", required=False, default="default")
args = parser.parse_args()

outdir = Path(args.outdir+"/"+args.source+"_pipeline_info")
#source_pipeline = args.source
outdir.mkdir(parents=True, exist_ok=True)

# HTML table generator
def _make_versions_html(versions):
    html = [
        dedent("""\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
        """)
    ]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(f"""\
                    <tr>
                        <td><samp>{process if i == 0 else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                """)
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\n".join(html)

# Record software versions for this module
versions_this_module = {
    "python": platform.python_version(),
    "yaml": yaml.__version__,
}
versions_dict = { "dumpsoftwareversions": versions_this_module }

# Load collated versions from input file
with open(args.versions) as f:
    versions_by_process = yaml.load(f, Loader=yaml.BaseLoader)
    versions_by_process = {**versions_by_process, **versions_dict}

# Aggregate by module
versions_by_module = {}
for process, process_versions in versions_by_process.items():
    module = process.split(":")[-1]
    try:
        assert versions_by_module[module] == process_versions, (
            "Inconsistent versions across modules. Open an issue with nf-core if this happens."
        )
    except KeyError:
        versions_by_module[module] = process_versions

# Add workflow version block
versions_by_module["Workflow"] = {
    "Nextflow": args.nextflow_version,
    args.workflow_name: args.workflow_version,
}

# Prepare MultiQC data block
versions_mqc = {
    "id": "software_versions",
    "section_name": f"{args.workflow_name} Software Versions",
    "section_href": f"https://github.com/{args.workflow_name}",
    "plot_type": "html",
    "description": "are collected at run time from the software output.",
    "data": _make_versions_html(versions_by_module),
}

if args.source == "default":
    filename="software_versions.yml"
else:
    filename=args.source+"_software_versions.yml"

# Write output files
with open(outdir / filename, "w") as f:
    yaml.dump(versions_by_module, f, default_flow_style=False)

#with open(outdir / "software_versions_mqc.yml", "w") as f:
#    yaml.dump(versions_mqc, f, default_flow_style=False)

#with open(outdir / "CENTAR_versions.yml", "w") as f:
#    yaml.dump(versions_dict, f, default_flow_style=False)