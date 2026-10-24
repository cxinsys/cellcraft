"""
Build golden Snakefile cases for official plugins (#137).

Copies each official plugin template and the rules of its metadata.json
(backend/plugin/official, plugin branch release/plugins-v1.1) into
``templates/`` and writes ``cases.json`` with requests shaped like the
frontend sends them (Drawflow workflow_info / visualization request).

Usage (from backend/, plugin submodule on release/plugins-v1.1):
    python tests/fixtures/snakefile_golden/build_cases.py
    CELLCRAFT_BACKEND=<backend of the reference revision> \\
        python tests/fixtures/snakefile_golden/generate_golden.py
"""
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BACKEND = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OFFICIAL = os.path.join(BACKEND, "plugin", "official")

# 프론트엔드가 실제로 보내는 형태의 값 (FE 분석 기준)
H5AD_VALUES = {
    "UMAP lasso": "./user/john/workflow_1/ScatterPlot_4.json",
    "cell group": "leiden",
    "pseudotime": "dpt_pseudotime",
    "clusters": ["0", "1", "B cell"],
}
INPUT_FILES = {  # defaultValue -> (file_name, fileSource)
    "input.h5ad": ("pbmc_3k.h5ad", "user"),
    "geneList.txt": ("genes.txt", "user"),
    "TF.txt": ("human_TFs.txt", "shared"),
}


def analysis_case(name, metadata):
    io, rules = [], []
    for rule in metadata["rules"].values():
        params = []
        for p in rule.get("parameters", []):
            if p["type"] in ("inputFile", "optionalInputFile"):
                file_name, source = INPUT_FILES.get(p["defaultValue"], (None, "user"))
                io.append({"type": p["type"], "defaultValue": p["defaultValue"],
                           "file_name": file_name, "fileSource": source})
            elif p["type"] == "outputFile":
                io.append({"type": "outputFile", "defaultValue": p["defaultValue"], "activate": True})
            params.append({"name": p["name"], "type": p["type"],
                           "defaultValue": H5AD_VALUES.get(p["name"], p.get("defaultValue"))})
        rules.append({"name": rule["name"], "parameters": params})
    node = {"id": 2, "name": "algorithm", "class": "Algorithm", "outputs": {},
            "data": {"selectedPlugin": {"name": name, "source": "official"},
                     "selectedPluginInputOutput": io, "selectedPluginRules": rules}}
    return {"id": name, "plugin": name, "kind": "analysis",
            "workflow_info": {"drawflow": {"Home": {"data": {"2": node}}}}}


def visualization_cases(name, metadata):
    cases = []
    for rule in metadata["rules"].values():
        params = []
        for p in rule.get("parameters", []):
            param = {"type": p["type"], "name": p["name"], "defaultValue": p.get("defaultValue")}
            if p["type"] in ("inputFile", "optionalInputFile"):
                param["selectedFile"] = "Numlinks.sif" if p["defaultValue"] == "input.txt" else p["defaultValue"]
            params.append(param)
        cases.append({"id": f"{name}_{rule['name'].replace(' ', '_')}", "plugin": name, "kind": "visualization",
                      "request": {"id": 1, "current_node_id": "7", "algorithm_id": "2",
                                  "selectedVisualizationPlugin": {"name": name, "source": "official"},
                                  "selectedScript": {"name": rule["name"]},
                                  "selectedVisualizationParams": params}})
    return cases


def main():
    cases = []
    for name in sorted(os.listdir(OFFICIAL)):
        plugin_dir = os.path.join(OFFICIAL, name)
        meta_path = os.path.join(plugin_dir, "metadata.json")
        if not os.path.isfile(meta_path):
            continue
        with open(meta_path) as f:
            metadata = json.load(f)
        target = os.path.join(HERE, "templates", name)
        os.makedirs(target, exist_ok=True)
        shutil.copy(os.path.join(plugin_dir, "Snakefile"), os.path.join(target, "Snakefile"))
        with open(os.path.join(target, "metadata.json"), "w") as f:
            json.dump({"name": metadata["name"], "rules": metadata["rules"]}, f, indent=2, ensure_ascii=False)
        if name == "GRNViz":
            cases.extend(visualization_cases(name, metadata))
        else:
            cases.append(analysis_case(name, metadata))
    with open(os.path.join(HERE, "cases.json"), "w") as f:
        json.dump(cases, f, indent=2, ensure_ascii=False)
    print(f"{len(cases)} cases")


if __name__ == "__main__":
    sys.exit(main())
