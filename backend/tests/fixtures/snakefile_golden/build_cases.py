"""
Build golden Snakefile cases for official plugins (#137).

Copies each official plugin template (backend/plugin/official, plugin branch
release/plugins-v1.1) into ``templates/`` and writes ``cases.json`` with the
inputs the frontend sends (bare file names, metadata parameters). Golden outputs
are produced by ``generate_golden`` with the pre-#137 implementation and must
stay byte-identical afterwards.

Usage (from backend/): python tests/fixtures/snakefile_golden/build_cases.py
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
INPUT_FILES = {"input.h5ad": "pbmc_3k.h5ad", "geneList.txt": "genes.txt", "TF.txt": "human_TFs.txt"}
CONTEXT = {"user_name": "john", "workflow_id": "1", "algorithm_id": "2"}


def analysis_case(name, metadata):
    params = {}
    for rule in metadata["rules"].values():
        for p in rule.get("parameters", []):
            params[p["name"]] = H5AD_VALUES.get(p["name"], p.get("defaultValue"))
    user_input = {}
    for key, file_name in INPUT_FILES.items():
        source_dir = "tutorials" if key == "TF.txt" else "user/john/data"
        user_input[key] = f"{source_dir}/{file_name}"
    user_input.update(CONTEXT, plugin_name=name)
    return {"plugin": name, "kind": "analysis", "user_input": user_input, "plugin_params": params}


def visualization_cases(name, metadata):
    cases = []
    for rule in metadata["rules"].values():
        params, files = {}, {}
        for p in rule.get("parameters", []):
            if p["type"] in ("inputFile", "optionalInputFile"):
                files[p["defaultValue"]] = "Numlinks.sif" if p["defaultValue"] == "input.txt" else p["defaultValue"]
            elif p["type"] != "outputFile":
                params[p["name"]] = p.get("defaultValue")
        context = {"user_name": "john", "workflow_id": "1", "algorithm_id": "2", "visualization_id": "7",
                   "plugin_name": name, "visualization_name": rule["name"],
                   "result_filename": f"{rule['name']}_result.json"}
        cases.append({"plugin": name, "kind": "visualization", "rule": rule["name"],
                      "parameters": {**params, **context}, "file_mappings": files})
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
        os.makedirs(os.path.join(HERE, "templates", name), exist_ok=True)
        shutil.copy(os.path.join(plugin_dir, "Snakefile"), os.path.join(HERE, "templates", name, "Snakefile"))
        if name == "GRNViz":
            cases.extend(visualization_cases(name, metadata))
        else:
            cases.append(analysis_case(name, metadata))
    with open(os.path.join(HERE, "cases.json"), "w") as f:
        json.dump(cases, f, indent=2, ensure_ascii=False)
    print(f"{len(cases)} cases")


if __name__ == "__main__":
    sys.exit(main())
