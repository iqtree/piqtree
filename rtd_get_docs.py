# this file
# is directly used by .readthedocs.yaml
# it extracts the built docs from the github artefact
# created by the build_docs.yml github action
import os
import pathlib
import time
import zipfile

import requests

MAX_WAIT_TIME = 600  # Typically takes 5.5 minutes
POLL_INTERVAL = 120


def get_rtd_version_name() -> str:
    return os.environ.get("READTHEDOCS_VERSION_NAME", "").lower()


def get_github_token() -> str:
    token = os.environ.get("GITHUB_TOKEN")
    if token is None:
        token = os.environ.get("GITHUB_TOKEN_PRIVATE")
    if not token:
        msg = "GitHub token not found."
        raise OSError(msg)
    return token


def get_latest_run(workflow_filename: str, headers: dict) -> dict:
    url = f"https://api.github.com/repos/iqtree/piqtree/actions/workflows/{workflow_filename}/runs"

    response = requests.get(url, headers=headers, timeout=10)
    response.raise_for_status()

    runs = response.json()["workflow_runs"]
    if not runs:
        msg = f"No workflow runs found for: '{workflow_filename}'"
        raise ValueError(msg)

    return runs[0]


def wait_for_run_completion(run: dict, headers: dict) -> dict:
    run_id = run["id"]
    run_url = f"https://api.github.com/repos/iqtree/piqtree/actions/runs/{run_id}"

    waited = 0
    while waited < MAX_WAIT_TIME:
        response = requests.get(run_url, headers=headers, timeout=10)
        response.raise_for_status()
        run_status = response.json()

        status = run_status["status"]
        conclusion = run_status["conclusion"]

        if status == "completed":
            if conclusion != "success":
                msg = f"Latest workflow run failed with conclusion: '{conclusion}'"
                raise RuntimeError(msg)
            return run_status

        time.sleep(POLL_INTERVAL)
        waited += POLL_INTERVAL

    msg = "Timed out waiting for workflow run to complete."
    raise TimeoutError(msg)


def download_and_extract_artifact(run: dict, headers: dict) -> None:
    artifact_name = "piqtree-docs-html"
    artifacts_url = run["artifacts_url"]

    response = requests.get(artifacts_url, headers=headers, timeout=10)
    response.raise_for_status()

    artifacts = response.json()["artifacts"]

    artifact = next((a for a in artifacts if a["name"] == artifact_name), None)
    if artifact is None:
        msg = f"Artifact '{artifact_name}' not found in the run."
        raise ValueError(msg)

    download_url = artifact["archive_download_url"]
    response = requests.get(download_url, headers=headers, timeout=10)
    response.raise_for_status()

    out = pathlib.Path(f"{artifact_name}.zip")
    out.write_bytes(response.content)

    with zipfile.ZipFile(out, "r") as zip_ref:
        zip_ref.extractall("_readthedocs/html/")


def download_and_extract_docs() -> None:
    version = get_rtd_version_name()

    if version not in ("latest", "stable"):
        msg = f"Unexpected version '{version}' for readthedocs."
        raise ValueError(msg)

    workflow_filename = "release.yml" if version == "stable" else "build_docs.yml"

    headers = {"Authorization": f"token {get_github_token()}"}

    latest_run = get_latest_run(workflow_filename, headers)
    completed_run = wait_for_run_completion(latest_run, headers)
    download_and_extract_artifact(completed_run, headers)


if __name__ == "__main__":
    download_and_extract_docs()
