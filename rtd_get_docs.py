# this file
# is directly used by .readthedocs.yaml
# it extracts the built docs from the github artefact
# created by the build_docs.yml github action
import os
import pathlib
import zipfile

import requests


def download_and_extract_docs() -> None:
    token = os.environ.get("GITHUB_TOKEN")
    if token is None:
        token = os.environ.get("GITHUB_TOKEN_PRIVATE")

    headers = {"Authorization": f"token {token}"}
    api_url = "https://api.github.com/repos/iqtree/piqtree/actions/runs"
    response = requests.get(api_url, headers=headers, timeout=10)
    got = response.json()

    if got["status"] == "401":
        msg = "Request to GitHub API not authorised."
        raise requests.RequestException(msg)

    runs = got["workflow_runs"]
    latest_run = next(
        run
        for run in runs
        if run["name"] == "Build Docs" and run["conclusion"] == "success"
    )
    artifacts_url = latest_run["artifacts_url"]
    response = requests.get(artifacts_url, headers=headers, timeout=10)
    artifacts = response.json()["artifacts"]
    docs_artifact = next(
        artifact for artifact in artifacts if artifact["name"] == "piqtree-docs-html"
    )
    download_url = docs_artifact["archive_download_url"]
    response = requests.get(download_url, headers=headers, timeout=10)
    out = pathlib.Path("piqtree-docs-html.zip")
    out.write_bytes(response.content)
    with zipfile.ZipFile("piqtree-docs-html.zip", "r") as zip_ref:
        zip_ref.extractall("_readthedocs/html/")


if __name__ == "__main__":
    download_and_extract_docs()
