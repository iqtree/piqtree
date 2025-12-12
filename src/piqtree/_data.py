import pathlib
import zipfile

import requests

_data_files = {
    "mammal-orths.zip": "https://github.com/user-attachments/files/17806562/mammal-orths.zip",
    "brca1.fasta.gz": "https://github.com/user-attachments/files/17806563/brca1.fasta.gz",
    "example.phy.gz": "https://github.com/user-attachments/files/17806561/example.phy.gz",
    "example.tree.gz": "https://github.com/user-attachments/files/17821150/example.tree.gz",
}


def _inflate_zip(zip_path: pathlib.Path, output_dir: pathlib.Path) -> pathlib.Path:
    """Decompress the contents of a zip file to a named directory."""
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(output_dir)
    return output_dir


def _get_url(name: str) -> str:
    """URL for a data file."""
    if name not in _data_files:
        msg = f"Unknown data file: {name}"
        raise ValueError(msg)
    return _data_files[name]


def dataset_names() -> list[str]:
    """Return the names of available datasets."""
    return list(_data_files.keys())


def download_dataset(
    name: str,
    dest_dir: str | pathlib.Path,
    dest_name: str | None = None,
    *,
    inflate_zip: bool = True,
) -> pathlib.Path:
    """Download data files used in docs, requires an internet connection.

    Parameters
    ----------
    name : str
        Data set name, see `dataset_names()`.
    dest_dir : str | pathlib.Path
        Where to write a local copy.
    dest_name : str | None, optional
        Name of the file to write, if None (default) uses name.
    inflate_zip : bool, optional
        Whether to unzip archives, by default True.

    Returns
    -------
    pathlib.Path
        Path to the downloaded file

    Notes
    -----
    Only downloads if dest_dir / dest_name does not exist.

    """
    dest_dir = pathlib.Path(dest_dir)

    url = _get_url(name)
    dest_name = dest_name or name
    outpath = dest_dir / dest_name
    outpath.parent.mkdir(parents=True, exist_ok=True)
    if outpath.exists():
        return outpath

    response = requests.get(url, stream=True, timeout=20)
    block_size = 4096
    with outpath.open("wb") as out:
        for data in response.iter_content(block_size):
            out.write(data)

    if inflate_zip and outpath.suffix == ".zip":
        outpath = _inflate_zip(outpath, dest_dir)
    return outpath
