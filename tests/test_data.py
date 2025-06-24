import pathlib

from pytest_mock import MockerFixture

import piqtree


def test_dataset_names() -> None:
    names = piqtree.dataset_names()
    assert len(names) > 0


def test_download_dataset(tmp_path: pathlib.Path, mocker: MockerFixture) -> None:
    # Mock the requests response
    fake_content = b"fake dataset data"

    mock_response = mocker.Mock()
    mock_response.iter_content.return_value = [fake_content]
    mocker.patch("piqtree._data.requests.get", return_value=mock_response)

    dataset_name = "example.tree.gz"
    dest_path = piqtree.download_dataset(dataset_name, dest_dir=tmp_path)

    assert dest_path.exists()
    assert dest_path.read_bytes() == fake_content

    # Check correct arguments
    expected_url = piqtree._data._get_url(dataset_name)
    piqtree._data.requests.get.assert_called_once_with(  # type: ignore[attr-defined]
        expected_url,
        stream=True,
        timeout=20,
    )
