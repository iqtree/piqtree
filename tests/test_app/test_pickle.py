"""Apps must be pickleable to be able to be run with parallel=True"""

import pickle

from cogent3 import get_app, make_tree

from piqtree._app import _ALL_APP_NAMES


def test_pickle() -> None:
    app_args = {
        "piq_build_tree": ("JC",),
        "piq_fit_tree": (make_tree("(a,b,(c,d));"), "JC"),
    }
    for app_name in _ALL_APP_NAMES:
        app = get_app(app_name, *app_args.get(app_name, ()))
        assert len(pickle.dumps(app)) > 0
