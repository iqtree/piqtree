### ENH

- Added `fit_tree_hessian` which fits a tree and additionally returns the
  log-likelihood gradient (first derivatives w.r.t. branch lengths) and the
  full Hessian matrix at the fitted point. Results are captured in memory via
  a new `fit_tree_hessian` entry point in the IQ-TREE library — no
  `.mcmctree.hessian` file is written.
