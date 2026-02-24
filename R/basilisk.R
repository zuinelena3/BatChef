#' @importFrom basilisk BasiliskEnvironment
py_env <- BasiliskEnvironment(
  envname = "py_env",
  pkgname = "BatChef",
  packages = c(
    "python=3.12.3",
    "scanpy=1.11.4",
    "scipy=1.16.2",
    "anndata=0.12.2",
    "scanorama=1.7.4",
    "leidenalg=0.10.2",
    "python-igraph=0.11.9",
    "bbknn=1.6.0",
    "scvi-tools=1.4.0"
  )
)
