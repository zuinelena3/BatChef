#' @importFrom basilisk BasiliskEnvironment
py_env <- BasiliskEnvironment(envname = "py_env", pkgname = "BatChef",
                              channels = c("conda-forge", "bioconda"),
                              packages = c("python==3.10.12",
                                           "scanorama==1.7.4",
                                           "scanpy==1.10.1",
                                           "python-igraph==0.11.4",
                                           "leidenalg==0.10.2",
                                           "scikit-learn==1.4.2",
                                           "anndata==0.10.7", "bbknn==1.6.0",
                                           "scvi-tools==1.1.2",
                                           "lightning==2.0.9.post0",
                                           "openssl==3.2.0"))
