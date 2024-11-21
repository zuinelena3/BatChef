#' Yaml file convertion
#'
#' @param path yml name file
#'
#' @importFrom yaml read_yaml
#'
yml_convert <- function(path){
  yaml_txt <- yaml::read_yaml(file = path)

  channels <- yaml_txt$channels

  deps <- unlist(yaml_txt$dependencies)

  python <- deps[startsWith(deps,"python")]
  python[1] <- gsub(paste0("=", sub(".*=", "", python[1])), "", python[1])
  python <- unname(python)

  deps <- deps[!startsWith(deps, "python")]
  pip_pkgs <- unname(deps[startsWith(names(deps),"pip")])
  nonpip_pkgs <- unname(deps[!startsWith(names(deps),"pip")])

  ll <- list(channels = channels, pkg = c(python, nonpip_pkgs), pip = pip_pkgs)
  return(ll)
}

#' @importFrom basilisk BasiliskEnvironment
#'
scvi_env <- BasiliskEnvironment(envname = "scvi_env", pkgname = "BatChef",
                           channels = yml_convert("inst/scvi_env.yml")$channels,
                           packages = yml_convert("inst/scvi_env.yml")$pkg,
                           pip = yml_convert("inst/scvi_env.yml")$pip)
