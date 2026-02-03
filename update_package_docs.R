getwd()
#setwd(dirname(getwd()))

# Create or update package metadata
usethis::use_description(fields = list(
  Package = "R4EwE",
  Title = "R for Ecopath with Ecosim (R4EwE)",
  Description = "Utilities, workflows, and helpers for Ecopath with Ecosim users.",
  "Authors@R" = utils::person("David", "Chagaris", email = "dchagaris@ufl.edu", role = c("aut", "cre")),
  License = "MIT + file LICENSE",
  Language = "en-US",
  Encoding = "UTF-8",
  Roxygen = "list(markdown = TRUE)"
))
usethis::use_mit_license("David Chagaris")
usethis::use_readme_rmd()

usethis::use_package("fields")
usethis::use_package("colorRamps")
usethis::use_package("ncdf4")
usethis::use_package("terra")
usethis::use_package("raster")
usethis::use_package("rvest")
usethis::use_package("httr")
usethis::use_package("digest")
usethis::use_package("doParallel")
usethis::use_package("reshape2")
usethis::use_package("gifski")
usethis::use_package("maps")
usethis::use_package("viridis")
usethis::use_package("foreach")
usethis::use_package("stats")

devtools::document()
devtools::check(vignettes=FALSE)

devtools::build()

install.packages("../R4EwE_0.0.0.9000.tar.gz", repos=NULL, type='source')
library('R4EwE')
?fn.runEwE
?fn.pull_cefi

usethis::use_readme_rmd()
devtools::build_readme()
