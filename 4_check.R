setwd("~/git/snfa-devel")
devtools::check("../snfa")

# setwd("~/git/snfa")
# usethis::use_travis()

setwd("~/git/snfa-devel")
system("R CMD check --as-cran snfa_0.0.1.tar.gz")
