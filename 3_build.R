setwd("~/git/snfa-scripts")

devtools::build("../snfa", ".")
system("R CMD Rd2pdf ../snfa --force --output=snfa.pdf --no-preview .")
