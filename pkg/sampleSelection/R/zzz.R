.onAttach <- function( libname, pkgname ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'sampleSelection' package as:\n",
         "Toomet, Ott and Henningsen, Arne (2008). ",
        "Sample Selection Models in R: Package sampleSelection. ",
        "Journal of Statistical Software 27(7). ",
        "DOI 10.18637/jss.v027.i07.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'sampleSelection' package, ",
         "please use its R-Forge site: ",
         "https://r-forge.r-project.org/projects/sampleselection/" ),
      domain = NULL,  appendLF = TRUE )
}
