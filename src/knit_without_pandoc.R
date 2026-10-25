# Used by this report's YAML knit function in VS Code or RStudio.
# Manual alternative: source('src/knit_without_pandoc.R');
# knit_without_pandoc('src/test_rd_identification.Rmd')
knit_without_pandoc <- function(inputFile, ...) {
  for (pkg in c('knitr', 'markdown')) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop('Install the required R package first: install.packages("', pkg, '")')
    }
  }
  inputFile <- normalizePath(inputFile, winslash = '/', mustWork = TRUE)
  root <- dirname(inputFile)
  while (!file.exists(file.path(root, 'data', 'nhis_clean_full.csv'))) {
    parent <- dirname(root)
    if (identical(parent, root)) stop('Cannot locate project data.')
    root <- parent
  }
  out <- file.path(root, 'output', 'rd_identification')
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  old_opts <- knitr::opts_knit$get()
  old_chunk <- knitr::opts_chunk$get()
  on.exit({setwd(old_wd); knitr::opts_knit$restore(old_opts)
           knitr::opts_chunk$restore(old_chunk)}, add = TRUE)
  setwd(out)
  knitr::opts_knit$set(root.dir = out)
  stem <- tools::file_path_sans_ext(basename(inputFile))
  md <- tempfile(pattern=paste0(stem, '-'), tmpdir=out, fileext='.md')
  on.exit(if(file.exists(md)) file.remove(md), add=TRUE)
  # Embed plots in the HTML; figure intermediates stay in R's temporary folder.
  knitr::opts_chunk$set(fig.path=paste0(gsub('\\\\','/',tempdir()),'/',stem,'-'))
  html <- file.path(out, paste0(stem, '.html'))
  knitr::knit(inputFile, output = md, envir = new.env(parent = globalenv()),
              encoding = 'UTF-8')
  markdown::markdownToHTML(md, output = html,
    options = c('toc', 'base64_images'))
  message('Output created: ', html)
  invisible(html)
}
