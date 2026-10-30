# Regenerate docs/detailed_overview.md, the GitHub-readable version of
# docs/detailed_overview.Rmd. Run from the repository root after installing the
# package (R CMD INSTALL .). JAGS is needed for the backend-agreement sections.
#
# To read the document locally instead, knit docs/detailed_overview.Rmd to HTML.

render_md <- function() {
  rmarkdown::render(
    "docs/detailed_overview.Rmd",
    output_format = rmarkdown::github_document(number_sections = TRUE,
                                               html_preview = FALSE),
    output_file = "detailed_overview.md"
  )
}

# Move the figures into docs/figures/full, point the links at their raw GitHub
# copies (so the Markdown also displays outside GitHub), and use an ATX title.
tidy_md <- function(md = "docs/detailed_overview.md") {
  src <- "docs/detailed_overview_files/figure-gfm"
  dst <- "docs/figures/full"
  if (dir.exists(src)) {
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(dst, full.names = TRUE))
    file.copy(list.files(src, full.names = TRUE), dst)
    unlink("docs/detailed_overview_files", recursive = TRUE)
  }

  txt <- readLines(md, encoding = "UTF-8")
  txt <- gsub("](detailed_overview_files/figure-gfm/",
              "](https://raw.githubusercontent.com/dsmi313/cgrcBayes/main/docs/figures/full/",
              txt, fixed = TRUE)
  sub <- rmarkdown::yaml_front_matter("docs/detailed_overview.Rmd")$subtitle
  if (grepl("^=+$", txt[2])) {
    txt <- c(paste("#", txt[1]), "", if (length(sub)) c(paste("###", sub), ""),
             paste0("*", txt[3], "*"), txt[-(1:3)])
  }
  writeLines(txt, md, useBytes = TRUE)
}

if (sys.nframe() == 0L) {
  render_md()
  tidy_md()
}
