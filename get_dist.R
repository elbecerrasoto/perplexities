#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Globals ----

HITS <- "hits_flanks.tsv"
TAX <- "results/absence_presence.tsv"
REGS <- "regs.tsv"
ISCAN <- "results/iscan.tsv"

names(tax)

# Helpers ----

hits <- read_tsv(HITS)
tax <- read_tsv(TAX)
regs <- read_tsv(REGS)
iscan <- read_tsv(ISCAN) |>
  select(pid, length) |>
  group_by(pid) |>
  reframe(length = first(length))

# Unassigned group
otherG <- length(regs$reg) + 1

# wrangle hits
hits <- hits |>
  filter(lengtho_ext == 7) |>
  distinct(archIPR_ext, .keep_all = TRUE)

hits <- left_join(hits, iscan, join_by(hit == pid))
hits <- left_join(hits, tax, join_by(genome))

# Add curated

# annotate curated groups
semicolon_archPF <- hits$archPF |>
  str_replace_all("\\|", ";")

groups_with_members <- map(
  regs$reg,
  \(r) which(str_detect(
    semicolon_archPF,
    r
  ))
)


str_detect_regs <- function(x, regs) {
  stopifnot("Only vectorized over regexes." = length(x) == 1)
  map_lgl(regs, \(reg) str_detect(x, reg))
}

f <- partial(str_detect_regs, regs = regs$reg)

members_with_groups <- map(semicolon_archPF, \(x) which(f(x)))
ngroups <- map_int(members_with_groups, length)
stopifnot("Non-mutual exclusive groups" = all(!n_groups > 1))

nzeros <- length(which(ngroups == 0))
members_with_groups[ngroups == 0] <- as.list(rep(otherG, nzeros))

stopifnot("Some hits don't have curated group." = all(map_int(members_with_groups, length) == 1))

hits$curated <- unlist(members_with_groups)

hits$curated
