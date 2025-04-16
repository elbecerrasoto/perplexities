#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Globals ----

HITS <- "hits_flanks.tsv"
TAX <- "results/absence_presence.tsv"
REGS <- "regs.tsv"
ISCAN <- "results/iscan.tsv"

WEIGHTS <- c(s = 1, i = 2 / 3, d = 2 / 3, t = 1 / 3)


# Helpers ----

arch_length <- function(arch) {
  if (is.na(arch)) {
    0L
  } else {
    length(str_split_1(arch, pattern = "\\|"))
  }
}

reduce_arch <- function(arch) {
  if (is.na(arch)) {
    return(NA)
  }

  arch_vec <- str_split_1(arch, pattern = "\\|")

  out <- vector(mode = "character", length = length(arch_vec))

  last_dom <- arch_vec[1]
  out[1] <- last_dom
  out_idx <- 1

  for (i in seq_along(arch_vec)) {
    current <- arch_vec[i]
    if (current == last_dom) {
      next
    } else {
      out_idx <- out_idx + 1
      last_dom <- current
      out[out_idx] <- current
    }
  }
  str_flatten(out[out != ""], collapse = "|")
}


one_lettercode <- function(doms) {
  OFFSET <- 33 # avoids non-printable
  INT_LEN <- 6
  LEAD_CHAR <- "IPR"

  pfam_chars <- str_extract(doms, "\\d+")
  stopifnot("Bad ID." = all(str_length(pfam_chars) == INT_LEN))

  pfam_ints <- as.integer(pfam_chars)
  stopifnot("Some extracted IDs are non-numeric." = all(!is.na(pfam_ints)))

  # ?intToUtf8
  # The code points in the surrogate-pair range
  # 0xD800 to 0xDFFF are prohibited in UTF-8 and so are regarded
  # as invalid by utf8ToInt and by default by intToUtf8.

  avoid_invalid <- (0xDFFF - 0xD800) + 1

  target_ints <- pfam_ints + OFFSET
  mask <- target_ints >= 0xD800

  target_ints[mask] <- target_ints[mask] + avoid_invalid

  stopifnot("Unicode points out of range." = all((target_ints) <= 0x10FFFF))
  pfam_codes <- strsplit(intToUtf8(target_ints), "")[[1]]

  stopifnot("Conversion to utf-8 failed." = length(pfam_codes) == length(doms))
  names(pfam_codes) <- doms
  pfam_codes
}


# Main ----

hits <- read_tsv(HITS)
tax <- read_tsv(TAX)
regs <- read_tsv(REGS)
iscan <- read_tsv(ISCAN) |>
  select(pid, length) |>
  group_by(pid) |>
  reframe(length = first(length))

# Unclassified group
otherG <- 0

# wrangle hits
hits <- hits |>
  filter(lengtho_ext == 7, !is.na(archIPR_ext)) |>
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
stopifnot("Non-mutual exclusive groups" = all(!ngroups > 1))

nzeros <- length(which(ngroups == 0))
members_with_groups[ngroups == 0] <- as.list(rep(otherG, nzeros))

stopifnot("Some hits don't have curated group." = all(map_int(members_with_groups, length) == 1))

hits$curated <- unlist(members_with_groups)

# Add reduced Groups

hits <- hits |>
  mutate(
    ARCH = map_chr(archIPR_ext, reduce_arch),
    LARCH = map_int(archIPR_ext, arch_length)
  )
