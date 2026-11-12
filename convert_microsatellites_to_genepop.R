library(readxl)
library(dplyr)
library(stringr)
library(purrr)

# 1. Read microsatellite data
dat <- read_excel("microsatellites.xlsx", sheet = "microsatellites")

# 2. Define individual and population columns
id_col  <- "Individual"
pop_col <- "Locality"

# 3. Identify allele columns
# Each locus is represented by two consecutive columns: allele 1 and allele 2
allele_cols <- names(dat)[-(1:2)]

# 4. Extract locus names from the first column of each allele pair
original_locus_names <- allele_cols[seq(1, length(allele_cols), by = 2)]

# 5. Use generic locus names, as in the reference GENEPOP file
locus_names <- paste0("LOC", seq_along(original_locus_names))

# 6. Function to recode alleles within each locus as two-digit codes
recode_locus <- function(a1, a2) {
  
  alleles <- c(a1, a2)
  alleles <- alleles[!is.na(alleles) & alleles != -9]
  
  # Unique allele sizes, sorted increasingly
  allele_levels <- sort(unique(alleles))
  
  # Recode each observed allele as 01, 02, 03, ...
  allele_codes <- setNames(
    str_pad(seq_along(allele_levels), width = 2, pad = "0"),
    allele_levels
  )
  
  recode_one <- function(x) {
    if (is.na(x) || x == -9) {
      return("00")
    } else {
      return(allele_codes[as.character(x)])
    }
  }
  
  g1 <- map_chr(a1, recode_one)
  g2 <- map_chr(a2, recode_one)
  
  # Missing genotypes are coded as 0000
  genotype <- ifelse(g1 == "00" | g2 == "00",
                     "0000",
                     paste0(g1, g2))
  
  genotype
}

# 7. Convert all loci to four-digit GENEPOP genotypes
geno_genepop <- map_dfc(seq(1, length(allele_cols), by = 2), function(i) {
  
  a1 <- dat[[allele_cols[i]]]
  a2 <- dat[[allele_cols[i + 1]]]
  
  locus_index <- (i + 1) / 2
  
  tibble(
    !!locus_names[locus_index] := recode_locus(a1, a2)
  )
})

# 8. Add numeric population codes
pop_key <- dat %>%
  distinct(.data[[pop_col]]) %>%
  mutate(pop_code = row_number())

dat_codes <- dat %>%
  left_join(pop_key, by = pop_col)

# 9. Build final GENEPOP data frame
genepop_df <- bind_cols(
  dat_codes[, c(id_col, pop_col, "pop_code")],
  geno_genepop
)

# 10. Write GENEPOP file
out_file <- "microsatellites_genepop.gen"

con <- file(out_file, open = "w")

# Title line
writeLines(
  paste0("gtypes created on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  con
)

# Locus names
writeLines(paste(locus_names, collapse = ", "), con)

# Write populations
for (pop in unique(genepop_df$pop_code)) {
  
  writeLines("POP", con)
  
  sub <- genepop_df %>%
    filter(pop_code == pop)
  
  for (i in seq_len(nrow(sub))) {
    
    geno_line <- paste(as.character(sub[i, locus_names]), collapse = " ")
    
    writeLines(
      paste0(pop, " , ", geno_line),
      con
    )
  }
}

close(con)

# 11. Save population code key
write.csv(
  pop_key,
  file = "population_codes.csv",
  row.names = FALSE
)

# 12. Check first lines of the GENEPOP file
cat(readLines(out_file, n = 30), sep = "\n")