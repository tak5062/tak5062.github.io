# Codon table (standard)
codon_table <- c(
  "ATA"="I", "ATC"="I", "ATT"="I", "ATG"="M",
  "ACA"="T", "ACC"="T", "ACG"="T", "ACT"="T",
  "AAC"="N", "AAT"="N", "AAA"="K", "AAG"="K",
  "AGC"="S", "AGT"="S", "AGA"="R", "AGG"="R",
  "CTA"="L", "CTC"="L", "CTG"="L", "CTT"="L",
  "CCA"="P", "CCC"="P", "CCG"="P", "CCT"="P",
  "CAC"="H", "CAT"="H", "CAA"="Q", "CAG"="Q",
  "CGA"="R", "CGC"="R", "CGG"="R", "CGT"="R",
  "GTA"="V", "GTC"="V", "GTG"="V", "GTT"="V",
  "GCA"="A", "GCC"="A", "GCG"="A", "GCT"="A",
  "GAC"="D", "GAT"="D", "GAA"="E", "GAG"="E",
  "GGA"="G", "GGC"="G", "GGG"="G", "GGT"="G",
  "TCA"="S", "TCC"="S", "TCG"="S", "TCT"="S",
  "TTC"="F", "TTT"="F", "TTA"="L", "TTG"="L",
  "TAC"="Y", "TAT"="Y", "TAA"="*", "TAG"="*", "TGA"="*",
  "TGC"="C", "TGT"="C", "TGG"="W"
)

aa_properties <- list(
  polar = c("S", "T", "Y", "N", "Q", "C"),
  nonpolar = c("A", "V", "L", "I", "P", "W", "F", "M", "G"),
  acidic = c("D", "E"),
  basic  = c("R", "H", "K")
)

# Reverse complement function
reverse_complement <- function(seq) {
  comp <- c(A="T", T="A", C="G", G="C")
  paste(rev(comp[strsplit(seq, "")[[1]]]), collapse = "")
}

# Function to find ORFs in a strand
find_orfs <- function(seq, strand = "forward") {
  orfs <- list()
  stop_codons <- c("TAA", "TAG", "TGA")
  
  for (frame in 0:2) {
    codons <- substring(seq, seq(1+frame, nchar(seq)-2, 3), seq(3+frame, nchar(seq), 3))
    i <- 1
    while (i <= length(codons)) {
      if (codons[i] == "ATG") {
        for (j in (i+1):length(codons)) {
          if (codons[j] %in% stop_codons) {
            codon_seq <- codons[i:(j-1)]
            aa_seq <- sapply(codon_seq, function(codon) codon_table[[codon]])
            orfs[[length(orfs)+1]] <- list(
              frame = frame + 1,
              strand = strand,
              codons = codon_seq,
              aa_seq = aa_seq
            )
            i <- j
            break
          }
        }
      }
      i <- i + 1
    }
  }
  return(orfs)
}

# Analyze amino acid sequence
analyze_aa_seq <- function(aa_seq) {
  property_counts <- sapply(aa_properties, function(group) sum(aa_seq %in% group))
  aa_counts <- table(factor(aa_seq, levels = unique(unlist(strsplit(paste(codon_table, collapse=""), "")))))
  annotated <- sapply(aa_seq, function(aa) {
    prop <- names(Filter(function(g) aa %in% g, aa_properties))
    if (length(prop) == 0) prop <- "unknown"
    paste(aa, "(", prop, ")", sep = "")
  })
  list(
    properties = property_counts,
    counts = aa_counts,
    annotated = annotated
  )
}

# Main Function
# Get file of DNA sequence
input_file <- "input_dna.txt"

if (!file.exists(input_file)) {
  stop("Input DNA file not found: ", input_file)
}

# Read file, remove newlines and spaces, and convert to uppercase
dna_seq <- toupper(gsub("[^ATCG]", "", paste(readLines(input_file), collapse = "")))

# Optional sanity check
cat("Loaded DNA sequence of length", nchar(dna_seq), "bp\n")


# Get ORFs from forward and reverse strands
forward_orfs <- find_orfs(dna_seq, "forward")
reverse_orfs <- find_orfs(reverse_complement(dna_seq), "reverse")

all_orfs <- c(forward_orfs, reverse_orfs)

# Output to file
output_file <- "full_orf_amino_acid_analysis.txt"
sink(output_file)

cat("=== Complete DNA ORF and Amino Acid Analysis Report ===\n")
cat("DNA Sequence Length:", nchar(dna_seq), "bp\n")
cat("Number of ORFs Found:", length(all_orfs), "\n\n")

for (i in seq_along(all_orfs)) {
  orf <- all_orfs[[i]]
  analysis <- analyze_aa_seq(orf$aa_seq)
  
  cat(sprintf("== ORF #%d (%s strand, frame %d) ==\n", i, orf$strand, orf$frame))
  cat("Amino Acid Sequence:\n", paste(orf$aa_seq, collapse = " "), "\n\n")
  cat("Annotated Sequence:\n", paste(analysis$annotated, collapse = " - "), "\n\n")
  
  cat("Amino Acid Property Counts:\n")
  for (prop in names(analysis$properties)) {
    cat(sprintf("  %s: %d\n", toupper(prop), analysis$properties[[prop]]))
  }
  
  cat("\nAmino Acid Counts:\n")
  for (aa in names(analysis$counts)) {
    cat(sprintf("  %s: %d\n", aa, analysis$counts[[aa]]))
  }
  cat("\n-------------------------------\n\n")
}

sink()
cat("Analysis complete. Report written to:", output_file, "\n")
