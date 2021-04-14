# internal
#' @title add_abundance_data
#'
#' @param transcriptDB 
#' @param abundance_filename 
#'
#' @return
#' @export
#' 
#' @importFrom readr read_lines read_tsv
#' @importFrom stringr str_detect
#' @importFrom dplyr mutate select left_join if_else
#' @importFrom magrittr set_colnames
#'
#' @examples
add_abundance_data <- function(transcriptDB, abundance_filename) {
  message("Loading abundances...")
  # check if file contains header
  head <- read_lines(abundance_filename, n_max = 20)
  has_header <- any(str_detect(head, "pbid[[:space:]]count_fl"))

  if (has_header){
    abundances <-
      read_tsv(
        file = abundance_filename,
        skip =
          which(
            str_detect(
              string = head,
              pattern =  "pbid[[:space:]]count_fl")
            ) - 1
        )
  } else {
    abundances <- read_tsv(abundance_filename)
  }
  
  abundances <- 
    abundances %>%
    select(c(1,2)) %>%
    set_colnames(c("PBID", "FL_reads"))
  
  # abundances$RawAbundance <- abundances$FL_reads
  # Isoseq is set to minimum 2 FL reads, so add missing 1-read counts
  # should do this better in the future
  transcriptDB <-
    left_join(
      x = transcriptDB,
      y = abundances,
      by = c("PBID")
      ) %>%
    mutate(
      FL_reads = if_else(condition = is.na(FL_reads), true = -1, false = FL_reads)
    )

  transcriptDB
}

# internal
get_gff_data <- function(gff_filepath) {
  # GFF2-compatible only
  message("Loading GFF data...")
  gff_data <- utils::read.table(gff_filepath, header = F, stringsAsFactors = F)
  # keep only the informative columns
  gff_data <- gff_data[,c(1,3:5,7,13)]
  names(gff_data) <- c("Chromosome", "TranscriptExon", "Start", "End",
                       "Strand", "PBID")
  return(gff_data)
}



# internal
make_ORF_db <- function(transcriptDB, aa_cutoff = 250) {
  if (!("ORF" %in% colnames(transcriptDB))) return(NULL)
  # TODO: optimize speed here
  orf_db <- transcriptDB[!is.na(transcriptDB$ORF),
                         c("PBID", "Chromosome", "Gene", "FL_reads", "ORF", "ORFLength")]
  counts <- as.data.frame(table(orf_db$ORF))
  duplicate_proteins <- as.character(counts$Var1[counts$Freq > 1])

  for (i in 1:length(duplicate_proteins)) {
    duplicates <- orf_db[orf_db$ORF == duplicate_proteins[i], ]

    PBIDs <- list(duplicates$PBID)
    merged <- data.frame(ORF        = duplicates$ORF[1],
                         ORFLength  = duplicates$ORFLength[1],
                         Gene       = sort(duplicates$Gene)[1],
                         FL_reads   = sum(as.numeric(duplicates$FL_reads)),
                         Chromosome = duplicates$Chromosome[1])
    merged$PBID <- PBIDs
    orf_db <- orf_db[which(orf_db$ORF != duplicate_proteins[i]), ]
    orf_db <- rbind(orf_db, merged)
  }
  if (max(orf_db$ORFLength) < aa_cutoff) return(NULL)
  # note minimum length cutoff
  orf_db <- orf_db[which(orf_db$ORFLength >= aa_cutoff),]
  orf_db <- calc_abundances(orf_db, genes_known = T)
  return(orf_db)
}

# internal
calc_abundances <- function(df, genes_known = F) {
  # here df is either a TranscriptDB or an OrfDB
  # TODO: optimize speed here
  new_db <- data.frame(matrix(ncol = ncol(df) + 2, nrow = 0))
  colnames(new_db) <- c(names(df), "PercAbundance", "CumPercAbundance")
  if (genes_known) {
    for (gene in unique(df$Gene)) {
      one_gene <- df[df$Gene == gene, ]
      one_gene$PercAbundance <- one_gene$FL_reads / sum(one_gene$FL_reads)
      one_gene <- one_gene[rev(order(one_gene$FL_reads)), ]
      one_gene$CumPercAbundance <- cumsum(one_gene$PercAbundance)
      new_db <- rbind(new_db, one_gene)
    }
  } else {
    # if genes are not known, then PBID prefixes are used as unique gene IDs
    # TODO: apply similar logic in other functions to give RawDBs functionality
    for (pref in unique(df$Prefix)) {
      one_gene <- df[df$Prefix == pref, ]
      one_gene$PercAbundance <- one_gene$FL_reads / sum(one_gene$FL_reads)
      one_gene <- one_gene[rev(order(one_gene$FL_reads)), ]
      one_gene$CumPercAbundance <- cumsum(one_gene$PercAbundance)
      new_db <- rbind(new_db, one_gene)
    }
  }
  return(new_db)
}

# internal
add_gene_names_to_transcripts <- function(transcriptDB, geneDB) {
  # only adds a gene name to a transcript if it doesn't already have one
  # so multiple sources can be used for inputting gene names iteratively
  if (is.null(transcriptDB$Gene)) {
    transcriptDB$Gene <- rep(NA, nrow(transcriptDB))
  }
  for (i in seq_len(length(transcriptDB$Prefix))) {
    name <- geneDB$Name[which(geneDB$ID == transcriptDB$Prefix[i])]
    if (is.na(transcriptDB$Gene[i])) {
      transcriptDB$Gene[i] <- ifelse(length(name) > 0, name, NA)
    }
  }
  return(transcriptDB)
}

# internal
add_gff_info <- function(transcriptDB, gffDB) {
  counts_table <- data.frame(table(gffDB$PBID))
  # subtracting 1 to not count "transcript" lines
  counts <- data.frame(PBID = as.character(counts_table$Var1),
                       ExonCount = counts_table$Freq - 1, stringsAsFactors = F)
  transcriptDB <- transcriptDB %>% dplyr::left_join(counts, by = "PBID")
  gff_transcript_lines <- gffDB[gffDB$TranscriptExon == "transcript", ]
  chromosomes <- data.frame(PBID = gff_transcript_lines$PBID,
                            Chromosome = gff_transcript_lines$Chromosome,
                            stringsAsFactors = F)
  return(transcriptDB %>% dplyr::left_join(chromosomes, by = "PBID"))
}

# internal
add_transcript_lengths <- function(transcriptDB, geneDB) {
  geneDB$AvgTranscriptLength <- sapply(geneDB$Name, function(x) {
    mean(transcriptDB$TranscriptLength[transcriptDB$Gene == x], na.rm = T)
  })
  geneDB$MaxTranscriptLength <- sapply(geneDB$Name, function(x) {
    max(transcriptDB$TranscriptLength[transcriptDB$Gene == x], na.rm = T)
  })
  # this is a weighted average, transcript != read
  geneDB$AvgReadLength <- sapply(geneDB$Name, function(x) {
    stats::weighted.mean(transcriptDB$TranscriptLength[transcriptDB$Gene == x],
                         transcriptDB$FL_reads[transcriptDB$Gene == x], na.rm = T)
  })
  return(geneDB)
}

# internal
add_isoform_info <- function(transcriptDB, geneDB) {
  geneDB$Exons <- sapply(geneDB$Name, function(x) {
    max(transcriptDB$ExonCount[transcriptDB$Gene == x])
  })
  geneDB$N50 <- sapply(geneDB$Name, function(x) {
    sum(transcriptDB$CumPercAbundance[transcriptDB$Gene == x] < 0.5) + 1
  })
  geneDB$N75 <- sapply(geneDB$Name, function(x) {
    sum(transcriptDB$CumPercAbundance[transcriptDB$Gene == x] < 0.75) + 1
  })
  if (length(unique(transcriptDB$Prefix)) == 1) {
    counts <- data.frame(ID = transcriptDB$Prefix[1],
                         Isoforms = nrow(transcriptDB))
  } else {
    counts <- data.frame(table(transcriptDB$Prefix))
  }
  names(counts) <- c("ID", "Isoforms")
  counts$ID <- as.character(counts$ID)
  return(geneDB %>% dplyr::left_join(counts, by = "ID"))
}

# internal
add_ORF_info <- function(transcriptDB, geneDB, orfDB = NULL) {
  if ("ORFLength" %in% colnames(transcriptDB)) {
    geneDB$AvgORFLength <- sapply(geneDB$Name, function(x) {
      mean(transcriptDB$ORFLength[transcriptDB$Gene == x], na.rm = T)
    })
    geneDB$MaxORFLength <- sapply(geneDB$Name, function(x) {
      max(transcriptDB$ORFLength[transcriptDB$Gene == x], na.rm = T)
    })
  }
  if (!is.null(orfDB)) {
    if (length(unique(orfDB$Gene)) == 1) {
      counts <- data.frame(Name = c(unique(orfDB$Gene)),
                           UniqueORFs = c(length(orfDB$Gene)))
    } else {
      counts <- data.frame(table(orfDB$Gene))
      names(counts) <- c("Name", "UniqueORFs")
    }
    counts$Name <- as.character(counts$Name)
    geneDB <- geneDB %>% dplyr::left_join(counts, by = "Name")
    geneDB$ORFN50 <- sapply(geneDB$Name, function(x) {
      sum(orfDB$CumPercAbundance[orfDB$Gene == x] < 0.5) + 1
    })
    geneDB$ORFN75 <- sapply(geneDB$Name, function(x) {
      sum(orfDB$CumPercAbundance[orfDB$Gene == x] < 0.75) + 1
    })
    geneDB$ORFsPerIsoform <- geneDB$UniqueORFs / geneDB$Isoforms
  }
  return(geneDB)
}

# internal
add_diversity_indeces <- function(geneDB, indeces = c("Gini", "Shannon"),
                                  transcriptDB = NULL, orfDB = NULL) {
  if (is.null(transcriptDB) & is.null(orfDB)) {
    stop("Error: no database provided.")
    return(NULL)
  }
  if (is.null(orfDB)) {
    DB <- transcriptDB
    Shannon_str <- "Shannon"
    Gini_str <- "Gini"
  } else {
    DB <- orfDB
    Shannon_str <- "ORFShannon"
    Gini_str <- "ORFGini"
  }

  diversities <- data.frame(Name = geneDB$Name, stringsAsFactors = FALSE)

  if ("Gini" %in% indeces & requireNamespace("ineq", quietly = TRUE)) {
    diversities[Gini_str] <- sapply(diversities$Name, function(gene) {
      ineq::ineq(DB$FL_reads[DB$Gene == gene], type = "Gini")
    })
  }
  if ("Shannon" %in% indeces & requireNamespace("vegan", quietly = TRUE)) {
    diversities[Shannon_str] <- c(sapply(diversities$Name, function(gene) {
      vegan::diversity(DB$FL_reads[DB$Gene == gene])
    }))
  }
  # only merge if packages loaded ok and at least one thing was calculated
  if (dim(diversities)[2] > 1) {
    geneDB <- geneDB %>% dplyr::left_join(diversities, by = "Name")
  }
  return(geneDB)
}


#' @title is_tsv_format
#'
#' @param filename 
#' @param PBIDs 
#'
#' @return
#' @export
#' 
#' @importFrom readr read_lines
#' @importFrom stringr str_detect str_split
#' @importFrom purrr map
#'
#' @examples
is_tsv_format <- function(filename, PBIDs = T) {
  if (PBIDs) {
    read_lines(
      file = filename, 
      n_max = 1
      ) %>%
    str_detect(
      pattern = "PB\\.[0-9]+\\.[0-9]+[[:space:]][acgtnACGTN]+$"
      )
  } else {
    read_lines(
      file = filename,
    ) %>%
    str_split(pattern = "\t") %>%
    map(length) %>%
    unlist() %>%
    any(. == 1)
  }
}

# internal
get_orf_sequences <- function(orf_filename) {
  if (!is_tsv_format(orf_filename)) {
    orfDB <- read_fasta_file(tmpfile, get_prefix = T)
  } else {
    orfDB <- read_tsv_file(orf_filename, get_prefix = T)
  }
  colnames(orfDB) <- c("PBID", "Transcript", "TranscriptLength", "Prefix")
  orfDB
}

# internal
#' @title add_ORF_sequences
#'
#' @param transcriptDB 
#' @param ORF_filename 
#'
#' @return
#' @export
#' 
#' @importFrom dplyr left_join
#'
#' @examples
add_ORF_sequences <- function(transcriptDB, ORF_filename) {
  # the ORF file comes from SQANTI's ORF predictions
  message("Loading ORFs...")
  if (!is_tsv_format(ORF_filename)) {
    ORFs <- read_fasta_file(ORF_filename, get_prefix = F)
  } else {
    ORFs <- read_tsv_file(ORF_filename, get_prefix = F)
  }
  colnames(ORFs) <- c("PBID", "ORF", "ORFLength")


  left_join(
    x = transcriptDB,
    y = ORFs,
    by = "PBID"
    )
}

# internal

#' @title get_transcript_sequences
#'
#' @param transcript_filename 
#'
#' @return
#' @export
#' 
#' @importFrom stringr str_detect str_extract
#' @importFrom magrittr set_colnames
#'
#' @examples
get_transcript_sequences <- function(transcript_filename) {
  message("Loading sequences...")
  if (!str_detect(string = transcript_filename, pattern = "tsv$")) {
    transcriptDB <- read_fasta_file(
      filename = transcript_filename,
      get_prefix = TRUE
    ) 
  } else {
    transcriptDB <-
      read_tsv_file(
        tsv_filename = transcript_filename,
        get_prefix = T
      )
  }
  
  set_colnames(
    x = transcriptDB,
    value = c("PBID", "Transcript", "TranscriptLength", "Prefix"))
}

#' @title read_tsv_file
#'
#' @param tsv_filename 
#' @param get_prefix 
#'
#' @return
#' @export
#' 
#' @importFrom readr read_tsv
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#'
#' @examples
read_tsv_file <- function(tsv_filename, get_prefix = T) {
  seqDB <-
    read_tsv(
      file = tsv_filename,
      skip = 1,
      col_names = c("PBID", "sequence")
      ) %>%
    mutate(
      sequence = toupper(sequence),
      seq_len = nchar(sequence)
    )
  if (get_prefix) {
    seqDB <-
      mutate(
        .data = seqDB,
        Prefix = str_extract(x, "PB\\.[0-9]*")
      )
  }
  seqDB
}

#' @title read_fasta_file
#'
#' @description convert fasta/fastq to tsv format for R-compatible reading in
#' example tsv file line:
#' PB.1.1   ACGTCATCATCAT
#'
#' @param filename
#' @param ORFs
#'
#' @importFrom seqinr read.fasta
#' @importFrom purrr map 
#' @importFrom tidyr unnest
#' @importFrom tibble enframe
#' @importFrom readr write_tsv
#' @importFrom magrittr %>%
#'
read_fasta_file <- function(filename, get_prefix = F) {
  fasta_file <-
    read.fasta(
      filename,
      as.string = TRUE
      ) %>%
    map(
      .f = as.character
    ) %>%
    enframe() %>%
    unnest(col = c(value)) %>%
    mutate(
      value = toupper(value),
      seq_len = nchar(value)
      )
  
  if (get_prefix){
    fasta_file <-
      mutate(
        .data = fasta_file,
        Prefix = str_extract(string = name, pattern = "^PB\\.[[:digit:]]")
      )
  }
  
  fasta_file
}
