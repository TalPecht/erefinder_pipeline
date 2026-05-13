# EREfinder pipeline
# Public functions are kept independent. Shared internals are prefixed with `.`.

.ERE_MATCH_COLS <- c(
  "Perfect_canonical_ERE",
  "perfect_half_site_1bp_sub",
  "perfect_half_site_2bp_sub",
  "perfect_half_site_3bp_sub"
)

.ERE_NUMERIC_COLS <- c(
  "genome_position",
  "N",
  "mean_Kd_inverse",
  .ERE_MATCH_COLS,
  "anymatch"
)

.STANDARD_ENSEMBL_CHR <- c(as.character(1:22), "X", "Y", "MT")
.STANDARD_UCSC_CHR <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

`%||%` <- function(x, y) if (is.null(x)) y else x

.require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.clean_gene_vector <- function(x, name = "genes") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) stop("No valid ", name, " provided.", call. = FALSE)
  x
}

.normalize_genetype <- function(genetype) {
  genetype <- toupper(genetype)
  if (!genetype %in% c("SYMBOL", "ENSEMBL")) {
    stop("genetype must be 'SYMBOL' or 'ENSEMBL'", call. = FALSE)
  }
  genetype
}

.gene_key_col <- function(genetype) {
  if (.normalize_genetype(genetype) == "SYMBOL") "SYMBOL" else "gene_id"
}

.make_mart <- function(mart = NULL, mirror = "useast", dataset = "hsapiens_gene_ensembl") {
  .require_pkgs("biomaRt")
  if (!is.null(mart)) return(mart)
  
  biomaRt::useEnsembl(
    biomart = "genes",
    dataset = dataset,
    mirror = mirror
  )
}

.default_hg38_genome <- function() {
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop(
      "Please install the default genome with:\n",
      "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')",
      call. = FALSE
    )
  }
  
  BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
}

.empty_presence_summary <- function(n_expected = 0L) {
  data.frame(
    n_windows = 0L,
    n_genes_scanned = n_expected,
    n_genes_observed_in_output = 0L,
    n_genes_with_anymatch = 0L,
    n_genes_with_canonical_ere = 0L,
    prop_genes_with_anymatch = if (!is.na(n_expected) && n_expected > 0) 0 else NA_real_,
    prop_genes_with_canonical_ere = if (!is.na(n_expected) && n_expected > 0) 0 else NA_real_,
    total_anymatch = 0,
    total_perfect_canonical_ERE = 0,
    mean_anymatch = 0,
    mean_perfect_canonical_ERE = 0,
    stringsAsFactors = FALSE
  )
}

.empty_class_counts <- function(n_expected = 0L) {
  data.frame(
    match_type = .ERE_MATCH_COLS,
    n_genes_scanned = n_expected,
    n_genes_with_hit = 0L,
    prop_genes_with_hit = if (!is.na(n_expected) && n_expected > 0) 0 else NA_real_,
    stringsAsFactors = FALSE
  )
}

.safe_run_erefinder <- function(seqs, label, ...) {
  tryCatch(
    run_EREfinder(seqs = seqs, out_prefix = label, ...),
    error = function(e) {
      warning("EREfinder failed for ", label, ": ", e$message, call. = FALSE)
      list(
        results = NULL,
        summary = .empty_presence_summary(0L),
        class_counts = .empty_class_counts(0L),
        run_info = list(error = e$message)
      )
    }
  )
}

# Function: get_one_promoter_sequence 

#' Retrieve one promoter sequence per gene from Ensembl
#'
#' Queries Ensembl BioMart for transcript annotations and retrieves one promoter
#' region per gene using the canonical transcript when available. Promoter
#' sequences are extracted from the supplied genome object.
#'
#' @param genes Character vector of gene symbols or Ensembl gene IDs.
#' @param genetype Either `"SYMBOL"` or `"ENSEMBL"`.
#' @param upstream Number of bases upstream of the transcription start site.
#' @param downstream Number of bases downstream of the transcription start site.
#' @param genome Optional BSgenome object. If NULL, hg38 is used.
#' @param dataset Ensembl dataset name. Defaults to human genes.
#' @param unique_names Logical; whether to generate unique FASTA sequence names.
#' @param mart Optional BioMart object. If NULL, a new Ensembl mart is created.
#' @param mirror Ensembl mirror to use when creating a BioMart connection.
#'
#' @return A list with `seqs`, a DNAStringSet of promoter sequences, and `meta`,
#' a data.frame of promoter metadata.
#'
#' @export
get_one_promoter_sequence <- function(
    genes,
    genetype = "SYMBOL",
    upstream = 1000,
    downstream = 200,
    genome = NULL,
    dataset = "hsapiens_gene_ensembl",
    unique_names = TRUE,
    mart = NULL,
    mirror = "useast"
) {
  .require_pkgs(c("biomaRt", "GenomicRanges", "GenomeInfoDb", "BSgenome", "IRanges"))
  
  if (is.null(genome)) genome <- .default_hg38_genome()
  
  genetype <- .normalize_genetype(genetype)
  genes <- .clean_gene_vector(genes)
  mart <- .make_mart(mart = mart, mirror = mirror, dataset = dataset)
  
  attrs <- c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "chromosome_name",
    "strand",
    "transcription_start_site",
    "transcript_is_canonical"
  )
  
  bm <- biomaRt::getBM(
    attributes = attrs,
    filters = if (genetype == "SYMBOL") "hgnc_symbol" else "ensembl_gene_id",
    values = genes,
    mart = mart
  )
  
  if (!nrow(bm)) stop("No Ensembl records found for the supplied genes.", call. = FALSE)
  
  if (ncol(bm) != 7) { # verify all columns were returned before renaming
    stop("Unexpected number of columns returned from biomaRt.", call. = FALSE)
  }
  
  names(bm) <- c(
    "SYMBOL",
    "gene_id",
    "tx_name",
    "chromosome_name",
    "strand",
    "tss",
    "transcript_is_canonical"
  )
  
  
  required_cols <- c("gene_id", "tx_name", "chromosome_name", "tss") # remove rows missing any of the fields
  
  keep_complete <- Reduce(`&`, lapply(required_cols, function(col) {
    x <- bm[[col]]
    !is.na(x) & nzchar(trimws(as.character(x)))
  }))
  
  bm <- bm[keep_complete, , drop = FALSE] # keep only rows with all attributes present
  
  if (genetype == "SYMBOL") { 
    bm <- bm[!is.na(bm$SYMBOL) & nzchar(bm$SYMBOL), , drop = FALSE] # Keep only rows where SYMBOL is not NA and not an empty string.
  }
  
  bm <- bm[bm$chromosome_name %in% .STANDARD_ENSEMBL_CHR, , drop = FALSE] # Keep only rows whose chromosome is in the allowed standard Ensembl chromosome list
  if (!nrow(bm)) stop("No usable transcript coordinates returned by BioMart.", call. = FALSE)
  
  dedup_col <- .gene_key_col(genetype) # defines which column is used to remove duplication (symbol or geneid)
  
  bm$canon_rank <- ifelse(is.na(bm$transcript_is_canonical), 0L, as.integer(bm$transcript_is_canonical)) # Create a new column called canon_rank. If transcript_is_canonical is missing, set it to 0.  Otherwise convert it to an integer.  So canonical transcripts get ranked higher
  
  bm <- bm[order(bm[[dedup_col]], -bm$canon_rank, bm$gene_id, bm$tx_name), , drop = FALSE] # sort rows - put cannonical (1) transcript first
  
  bm_one <- bm[!duplicated(bm[[dedup_col]]), , drop = FALSE] # Keep only the first row for each gene key
  
  prom_start <- ifelse(bm_one$strand == 1, bm_one$tss - upstream, bm_one$tss - downstream) # Calculate seq start coordinate: For genes on the + strand, seq starts upstream of the TSS, For genes on the - strand, “upstream” is in the opposite genomic direction
  
  prom_end <- ifelse(bm_one$strand == 1, bm_one$tss + downstream - 1, bm_one$tss + upstream - 1) # Calculate seq end coordinate
  
  gr <- GenomicRanges::GRanges( # Create a GRanges object, which stores genomic intervals.
    seqnames = paste0("chr", bm_one$chromosome_name), # Set chromosome names by adding "chr"
    ranges = IRanges::IRanges(start = prom_start, end = prom_end), # Set the genomic start and end coordinates
    strand = ifelse(bm_one$strand == 1, "+", "-") # Set strand as "+" or "-"
  )
  
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC" # set chromosome names are in UCSC style, like chr1, chr2, chrX
  
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)[GenomeInfoDb::seqlevels(gr)] # Attach chromosome length and genome metadata from the BSgenome object
  
  gr <- GenomicRanges::trim(gr) # Trim ranges that go outside chromosome boundaries.
  
  keep <- as.character(GenomicRanges::seqnames(gr)) %in% .STANDARD_UCSC_CHR & # building a logical vector saying which ranges to keep: in STANDARD_UCSC_CHR and ange must have positive width
    GenomicRanges::width(gr) > 0
  
  gr <- gr[keep] # Filter the seq ranges
  
  bm_one <- bm_one[keep, , drop = FALSE] # Filter the metadata rows the same way, so bm_one still matches gr.
  
  if (!length(gr)) stop("No valid seq ranges remain after standard chromosome filtering.", call. = FALSE) # If no promoter ranges remain, stop with an error.
  
  seqs <- BSgenome::getSeq(genome, gr) # Extract the actual DNA sequences from the genome for those promoter ranges.
  
  # build a metadata table describing each  sequence extracted
  meta <- data.frame(
    input_gene = bm_one[[dedup_col]],
    gene_id = bm_one$gene_id,
    SYMBOL = bm_one$SYMBOL,
    tx_name = bm_one$tx_name,
    transcript_is_canonical = bm_one$transcript_is_canonical,
    chr = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    strand = as.character(GenomicRanges::strand(gr)),
    promoter_width = GenomicRanges::width(gr),
    stringsAsFactors = FALSE
  )
  
  names(seqs) <- if (unique_names) {
    paste(meta$input_gene, meta$gene_id, meta$SYMBOL, meta$tx_name, sep = "|")
  } else {
    meta$input_gene
  }
  
  list(seqs = seqs, meta = meta)
}

#' Parse EREfinder output
#'
#' Parses the text output generated by EREfinder into a tidy data.frame.
#'
#' @param output_file Path to an EREfinder output file.
#'
#' @return A data.frame of EREfinder window-level results, or NULL if no valid
#' rows were parsed.
#'
#' @export
parse_EREfinder_output <- function(output_file) {
  lines <- readLines(output_file, warn = FALSE) # Reads the entire text file at output_file into a character vector called lines
  gene <- NA_character_ # Initializes the current gene name as missing
  rows <- vector("list", length(lines)) # Creates an empty list called rows. The list has the same length as the number of lines in the file. Each valid parsed data row will later be stored as one element of this list.
  n <- 0L # Initializes a counter called n.
  
  for (line in lines) { # Starts a loop over every line in the file.
    line <- trimws(line) # Removes leading and trailing whitespace from the current line
    if (!nzchar(line)) next # Skips empty lines
    
    if (grepl("^>", line)) { # Checks whether the line starts with >
      gene <- sub("^>", "", line) # Removes the leading > and stores the rest as the current gene name.
      next
    }
    
    if (grepl("genome_position", line, fixed = TRUE)) next # Skips header rows that contain the text "genome_position"
    
    parts <- strsplit(line, ",", fixed = TRUE)[[1]] # Splits the current line by commas.
    if (length(parts) < length(.ERE_NUMERIC_COLS)) next # Skips the line if it does not contain enough comma-separated fields.
    
    n <- n + 1L
    rows[[n]] <- data.frame(
      gene = gene,
      as.list(stats::setNames(as.numeric(parts[seq_along(.ERE_NUMERIC_COLS)]), .ERE_NUMERIC_COLS)), # Turns the named numeric vector into a list.
      stringsAsFactors = FALSE, # Prevents character columns from being converted into factors.
      check.names = FALSE # Prevents R from modifying column names.
    )
  }
  
  rows <- rows[seq_len(n)] # Keeps only the filled elements of rows.
  if (!length(rows)) { # Checks whether no valid rows were parsed.
    warning("No valid EREfinder data parsed.", call. = FALSE)
    return(NULL)
  }
  
  do.call(rbind, rows) # Combines all one-row data frames into a single data frame.
}

#' Run EREfinder on selected regions
#'
#' Writes sequences to a temporary FASTA file, runs the EREfinder executable, and
#' returns parsed results plus useful R summaries. By default, files are written
#' to `tempdir()` and deleted after the run.
#'
#' @param seqs A DNAStringSet of sequences.
#' @param gene_ids Optional clean gene IDs corresponding to `seqs`, usually
#' `promoters_obj$meta$input_gene`. If supplied, EREfinder's technical FASTA
#' names are mapped back to these gene IDs before summaries/statistics.
#' @param out_prefix Prefix for FASTA and output files.
#' @param erefinder_path Path to the EREfinder executable.
#' @param out_dir Directory for temporary output files. Defaults to `tempdir()`.
#' @param alphabet Alphabet mode passed to EREfinder.
#' @param quiet Logical; suppress command-line output.
#' @param cleanup Logical; remove temporary files after execution.
#' @param window Window size used by EREfinder.
#' @param step Step size used by EREfinder.
#' @param return_metadata Logical; if TRUE, return a list containing results,
#' summaries, class counts, and run information. If FALSE, return only the parsed
#' data.frame.
#'
#' @return If `return_metadata = TRUE`, a list with `results`, `summary`,
#' `class_counts`, and `run_info`. If `return_metadata = FALSE`, a data.frame or
#' NULL.
#'
#' @export
run_EREfinder <- function(
    seqs,
    gene_ids = NULL,
    out_prefix = paste0("erefinder_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    erefinder_path,
    out_dir = tempdir(),
    alphabet = "a",
    quiet = TRUE,
    cleanup = TRUE,
    window = 200,
    step = 50,
    return_metadata = TRUE
) {
  .require_pkgs(c("Biostrings", "utils"))
  
  if (!length(seqs)) {
    warning("No sequences supplied to run_EREfinder().", call. = FALSE)
    return(NULL)
  }
  
  if (is.null(names(seqs)) || any(!nzchar(names(seqs)))) {
    names(seqs) <- paste0("seq_", seq_along(seqs))
  }
  
  if (!is.null(gene_ids)) {
    gene_ids <- as.character(gene_ids)
    if (length(gene_ids) != length(seqs)) {
      stop("gene_ids must be NULL or have the same length as seqs.", call. = FALSE)
    }
  }
  
  expected_genes <- if (is.null(gene_ids)) names(seqs) else unique(gene_ids)
  
  if (missing(erefinder_path) || !nzchar(erefinder_path)) {
    stop("erefinder_path is required.", call. = FALSE)
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fasta_file <- file.path(out_dir, paste0(out_prefix, ".fa"))
  output_file <- file.path(out_dir, paste0(out_prefix, "_output.csv"))
  
  if (cleanup) {
    on.exit(unlink(c(fasta_file, output_file)), add = TRUE)
  }
  
  Biostrings::writeXStringSet(seqs, fasta_file)
  
  args <- c(
    "-i", fasta_file,
    "-o", output_file,
    "-a", alphabet,
    "-w", as.character(window),
    "-d", as.character(step)
  )
  
  cmd_output <- system2(
    erefinder_path,
    args,
    stdout = if (quiet) TRUE else "",
    stderr = if (quiet) TRUE else ""
  )
  
  exit_status <- attr(cmd_output, "status")
  if (is.null(exit_status)) exit_status <- 0L
  
  run_info <- list(
    status = exit_status,
    command = paste(shQuote(erefinder_path), paste(shQuote(args), collapse = " ")),
    stdout_stderr = cmd_output,
    n_sequences = length(seqs),
    sequence_names = names(seqs),
    gene_ids = gene_ids,
    window = window,
    step = step,
    alphabet = alphabet,
    fasta_file = if (cleanup) NA_character_ else fasta_file,
    output_file = if (cleanup) NA_character_ else output_file,
    cleanup = cleanup
  )
  
  if (!file.exists(output_file)) {
    warning("EREfinder output file not found. Check input, executable path, and arguments.", call. = FALSE)
    
    if (return_metadata) {
      return(list(
        results = NULL,
        summary = summarise_ere_hits_presence(NULL, expected_genes = expected_genes),
        class_counts = count_ere_gene_presence_by_class(NULL, expected_genes = expected_genes),
        run_info = run_info
      ))
    }
    
    return(NULL)
  }
  
  results <- parse_EREfinder_output(output_file)
  
  if (!is.null(results) && !is.null(gene_ids)) {
    name_map <- data.frame(
      sequence_name = names(seqs),
      gene_id_clean = gene_ids,
      stringsAsFactors = FALSE
    )
    
    match_idx <- match(results$gene, name_map$sequence_name)
    results$sequence_name <- results$gene
    results$gene <- name_map$gene_id_clean[match_idx]
    
    unmapped <- is.na(results$gene)
    if (any(unmapped)) {
      warning(
        sum(unmapped),
        " EREfinder rows could not be mapped back to gene_ids; keeping original FASTA names for those rows.",
        call. = FALSE
      )
      results$gene[unmapped] <- results$sequence_name[unmapped]
    }
  }
  
  if (!return_metadata) return(results)
  
  list(
    results = results,
    summary = summarise_ere_hits_presence(results, expected_genes = expected_genes),
    class_counts = count_ere_gene_presence_by_class(results, expected_genes = expected_genes),
    run_info = run_info
  )
}

#' Subset promoter sequences by gene set
#'
#' Filters a promoter object returned by `get_one_promoter_sequence()` to retain
#' only selected genes.
#'
#' @param promoters_obj A promoter object containing `seqs` and `meta`.
#' @param genes Character vector of genes to retain.
#' @param genetype Either `"SYMBOL"` or `"ENSEMBL"`.
#'
#' @return A filtered promoter object containing `seqs` and `meta`, or NULL if
#' no genes matched.
#'
#' @export
subset_promoters_by_genes <- function(promoters_obj, genes, genetype = "SYMBOL") {
  genes <- .clean_gene_vector(genes)
  key_col <- .gene_key_col(genetype)
  
  # Check that promoters_obj has the expected structure.
  if (!is.list(promoters_obj) || is.null(promoters_obj$meta) || is.null(promoters_obj$seqs)) {
    stop("promoters_obj must be a list with `seqs` and `meta`.", call. = FALSE)
  }
  
  # Check whether the chosen gene ID column exists in the metadata table.
  if (!key_col %in% names(promoters_obj$meta)) {
    stop("Column `", key_col, "` not found in promoters_obj$meta.", call. = FALSE)
  }
  
  # Find the row positions in promoters_obj$meta where the gene column matches one of the requested genes.
  keep_idx <- which(promoters_obj$meta[[key_col]] %in% genes)
  if (!length(keep_idx)) return(NULL) # If there are zero matching promoter sequences, return nothing.
  
  # Return a filtered version of promoters_obj.
  list(
    seqs = promoters_obj$seqs[keep_idx],
    meta = promoters_obj$meta[keep_idx, , drop = FALSE] # keeps the result as a data frame, even if only one row matches.
  )
}

#' Summarise EREfinder hits
#'
#' Produces a compact summary of EREfinder results across all windows.
#'
#' @param df A data.frame returned by `run_EREfinder()` or the `results` element
#' from a `run_EREfinder()` result list.
#'
#' @return A one-row data.frame containing counts and averages of ERE hits.
#'
#' @export
summarise_ere_hits <- function(df) {
  
  # Checks whether there is no usable data.
  if (is.null(df) || !nrow(df)) {
    return(data.frame(
      n_windows = 0L,
      n_genes = 0L,
      total_anymatch = 0,
      total_perfect_canonical_ERE = 0,
      mean_anymatch = NA_real_,
      mean_perfect_canonical_ERE = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    n_windows = nrow(df),
    n_genes = length(unique(df$gene)),
    total_anymatch = sum(df$anymatch, na.rm = TRUE),
    total_perfect_canonical_ERE = sum(df$Perfect_canonical_ERE, na.rm = TRUE),
    mean_anymatch = mean(df$anymatch, na.rm = TRUE),
    mean_perfect_canonical_ERE = mean(df$Perfect_canonical_ERE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Count gene-level ERE presence by match class
#'
#' Counts how many genes contain at least one positive EREfinder window within
#' each ERE match class.
#'
#' @param df A data.frame returned by `run_EREfinder()` or the `results` element
#' from a `run_EREfinder()` result list.
#' @param expected_genes Optional character vector of expected genes used as the
#' denominator for proportions.
#'
#' @return A data.frame with one row per ERE match class.
#'
#' @export
count_ere_gene_presence_by_class <- function(df, expected_genes = NULL) {
  n_expected <- if (is.null(expected_genes)) 0L else length(unique(expected_genes))
  if (is.null(df) || !nrow(df)) return(.empty_class_counts(n_expected))
  
  scanned_genes <- if (is.null(expected_genes)) unique(df$gene) else unique(expected_genes)
  n_scanned <- length(scanned_genes)
  
  out <- lapply(.ERE_MATCH_COLS, function(col) {
    genes_with_hit <- unique(df$gene[!is.na(df[[col]]) & df[[col]] > 0])
    data.frame(
      match_type = col,
      n_genes_scanned = n_scanned,
      n_genes_with_hit = length(genes_with_hit),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  out$prop_genes_with_hit <- if (n_scanned > 0) out$n_genes_with_hit / n_scanned else NA_real_
  out
}

#' Summarise gene-level ERE presence
#'
#' Produces a broad summary of ERE presence across all genes, including any
#' ERE-like matches and canonical ERE matches.
#'
#' @param df A data.frame returned by `run_EREfinder()` or the `results` element
#' from a `run_EREfinder()` result list.
#' @param expected_genes Optional character vector of expected genes used as the
#' denominator for proportions.
#'
#' @return A one-row data.frame summarising ERE presence statistics.
#'
#' @export
summarise_ere_hits_presence <- function(df, expected_genes = NULL) {
  n_expected <- if (is.null(expected_genes)) 0L else length(unique(expected_genes))
  if (is.null(df) || !nrow(df)) return(.empty_presence_summary(n_expected))
  
  observed_genes <- unique(df$gene)
  scanned_genes <- if (is.null(expected_genes)) observed_genes else unique(expected_genes)
  n_scanned <- length(scanned_genes)
  
  genes_any <- unique(df$gene[!is.na(df$anymatch) & df$anymatch > 0])
  genes_canonical <- unique(df$gene[!is.na(df$Perfect_canonical_ERE) & df$Perfect_canonical_ERE > 0])
  
  n_any <- length(intersect(genes_any, scanned_genes))
  n_canonical <- length(intersect(genes_canonical, scanned_genes))
  
  data.frame(
    n_windows = nrow(df),
    n_genes_scanned = n_scanned,
    n_genes_observed_in_output = length(observed_genes),
    n_genes_with_anymatch = n_any,
    n_genes_with_canonical_ere = n_canonical,
    prop_genes_with_anymatch = if (n_scanned > 0) n_any / n_scanned else NA_real_,
    prop_genes_with_canonical_ere = if (n_scanned > 0) n_canonical / n_scanned else NA_real_,
    total_anymatch = sum(df$anymatch, na.rm = TRUE),
    total_perfect_canonical_ERE = sum(df$Perfect_canonical_ERE, na.rm = TRUE),
    mean_anymatch = mean(df$anymatch, na.rm = TRUE),
    mean_perfect_canonical_ERE = mean(df$Perfect_canonical_ERE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

#' Calculate an empirical permutation p-value
#'
#' Compares an observed statistic against a permutation distribution using an
#' empirical p-value calculation.
#'
#' @param obs Observed statistic.
#' @param perm Numeric vector of permutation statistics.
#' @param alternative One of `"greater"`, `"less"`, or `"two.sided"`.
#'
#' @return A one-row data.frame containing observed value, permutation summary,
#' empirical p-value, and fold enrichment.
#'
#' @export
calc_empirical_test <- function(obs, perm, alternative = "greater") {
  if (!alternative %in% c("greater", "less", "two.sided")) {
    stop("alternative must be 'greater', 'less', or 'two.sided'", call. = FALSE)
  }
  
  perm <- perm[!is.na(perm)]
  
  if (length(obs) != 1 || is.na(obs) || !length(perm)) {
    return(data.frame(
      observed = obs,
      mean_perm = NA_real_,
      median_perm = NA_real_,
      n_perm = length(perm),
      p_value = NA_real_,
      fold_enrichment_vs_mean_perm = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  mean_perm <- mean(perm)
  
  p_value <- switch(
    alternative,
    greater = (sum(perm >= obs) + 1) / (length(perm) + 1),
    less = (sum(perm <= obs) + 1) / (length(perm) + 1),
    two.sided = (sum(abs(perm - mean_perm) >= abs(obs - mean_perm)) + 1) / (length(perm) + 1)
  )
  
  data.frame(
    observed = obs,
    mean_perm = mean_perm,
    median_perm = stats::median(perm),
    n_perm = length(perm),
    p_value = p_value,
    fold_enrichment_vs_mean_perm = if (!is.na(mean_perm) && mean_perm != 0) obs / mean_perm else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Run the full ERE promoter enrichment pipeline
#'
#' Retrieves promoter sequences, runs EREfinder on genes of interest and a
#' background universe, performs optional permutations, and calculates empirical
#' enrichment statistics.
#'
#' @param genes_of_interest Character vector of genes of interest.
#' @param universe Character vector defining the background gene universe.
#' @param genetype Either `"SYMBOL"` or `"ENSEMBL"`.
#' @param seed Random seed for permutation reproducibility.
#' @param n_perm Number of permutations.
#' @param run_permutations Logical; whether to run permutation testing.
#' @param upstream Number of upstream promoter bases.
#' @param downstream Number of downstream promoter bases.
#' @param genome Optional BSgenome object. If NULL, hg38 is used.
#' @param dataset Ensembl dataset name. Defaults to human genes.
#' @param unique_names Logical; whether FASTA names should be unique.
#' @param mirror Ensembl mirror used for BioMart queries.
#' @param mart Optional BioMart object.
#' @param erefinder_path Path to the EREfinder executable.
#' @param out_dir Directory for EREfinder temporary files. Defaults to `tempdir()`.
#' @param out_prefix Prefix for generated files.
#' @param alphabet Alphabet mode passed to EREfinder.
#' @param quiet Logical; suppress command-line output.
#' @param cleanup Logical; remove temporary files.
#' @param window EREfinder window size.
#' @param step EREfinder step size.
#' @param show_progress Logical; show progress bar during permutations.
#' @param p.val.alternative Alternative hypothesis used for empirical p-values.
#'
#' @return A nested list containing input information, filtered genes, promoter
#' objects, EREfinder results, run metadata, summaries, and p-values.
#'
#' @export
run_ere_pipeline <- function(
    genes_of_interest,
    universe,
    genetype = "SYMBOL",
    seed = 42,
    n_perm = 1000,
    run_permutations = TRUE,
    upstream = 10000,
    downstream = 2000,
    genome = NULL,
    dataset = "hsapiens_gene_ensembl",
    unique_names = TRUE,
    mirror = "useast",
    mart = NULL,
    erefinder_path,
    out_dir = tempdir(),
    out_prefix = "ere_pipeline",
    alphabet = "a",
    quiet = TRUE,
    cleanup = TRUE,
    window = 1000,
    step = 1000,
    show_progress = TRUE,
    p.val.alternative = "greater"
) {
  .require_pkgs(c("biomaRt", "BSgenome", "Biostrings", "GenomicRanges", "GenomeInfoDb", "IRanges"))
  
  genetype <- .normalize_genetype(genetype)
  universe <- .clean_gene_vector(universe, "universe genes")
  genes_of_interest <- .clean_gene_vector(genes_of_interest, "genes_of_interest")
  
  if (is.null(genome)) genome <- .default_hg38_genome()
  
  old_timeout <- getOption("timeout")
  old_biomart_cache <- Sys.getenv("BIOMART_CACHE", unset = NA_character_)
  
  options(timeout = max(300, old_timeout %||% 60))
  Sys.setenv(BIOMART_CACHE = tempdir())
  
  on.exit({
    options(timeout = old_timeout)
    if (is.na(old_biomart_cache)) {
      Sys.unsetenv("BIOMART_CACHE")
    } else {
      Sys.setenv(BIOMART_CACHE = old_biomart_cache)
    }
  }, add = TRUE)
  
  mart <- .make_mart(mart = mart, mirror = mirror, dataset = dataset)
  
  promoters_universe_raw <- get_one_promoter_sequence(
    genes = universe,
    genetype = genetype,
    upstream = upstream,
    downstream = downstream,
    genome = genome,
    dataset = dataset,
    unique_names = unique_names,
    mart = mart,
    mirror = mirror
  )
  
  universe_key <- .gene_key_col(genetype)
  valid_universe <- unique(promoters_universe_raw$meta[[universe_key]])
  valid_universe <- valid_universe[!is.na(valid_universe) & nzchar(valid_universe)]
  
  goi_valid <- intersect(genes_of_interest, valid_universe)
  if (!length(goi_valid)) {
    stop("None of the genes_of_interest were found in the valid promoter universe.", call. = FALSE)
  }
  
  promoters_universe <- subset_promoters_by_genes(promoters_universe_raw, valid_universe, genetype)
  promoters_goi <- subset_promoters_by_genes(promoters_universe_raw, goi_valid, genetype)
  
  ere_args <- list(
    erefinder_path = erefinder_path,
    out_dir = out_dir,
    alphabet = alphabet,
    quiet = quiet,
    cleanup = cleanup,
    window = window,
    step = step,
    return_metadata = TRUE
  )
  
  ere_goi_obj <- do.call(.safe_run_erefinder, c(
    list(
      seqs = promoters_goi$seqs,
      gene_ids = promoters_goi$meta$input_gene,
      label = paste0(out_prefix, "_goi")
    ),
    ere_args
  ))
  
  ere_universe_obj <- do.call(.safe_run_erefinder, c(
    list(
      seqs = promoters_universe$seqs,
      gene_ids = promoters_universe$meta$input_gene,
      label = paste0(out_prefix, "_universe")
    ),
    ere_args
  ))
  
  ere_goi <- ere_goi_obj$results
  ere_universe <- ere_universe_obj$results
  
  summary_goi <- ere_goi_obj$summary
  summary_universe <- ere_universe_obj$summary
  class_counts_goi <- ere_goi_obj$class_counts
  class_counts_universe <- ere_universe_obj$class_counts
  
  perm_results <- NULL
  perm_summaries <- NULL
  perm_class_counts <- NULL
  perm_run_info <- NULL
  
  if (run_permutations) {
    set.seed(seed)
    n_test <- length(goi_valid)
    
    perm_results <- vector("list", n_perm)
    perm_summaries <- vector("list", n_perm)
    perm_class_counts <- vector("list", n_perm)
    perm_run_info <- vector("list", n_perm)
    
    if (show_progress) {
      pb <- utils::txtProgressBar(min = 0, max = n_perm, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    
    for (i in seq_len(n_perm)) {
      perm_genes <- sample(valid_universe, n_test, replace = FALSE)
      
      perm_obj <- tryCatch({
        perm_promoters <- subset_promoters_by_genes(promoters_universe_raw, perm_genes, genetype)
        
        if (is.null(perm_promoters) || nrow(perm_promoters$meta) != length(perm_genes)) {
          stop(
            "Permutation ", i, " mismatch: ",
            length(perm_genes), " genes sampled, but ",
            if (is.null(perm_promoters)) 0 else nrow(perm_promoters$meta),
            " promoter entries were returned."
          )
        }
        
        perm_ere_obj <- do.call(.safe_run_erefinder, c(
          list(
            seqs = perm_promoters$seqs,
            gene_ids = perm_promoters$meta$input_gene,
            label = paste0(out_prefix, "_perm_", i)
          ),
          ere_args
        ))
        
        list(
          run = i,
          genes = perm_genes,
          ere = perm_ere_obj$results,
          summary = perm_ere_obj$summary,
          class_counts = perm_ere_obj$class_counts,
          run_info = perm_ere_obj$run_info
        )
      }, error = function(e) {
        message("Permutation ", i, " failed: ", e$message)
        list(
          run = i,
          genes = perm_genes,
          ere = NULL,
          summary = .empty_presence_summary(NA_integer_),
          class_counts = transform(
            .empty_class_counts(NA_integer_),
            n_genes_with_hit = NA_integer_,
            prop_genes_with_hit = NA_real_
          ),
          run_info = list(error = e$message),
          error = e$message
        )
      })
      
      perm_results[[i]] <- perm_obj
      perm_summaries[[i]] <- cbind(run = i, perm_obj$summary)
      perm_class_counts[[i]] <- cbind(run = i, perm_obj$class_counts)
      perm_run_info[[i]] <- perm_obj$run_info
      
      if (show_progress) utils::setTxtProgressBar(pb, i)
    }
    
    perm_summaries <- do.call(rbind, perm_summaries)
    perm_class_counts <- do.call(rbind, perm_class_counts)
  }
  
  empirical_p_anymatch <- NA_real_
  empirical_p_canonical <- NA_real_
  
  if (run_permutations && !is.null(perm_summaries) && nrow(perm_summaries)) {
    empirical_p_anymatch <- calc_empirical_test(
      obs = summary_goi$n_genes_with_anymatch,
      perm = perm_summaries$n_genes_with_anymatch,
      alternative = p.val.alternative
    )
    
    empirical_p_canonical <- calc_empirical_test(
      obs = summary_goi$n_genes_with_canonical_ere,
      perm = perm_summaries$n_genes_with_canonical_ere,
      alternative = p.val.alternative
    )
  }
  
  class_pvalues <- NULL
  
  if (run_permutations && !is.null(perm_class_counts) && nrow(perm_class_counts)) {
    class_pvalues <- do.call(rbind, lapply(seq_len(nrow(class_counts_goi)), function(j) {
      this_class <- class_counts_goi$match_type[j]
      
      test_res <- calc_empirical_test(
        obs = class_counts_goi$n_genes_with_hit[j],
        perm = perm_class_counts$n_genes_with_hit[perm_class_counts$match_type == this_class],
        alternative = p.val.alternative
      )
      
      cbind(class_counts_goi[j, , drop = FALSE], test_res)
    }))
  }
  
  list(
    input = list(
      genes_of_interest = genes_of_interest,
      universe = universe,
      genetype = genetype,
      dataset = dataset,
      n_perm = n_perm,
      run_permutations = run_permutations,
      upstream = upstream,
      downstream = downstream,
      window = window,
      step = step
    ),
    filtered = list(
      valid_universe = valid_universe,
      genes_of_interest_valid = goi_valid,
      genes_of_interest_missing_from_valid_universe = setdiff(genes_of_interest, goi_valid)
    ),
    promoters = list(
      universe = promoters_universe,
      genes_of_interest = promoters_goi
    ),
    erefinder = list(
      genes_of_interest = ere_goi,
      universe = ere_universe,
      permutations = perm_results
    ),
    erefinder_runs = list(
      genes_of_interest = ere_goi_obj$run_info,
      universe = ere_universe_obj$run_info,
      permutations = perm_run_info
    ),
    summaries = list(
      genes_of_interest = summary_goi,
      universe = summary_universe,
      permutations = perm_summaries,
      empirical_p_anymatch = empirical_p_anymatch,
      empirical_p_canonical = empirical_p_canonical
    ),
    class_specific = list(
      genes_of_interest = class_counts_goi,
      universe = class_counts_universe,
      permutations = perm_class_counts,
      pvalues = class_pvalues
    ),
    mart = mart
  )
}

# Visualization functions for ERE pipeline results

#' Plot permutation distribution
#'
#' Visualize the permutation distribution and observed statistic.
#'
#' @param observed Observed statistic.
#' @param perm Numeric vector of permutation statistics.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param bins Number of histogram bins.
#'
#' @return A ggplot object.
#' @export
plot_perm_distribution <- function(
    observed,
    perm,
    title = "Permutation distribution",
    xlab = "Statistic",
    bins = 30
) {
  .require_pkgs("ggplot2")
  
  perm <- perm[!is.na(perm)]
  
  df <- data.frame(value = perm)
  
  ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey70", color = "black") +
    ggplot2::geom_vline(xintercept = observed, color = "red", linewidth = 1.2) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = "Count"
    ) +
    ggplot2::theme_bw()
}


#' Plot class-specific gene hit counts
#'
#' Barplot showing how many genes contain at least one hit in each ERE class.
#'
#' @param class_counts Output from count_ere_gene_presence_by_class().
#'
#' @return A ggplot object.
#' @export
plot_class_counts <- function(class_counts) {
  .require_pkgs("ggplot2")
  
  ggplot2::ggplot(
    class_counts,
    ggplot2::aes(
      x = match_type,
      y = n_genes_with_hit
    )
  ) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(
      title = "Genes with ERE hits by class",
      x = "ERE match class",
      y = "Genes with at least one hit"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
}


#' Plot GOI vs universe comparison
#'
#' Compare proportions of genes with ERE hits between
#' genes of interest and the full universe.
#'
#' @param summary_goi Output from summarise_ere_hits_presence() for GOI.
#' @param summary_universe Output from summarise_ere_hits_presence() for universe.
#'
#' @return A ggplot object.
#' @export
plot_goi_vs_universe <- function(summary_goi, summary_universe) {
  .require_pkgs(c("ggplot2", "tidyr"))
  
  df <- data.frame(
    group = c("GOI", "Universe"),
    anymatch = c(
      summary_goi$prop_genes_with_anymatch,
      summary_universe$prop_genes_with_anymatch
    ),
    canonical = c(
      summary_goi$prop_genes_with_canonical_ere,
      summary_universe$prop_genes_with_canonical_ere
    )
  )
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = c(anymatch, canonical),
    names_to = "metric",
    values_to = "proportion"
  )
  
  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = group,
      y = proportion,
      fill = metric
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "GOI vs Universe ERE enrichment",
      x = NULL,
      y = "Proportion of genes"
    ) +
    ggplot2::theme_bw()
}


#' Plot observed vs permutation statistics
#'
#' Scatter plot comparing observed statistic against permutation values.
#'
#' @param observed Observed statistic.
#' @param perm Numeric vector of permutation statistics.
#' @param ylab Y-axis label.
#'
#' @return A ggplot object.
#' @export
plot_observed_vs_perm <- function(
    observed,
    perm,
    ylab = "Statistic"
) {
  .require_pkgs("ggplot2")
  
  perm <- perm[!is.na(perm)]
  
  df <- data.frame(
    run = seq_along(perm),
    value = perm
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = run, y = value)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(
      yintercept = observed,
      color = "red",
      linewidth = 1.2
    ) +
    ggplot2::labs(
      title = "Observed statistic vs permutations",
      x = "Permutation",
      y = ylab
    ) +
    ggplot2::theme_bw()
}


#' Plot top genes by anymatch score
#'
#' Summarize total anymatch score per gene and display the top genes.
#'
#' @param df EREfinder results data.frame.
#' @param top_n Number of top genes to display.
#'
#' @return A ggplot object.
#' @export
plot_top_genes_anymatch <- function(df, top_n = 20) {
  .require_pkgs(c("ggplot2", "dplyr"))
  
  gene_df <- df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      total_anymatch = sum(anymatch, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(total_anymatch)) |>
    dplyr::slice_head(n = top_n)
  
  gene_df$gene <- factor(gene_df$gene, levels = rev(gene_df$gene))
  
  ggplot2::ggplot(
    gene_df,
    ggplot2::aes(
      x = gene,
      y = total_anymatch
    )
  ) +
    ggplot2::geom_col(fill = "darkorange") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", top_n, " genes by anymatch score"),
      x = NULL,
      y = "Total anymatch"
    ) +
    ggplot2::theme_bw()
}


# Example usage:
#
# p1 <- plot_perm_distribution(
#   observed = res$summaries$genes_of_interest$n_genes_with_anymatch,
#   perm = res$summaries$permutations$n_genes_with_anymatch,
#   title = "Genes with any ERE hit"
# )
#
# p2 <- plot_class_counts(
#   res$class_specific$genes_of_interest
# )
#
# p3 <- plot_goi_vs_universe(
#   res$summaries$genes_of_interest,
#   res$summaries$universe
# )
#
# p4 <- plot_observed_vs_perm(
#   observed = res$summaries$genes_of_interest$n_genes_with_anymatch,
#   perm = res$summaries$permutations$n_genes_with_anymatch
# )
#
# p5 <- plot_top_genes_anymatch(
#   res$erefinder$genes_of_interest,
#   top_n = 15
# )

plot_perm_dist_labeled <- function(
    res,
    match_type = "perfect_half_site_1bp_sub",
    color_manual = NULL,
    title = "Observed and empirical ERE enrichment by motif class"
) {
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  
  perm_df_long <- res$class_specific$permutations
  obs_df_long <- res$class_specific$pvalues
  
  perm_plot <- perm_df_long |>
    dplyr::filter(match_type == !!match_type)
  
  obs_plot <- obs_df_long |>
    dplyr::filter(match_type == !!match_type) |>
    dplyr::mutate(
      label = paste0("p = ", signif(p_value, 3))
    )
  
  p <- ggplot2::ggplot(
    perm_plot,
    ggplot2::aes(x = match_type, y = n_genes_with_hit)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      ggplot2::aes(fill = match_type)
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      alpha = 0.1,
      size = 0.5
    ) +
    ggplot2::geom_point(
      data = obs_plot,
      ggplot2::aes(x = match_type, y = n_genes_with_hit),
      size = 3,
      color = "red",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = obs_plot,
      ggplot2::aes(
        x = match_type,
        y = n_genes_with_hit,
        label = label
      ),
      vjust = -2,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = match_type,
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
  
  if (!is.null(color_manual)) {
    p <- p + ggplot2::scale_fill_manual(values = color_manual)
  }
  
  p
}


plot_top_genes_tile <- function(
    df,
    match_type = "anymatch",
    top_n = 20,
    show_value = TRUE
) {
  
  requireNamespace("ggplot2")
  requireNamespace("dplyr")
  
  valid_cols <- c(
    "anymatch",
    "Perfect_canonical_ERE",
    "perfect_half_site_1bp_sub",
    "perfect_half_site_2bp_sub",
    "perfect_half_site_3bp_sub"
  )
  
  if (!match_type %in% valid_cols) {
    stop(
      "match_type must be one of: ",
      paste(valid_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  gene_df <- df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      n_hits = sum(.data[[match_type]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_hits)) |>
    dplyr::slice_head(n = top_n)
  
  gene_df$gene <- factor(
    gene_df$gene,
    levels = rev(gene_df$gene)
  )
  
  gene_df$match_type <- match_type
  
  p <- ggplot2::ggplot(
    gene_df,
    ggplot2::aes(
      x = match_type,
      y = gene,
      fill = n_hits
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    
    ggplot2::labs(
      title = paste0(
        "Top ",
        top_n,
        " genes by ",
        match_type
      ),
      x = NULL,
      y = NULL,
      fill = "Hits"
    ) +
    
    ggplot2::theme_bw()
  
  if (show_value) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = n_hits),
        size = 3
      )
  }
  
  p
}
