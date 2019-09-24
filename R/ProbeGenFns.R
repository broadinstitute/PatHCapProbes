#' Function to generate probes for PatH-Cap tool
#'
#' This is the main function to generate probes for hybrid-capture.
#' @param fasta_file Fasta file for the bacterial genome sequence
#' @param gff_file Gff file for the bacterial genome annotation
#' @param probes_file Path to text file where probes would be written
#' @param add_apt [Optional] Whether adapters would be added or not
#' @param pre_str [Optional] Sequence of prefix_adapter, required when add_apt is set to TRUE
#' @param suf_str [Optional] Sequence of suffix_adapter, required when add_apt is set to TRUE
#' @export
#' @examples
#' probeGen()

probeGen <- function(fasta_file, gff_file, probes_file, add_adapt = FALSE, pre_str = "", suf_str = "") {

    # Read the fasta file
    fasta_obj <- rtracklayer::import(fasta_file, format = "fasta")

    # read the gff file
    gff_tab1 <- readGFF(gff_file, tags=c("ID", "Parent", "Name", "locus_tag", "old_locus_tag"))
    gff_tab2 <- gff_tab1[gff_tab1$type == "CDS", ]
    gff_tab3 <- gff_tab2[, c("type", "start", "end", "strand", "Name", "locus_tag")]

    lval1 <- lapply(1:gff_tab3@nrows, get_probe_pos, ltab = gff_tab3, max_len = fasta_obj[[1]]@length)
    probe_lst <- unlist(lval1, recursive = FALSE)

    # Now get the sense and as sequences
    # start with probes under a single 

    file.create(probes_file)
    for (cds_name in names(probe_lst)) {
        #print(cds_name)
        probe_str_set <- probe_lst[[cds_name]]
        get_probe_sequences(cds_name, probe_str_set, fasta_obj[[1]], probes_file, add_adapt, pre_str, suf_str)
    }
} 
