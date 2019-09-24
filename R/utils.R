get_probe_pos <- function(lindex, ltab, max_len) {

    larr <- ltab[lindex, ]
    ltype <- larr[["type"]]
    start <- larr[["start"]]
    end <- larr[["end"]]
    lstrand <- larr[["strand"]]
    lname <- larr[["Name"]]
    llocus_tag <- larr[["locus_tag"]]

    # Always start from the 5-prime end, i.e. the start position and keep on adding
    probe_len <- 100
    probe_str_set <- c()

    #if (end - start <= probe_len) {
        #print("Small gene")
        #print(llocus_tag)
    #}

    feature_size <- end - start + 1
    for (lstart in seq(1, feature_size, probe_len)) {
        lend <- lstart + probe_len -1
        ab_start <- start + lstart -1
        ab_end <- start + lend -1
        if (ab_end <= max_len) {
            lprobe_str <- paste0(llocus_tag, "_", ltype, "_", "[", start, ":", 
                end, "]", "(", lstrand, ")", "[", lstart, ":", lend, "]", "___", ab_start, ":", ab_end)
            probe_str_set <- c(probe_str_set, lprobe_str)
            #print(lprobe_str)
        }
    }

    cds_str <- paste(llocus_tag, ltype, start, end, lstrand, sep = "_")    
    
    local_lst <- list()
    local_lst[[cds_str]] <- probe_str_set
    return (local_lst)
}

get_probe_sequences <- function(cds_name, probe_str_set, lfasta_obj1, probes_file, add_adapt, pre_str, suf_str) {

    # get strand
    parts <- strsplit(cds_name, "_")
    lstrand <- parts[[1]][5]

    sense_probe_count <- length(probe_str_set)
    
    as_probe_str_set <- NA
    if (lstrand == "+") {
        as_probe_str_set <- probe_str_set[seq(1, sense_probe_count, 2)]
    } else if (lstrand == "-") {
        as_probe_str_set <- probe_str_set[seq(sense_probe_count, 1, -2)]
    }

    for (lprobe_str in probe_str_set) {
        # Get the probe
        parts <- strsplit(lprobe_str, "___")[[1]]
        prob_head1 <- parts[[1]]
        probe_pos1 <- parts[[2]]
        
        parts2 <- strsplit(probe_pos1, ":")[[1]]

        lstart <- parts2[1]
        lend <- parts2[2]

        # Get the sequence
        lseq <- lfasta_obj1[lstart:lend]
        # Check if the strand is + or -
        lseq_s <- NA
        if (lstrand == '+') {
            lseq_s <- reverseComplement(lseq)
        } else {
            lseq_s <- lseq
        }

        sense_id <- paste0(">", prob_head1, "")

        sense_val1 <- as.character(lseq_s)
        sense_val <- NA
        if (add_adapt) {
            sense_val <- paste0(pre_str, sense_val1, suf_str)
        } else {
            sense_val <- sense_val1
        }

        #print(sense_id)
        #print(sense_val)
        write(sense_id, file = probes_file, append = TRUE)
        write(sense_val, file = probes_file, append = TRUE)
        
        if (lprobe_str %in% as_probe_str_set) {
            # Print lseq_as
            lseq_as <- reverseComplement(lseq_s)

            anti_sense_id <- paste0(">", prob_head1, "RC")

            anti_sense_val1 <- as.character(lseq_as)
            anti_sense_val <- NA
            if (add_adapt) {
                anti_sense_val <- paste0(pre_str, anti_sense_val1, suf_str)
            } else {
                anti_sense_val <- anti_sense_val1
            }

            #print(anti_sense_id)
            #print(anti_sense_val)
            write(anti_sense_id, file = probes_file, append = TRUE)
            write(anti_sense_val, file = probes_file, append = TRUE)
        }
        
    }
}
