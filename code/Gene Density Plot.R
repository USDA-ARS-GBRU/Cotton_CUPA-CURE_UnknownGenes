here::i_am("Figure 1 Script.R")

library(rtracklayer)
library(tidyverse)

args <- list("hirsutum" = c(gff = here::here("annot", "Ghirsutum_527_v2.1.gene.gff3.gz"),
                            data = here::here("annot", "Ghirsutum_527_v2.1.annotation_info.txt"),
                            unknowns = here::here("annot", "hirsutumoutputs2.csv")),
             "barbadense" = c(gff = here::here("annot", "Gbarbadense_526_v1.1.gene.gff3.gz"),
                              data = here::here("annot", "Gbarbadense_526_v1.1.annotation_info.txt"),
                              unknowns = here::here("annot", "barbadenseoutputs2.csv"))
             )
#####################
# Load the GFF3 file

make_plot <- function(gff_file,
                      data_file,
                      unknowns_file) {
  print(gff_file)
  gff <- import(gff_file, format = "gff3")
  
  data <- readr::read_delim(data_file)
  
  unknowns <- read.csv(unknowns_file)
  
  
  
  # Filter for gene features
  genes <- gff[gff$type == "gene"]
  
  # Extract start, end, and gene name (usually in the 'Name' or 'gene_name' attribute)
  gene_info <- data.frame(
    seqid = seqnames(genes),
    start = start(genes),
    end = end(genes),
    strand = strand(genes),
    gene_name = mcols(genes)$Name 
  )
  #################
  
  
  ##########
  # Assuming you already have `gene_info` from earlier:
  # gene_info should have columns: seqid, start, end, gene_name
  
  # Set bin size
  bin_size <- 10e5  # 1 million bases
  
  # Step 1: Get chromosome lengths (or define manually if unknown)
  # If not available, estimate from gene data
  chr_lengths <- gene_info |>
    group_by(seqid) |>
    summarize(max_pos = max(end))
  
  # Step 2: Create bins for each chromosome
  bins_list <- lapply(1:nrow(chr_lengths), function(i) {
    chr <- chr_lengths$seqid[i]
    chr_len <- chr_lengths$max_pos[i]
    starts <- seq(1, chr_len, by = bin_size)
    ends <- pmin(starts + bin_size - 1, chr_len)
    data.frame(
      seqid = chr,
      bin_start = starts,
      bin_end = ends,
      bin_id = paste0(chr, ":", starts, "-", ends)
    )
  })
  
  bins_df <- do.call(rbind, bins_list)
  
  # Step 3: Assign each gene to a bin based on its start position
  # We'll do a join-style operation to match gene starts to bins
  
  gene_bins <- gene_info %>%
    select(seqid, start, gene_name) %>%
    inner_join(bins_df, by = "seqid") %>%
    filter(start >= bin_start & start <= bin_end)
  
  # Step 4: Count genes per bin
  gene_counts <- gene_bins %>%
    group_by(bin_id, seqid, bin_start, bin_end) %>%
    summarize(gene_count = n(), .groups = "drop")
  
  # View the result
  head(gene_counts)
  
  gene_counts <- gene_counts %>% filter(!str_starts(bin_id, "s"))
  gene_counts <- gene_counts[order(gene_counts$seqid, gene_counts$bin_start), ]
  
  ggplot(gene_counts, aes(x = bin_start / 1e6, y = gene_count, group = seqid)) +
    geom_bar(stat = "identity") +
    # geom_line()+
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) +
    facet_wrap(~seqid, scales = "free_x") +
    xlab("Genomic Position (Mb)") +
    ylab("Gene Count") +
    ggtitle("Gene Counts per 10Mb Bin by Chromosome") +
    theme_minimal()
  ##################
  
  
  #Next up is adding the list of unknown genes in
  
  unknowns <- unknowns[2]
  unknowns <- unknowns %>% filter(str_starts(.[[1]],"G"))
  
  
  ####Okay we have a data frame that contains all of the unknown genes we need to repeat what we did earlier but instead use this list and not the entire genome

  unknown_gene_info <- gene_info %>%
    filter(gene_name %in% unknowns[[1]])
  
  unknown_gene_bins <- unknown_gene_info %>%
    inner_join(bins_df, by = "seqid") %>%
    filter(start >= bin_start & start <= bin_end)
  
  unknown_counts <- unknown_gene_bins %>%
    group_by(bin_id, seqid, bin_start, bin_end) %>%
    summarize(unknown_count = n(), .groups = "drop")
  
  
  plot_data <- gene_counts %>%
    left_join(unknown_counts, by = c("bin_id", "seqid", "bin_start", "bin_end")) %>%
    mutate(unknown_count = replace_na(unknown_count, 0))
  
  #labels <- as.data.frame(c("G. barbadense", "G.hirsutum"))
  g.sub <- ggplot(plot_data |> dplyr::filter(seqid %in% c("A01", "D01")), 
                  aes(x = bin_start / 1e6)) +
    geom_line(aes(y = gene_count), stat = "identity", fill = "gray70") +
    geom_line(aes(y = unknown_count, group = seqid), color = "blue") +
    ggh4x::facet_grid2(~seqid, scales = "free_x", space = "free_x") +
   xlab("") +
    theme_minimal()
      
  
  g.sub
  
  g.full <- ggplot(plot_data, aes(x = bin_start / 1e6)) +
    geom_line(aes(y = gene_count), stat = "identity", fill = "gray70") +
    geom_line(aes(y = unknown_count, group = seqid), color = "blue") +
    facet_wrap(~seqid, scales = "free_x") +
    xlab("Genomic Position (Mb)") +
    ylab("Gene Count") +
    theme_minimal()
  
  return(list("small" = g.sub,
              "big" = g.full))
}


plots <- lapply(args,
                \(x) make_plot(gff_file = x["gff"],
                               data_file = x["data"],
                               unknowns_file = x["unknowns"]))

done <- ggpubr::ggarrange(plots$barbadense[[1]]  + ggplot2::ylab("G. barbadense\n(Gene Count)"),
                          plots$hirsutum[[1]] +ggplot2::ylab("G. hirsutum\n(Gene Count)") + ggplot2::xlab("Genomic Position (Mb)"), ncol = 1, nrow = 2, labels = c("B", "C"))

ggsave("Figure1.png", done, scale = 1, dpi = 600)

#ggsave("Supplemental_barbadense", plots$barbadense[[2]], scale = 1, dpi = 600)

                    