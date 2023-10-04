#!/usr/bin/env Rscript

source(here::here("src/shared/R-helpers/get_settings.R"))

suppressPackageStartupMessages({
  library("digest")
  library("readr")
  library("httr")
  library("purrr")
  library("dplyr")
  library("ggplot2")
  library("igraph")
  library("tibble")
})

settings <- get_settings(omic = "proteomics", stage = "string_ppi", dataset_name = NULL)

dir.create(settings[["output_directory"]], recursive = TRUE, showWarnings = FALSE)

message("Loading differential expression of proteins...")
de_result_thresholds <- readRDS(settings[["de_result_thresholds_rds"]])

message("Reading UniProt to string id mapping...")
# https://string-db.org/mapping_files/uniprot/mouse.uniprot_2_string.2018.tsv.gz
uniprot2string <- readr::read_tsv(
  settings[["uniprot_to_string"]],
  col_names = c("NCBI_Organism", "UniProtAcc_UniProtName", "string_id", "number1", "number2"),
  col_types = cols(
    `NCBI_Organism` = col_character(),
    `UniProtAcc_UniProtName` = col_character(),
    `string_id` = col_character(),
    `number1` = col_double(),
    `number2` = col_double()
  )
) %>%
  tidyr::separate(
    col = "UniProtAcc_UniProtName",
    into = c("UniProtAcc", "UniProtName"),
    sep = "\\|"
  )


get_network <- function(string_ids, species = 9606, format = "image") {
  hash_value <- digest::digest(list(string_ids, species, format))
  hash_dir <- file.path(settings[["cache"]], "_string_db_cache")
  dir.create(hash_dir, showWarnings = FALSE, recursive = TRUE)
  hash_fn <- file.path(hash_dir, paste0(hash_value, ".rds"))
  if (file.exists(hash_fn)) {
    out <- readRDS(hash_fn)
  } else {
    Sys.sleep(1)
    out <- httr::POST(
      url = sprintf("https://string-db.org/api/%s/network", format),
      body = list(
        "identifiers" = paste(string_ids, collapse ="%0d"), # your protein
        "species" = species, # species NCBI identifier: e.g. 10090 mouse, 9606 human
        #"add_white_nodes" =  15, # add 15 white nodes to my protein
        "network_flavor" =  "confidence", # show confidence links
        "caller_identity" = "sergio.oller@zmnh.uni-hamburg.de" # your app name
      )
    )
    saveRDS(out, hash_fn)
  }
  out
}


message("Getting string images from string web service...")
uniprot_to_string <- uniprot2string %>% dplyr::select(UniProtName, string_id) %>% tibble::deframe()
stringdb_plots_dir <- file.path(settings[["output_directory"]], "stringdb_webplots")
dir.create(stringdb_plots_dir, showWarnings = FALSE, recursive = TRUE)
purrr::imap(de_result_thresholds, function(analysis, analysis_name) {
  message(sprintf(" - %s", analysis_name))
  purrr::imap(analysis, function(result, result_name) {
    message(sprintf("   + %s", result_name))
    purrr::imap(result, function(criterion, crit_name) {
      for (regulation in c("upregulated", "downregulated")) {
        message(sprintf("    + %s (%s)", crit_name, regulation))
        if (regulation == "upregulated") {
          uniprot_names <- criterion$UniProtName[criterion$log2FoldChange > 0]
        } else if (regulation == "downregulated") {
          uniprot_names <- criterion$UniProtName[criterion$log2FoldChange < 0]
        } else {
          stop("typo")
        }
        if (length(uniprot_names) == 0) {
          next
        }
        string_ids <- as.character(na.omit(uniprot_to_string[uniprot_names]))
        resp <- get_network(string_ids, species = settings[["ncbi_identifier"]])
        writeBin(
          resp$content,
          file.path(
            stringdb_plots_dir,
            sprintf("string-net_%s_%s_%s_%s.png",
                    analysis_name,
                    result_name,
                    crit_name,
                    regulation
            )
          )
        )
      }
    })
  })
})



message("Reading protein-protein interaction links...")
get_string_links <- function(cache_dirname, ppi_file, combined_score_threshold) {
  hash_value <- digest::digest(list(ppi_file, combined_score_threshold))
  dir.create(cache_dirname, showWarnings = FALSE, recursive = TRUE)
  cache_fn <- file.path(cache_dirname, paste0("ppi_links_", hash_value, ".rds"))
  if (file.exists(cache_fn)) {
    string_links <- readRDS(cache_fn)
  } else {
    string_links <- readr::read_delim(ppi_file,delim = " ") # 9606.protein.links.v11.0.txt.gz or "10090.protein.links.v11.0.txt.gz
    string_links <- string_links %>% dplyr::filter(combined_score > !!combined_score_threshold)
    saveRDS(string_links, cache_fn)
  }
  string_links
}
# https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens
string_links <- get_string_links(
  cache_dirname = file.path(settings[["cache"]], ".ppi_cache"),
  ppi_file = settings[["string_ppi_file"]],
  combined_score_threshold = settings[["string_combined_score_threshold"]]
)


cleanup_prot_names <- function(x) {
  purrr::map(x, function(prot_name) {
    if (endsWith(prot_name, "_MOUSE")) {
      stringr::str_to_title(gsub("_MOUSE", "", prot_name))
    } else if (endsWith(prot_name, "_HUMAN")) {
      stringr::str_to_upper(gsub("_HUMAN", "", prot_name))
    } else {
      prot_name
    }
  })
}

message("Generating ggraph networks...")
stringdb_ggraph_plots_dir <- file.path(settings[["output_directory"]], "stringdb_ggraph")
dir.create(stringdb_ggraph_plots_dir, showWarnings = FALSE, recursive = TRUE)

purrr::imap(de_result_thresholds, function(analysis, analysis_name) {
  purrr::imap(analysis, function(result, result_name) {
    purrr::imap(result, function(criterion, crit_name) {
      uniprot_names <- criterion$UniProtName
      if (length(uniprot_names) == 0) {
        return(NULL)
      }
      string_ids <- as.character(na.omit(uniprot_to_string[uniprot_names]))
      edges <- string_links %>%
        dplyr::select(protein1, protein2) %>%
        dplyr::filter(
          protein1 %in% string_ids,
          protein2 %in% string_ids
        )
      if (ncol(edges) == 0) {
        edges <- data.frame(protein1 = character(0L), protein2 = character(0L))
      }

      vertices <- tibble::tibble(string_id = string_ids)
      vertices <- dplyr::left_join(
        vertices,
        uniprot2string %>%
          dplyr::select(string_id, UniProtName),
        by = "string_id"
      )
      vertices <- dplyr::left_join(
        vertices,
        criterion,
        by = "UniProtName"
      )
      if (ncol(vertices) == 0) {
        vertices <- data.frame(
          string_id = character(0L),
          UniProtName = character(0L),
          log2foldChange=numeric(0L)
        )
      }

      print(head(edges))
      print(head(vertices))

      ppi_graph <- igraph::graph_from_data_frame(
        edges,
        directed = FALSE,
        vertices = vertices
      )
      saveRDS(
        ppi_graph,
        file = file.path(
          stringdb_ggraph_plots_dir,
          sprintf("%s-%s-%s_ppi_igraph.rds", analysis_name, result_name, crit_name)
        )
      )

      if (nrow(vertices) == 0) {
        return(NULL)
      }

      gplt <- ggraph::ggraph(ppi_graph, layout = 'igraph', algorithm = "kk") +
        ggraph::geom_edge_link(color = "gray") +
        ggraph::geom_node_point(aes(color = log2FoldChange), size = 3) +
        ggraph::geom_node_text(
          aes(label = cleanup_prot_names(UniProtName)),
          repel = TRUE
        ) +
        ggplot2::scale_color_gradient2(low = "darkgreen", mid = "gray", high = "red", midpoint = 0) +
        ggplot2::labs(
          title= "STRING Protein-Protein Interactions",
          subtitle = sprintf("%s - %s - %s" , analysis_name, result_name, crit_name)) +
        ggraph::theme_graph(base_family="sans")

      ggsave(
        plot = gplt,
        filename = file.path(
          stringdb_ggraph_plots_dir,
          sprintf("%s-%s-%s_ggraph_kk.png", analysis_name, result_name, crit_name))
      )

      gplt <- ggraph::ggraph(ppi_graph, layout = 'igraph', algorithm = "nicely") +
        ggraph::geom_edge_link(color = "gray") +
        ggraph::geom_node_point(aes(color = log2FoldChange), size = 3) +
        ggraph::geom_node_text(
          aes(label = cleanup_prot_names(UniProtName)),
          repel = TRUE
        ) +
        ggplot2::scale_color_gradient2(low = "darkgreen", mid = "gray", high = "red", midpoint = 0) +
        ggplot2::labs(
          title= "STRING Protein-Protein Interactions",
          subtitle = sprintf("%s - %s - %s" , analysis_name, result_name, crit_name)) +
        ggraph::theme_graph(base_family="sans")

      saveRDS(
        gplt,
        file = file.path(
          stringdb_ggraph_plots_dir,
          sprintf("%s-%s-%s_ggraph_nicely.rds", analysis_name, result_name, crit_name)
        )
      )

      ggsave(
        plot = gplt,
        filename = file.path(
          stringdb_ggraph_plots_dir,
          sprintf("%s-%s-%s_ggraph_nicely.png", analysis_name, result_name, crit_name)
        )
      )
    })
  })
})


