#────────────────────────  PACKAGES
suppressPackageStartupMessages({
  library(phyloseq); library(rsample); library(dplyr); library(purrr)
  library(ALDEx2);    library(ANCOMBC);    library(ggVennDiagram)
  library(RColorBrewer);    library(dplyr);    library(grid)
  
})

suppressPackageStartupMessages(library(stringr))

# ======== Manipulation functions ========


## SUBSET PHYLOSEQ DADO UN SPLIT
# TODO: included in final RMD
subset_ps <- function(ps, split) {
  samp <- analysis(split)$SampleID
  prune_samples(samp, ps)
}


## a simplified function
extracting_genus <- function(x) {
  # Removing prefixes including: __, f__, etc., _NA
  x <- gsub("^[a-z]__+", "", x)  # elimina prefijos como g__, f__, o c__, etc.
  
  # Optional: replace "" for NA if the string is empty
  x[x == ""] <- NA
  
  return(x)
}



sanitize_bacteria <- function(v) {
  v <- as.character(v)
  v <- trimws(v)
  v <- extracting_genus(v)   # your function: strips g__/f__/... but keeps value, even *_NA
  v[v == ""] <- NA                  # enforce NA for empty strings
  v
}

build_orientation_long <- function(res_list) {
  imap_dfr(res_list, ~{
    df <- as.data.frame(.x, stringsAsFactors = FALSE)
    if (nrow(df) == 0L) return(tibble(bacteria = character(0), metodo = .y, orientation = integer(0)))
    if (!"bacteria" %in% names(df)) df <- tibble::rownames_to_column(df, "bacteria")
    
    df$bacteria <- sanitize_bacteria(df$bacteria)
    df <- df[!is.na(df$bacteria), , drop = FALSE]
    
    # infer orientation if missing
    if (!"orientation" %in% names(df)) {
      cand <- intersect(names(df), c("log2FoldChange","logFC","effect","diff.btw",
                                     "lda","lda_score","LDA.score","beta","coef"))
      if (length(cand)) {
        sgn <- sign(df[[cand[1]]])
        df$orientation <- ifelse(is.na(sgn), NA_integer_,
                                 ifelse(sgn > 0, 1L, ifelse(sgn < 0, -1L, NA_integer_)))
      } else if ("enrich_group" %in% names(df)) {
        # fallback: 1 -> +1; 0 -> -1; else NA
        df$orientation <- ifelse(is.na(df$enrich_group), NA_integer_,
                                 ifelse(df$enrich_group == 1, 1L,
                                        ifelse(df$enrich_group == 0, -1L, NA_integer_)))
      } else {
        df$orientation <- NA_integer_
      }
    }
    
    df %>%
      transmute(bacteria, metodo = .y, orientation = as.integer(orientation)) %>%
      distinct(bacteria, metodo, .keep_all = TRUE)
  })
}


# Finding oral bacteria by each dataframe
oral_bacteria_finding <- function(df, comparison_name) {
  
  oral_genus <- c("Actinomyces", "Aggregatibacter", "Alloprevotella",
                  "Campylobacter", "Capnocytophaga", "Catonella",
                  "Desulfobulbus", "Dialister", "Filifactor",
                  "Fusobacterium", "Mogibacterium", "Parvimonas",
                  "Peptococcus", "Peptostreptococcus",
                  "Porphyromonas", "Prevotella", "Prevotella_7",
                  "Pseudoramibacter", "Tannerella", "Treponema", "Veillonella")
  
  df %>%
    dplyr::mutate(
      bacteria = as.character(as.vector(bacteria)),             
      enrich_group = as.character(enrich_group)                 
    ) %>%
    dplyr::filter(
      stringr::str_detect(bacteria, "g__"),                     
      stringr::str_remove(bacteria, ".*g__") %in% oral_genus    
    ) %>%
    dplyr::mutate(comparison = comparison_name)                
}




extract_clean_genus <- function(bacteria_vec) {
  if (is.null(bacteria_vec)) return(character(0))
  x <- as.character(bacteria_vec)
  # your taxonomy-to-genus function
  x <- extracting_genus(x)
  
  # keep last level, strip prefixes/suffixes
  x <- sub(".*\\|", "", x)
  x <- str_remove(x, "^[A-Za-z]+__")
  x <- str_remove(x, "_[A-Za-z]+__$")
  
  # drop empty/NA and Incertae_Sedis variants
  x[x == ""] <- NA_character_
  x <- x[!is.na(x)]
  x <- x[!(x %in% c("Incertae_Sedis", "Incertae_Sedis_NA"))]
  
  unique(x)
}

# TODO: check this function
fill_presence <- function(df_presencia, res_dfs_subgroup) {
  df_presencia$bacteria <- sanitize_bacteria(df_presencia$bacteria)
  df_presencia <- df_presencia[!is.na(df_presencia$bacteria), , drop = FALSE]
  
  for (nm in names(res_dfs_subgroup)) {
    df <- res_dfs_subgroup[[nm]]
    if (!is.data.frame(df) || !"bacteria" %in% names(df)) {
      df_presencia[[nm]] <- 0L
      next
    }
    set_bact <- unique(na.omit(df$bacteria))   # <-- critical to avoid NA in %in%
    df_presencia[[nm]] <- as.integer(df_presencia$bacteria %in% set_bact)
  }
  
  rownames(df_presencia) <- df_presencia$bacteria
  df_presencia
}

drop_incertae <- function(df) {
  keep <- !(df$bacteria %in% c("Incertae_Sedis_NA", "Incertae_Sedis"))
  df[keep, , drop = FALSE]
}

assert_no_na <- function(df, df_name) {
  if (anyNA(df)) {
    idx <- which(is.na(df), arr.ind = TRUE)
    bad <- unique(data.frame(
      row_index = idx[,1],
      bacteria  = df$bacteria[idx[,1]],
      column    = colnames(df)[idx[,2]],
      stringsAsFactors = FALSE
    ))
    message(df_name, " contains NAs at:")
    print(bad, row.names = FALSE)
    stop(df_name, " has NA values. Check upstream sanitization / inputs.")
    # To continue instead:
    # df[is.na(df)] <- 0L
    # return(invisible(TRUE))
  }
  invisible(TRUE)
}


# ancombc useful: 
last_informative_rank <- function(x,
                                  uninf = c("uncultured","unidentified","unclassified",
                                            "unassigned","unknown","norank","na")) {
  pick_one <- function(s) {
    if (is.na(s) || s == "") return(s)                       # keep as is
    if (!startsWith(s, "d__Bacteria")) return(s)             # only act on the long form
    m <- str_match_all(s, "_?([dpcfosg])__(.*?)(?=(?:_[dpcfosg]__|$))")[[1]]
    if (nrow(m) <= 1) return(s)                              # only one rank -> leave as is
    
    ranks <- m[, 2]
    names <- m[, 3]
    
    # Walk from deepest to shallowest; skip uninformative names
    for (i in seq.int(length(names), 1)) {
      nm  <- names[i]
      nm2 <- tolower(nm)
      bad <- nm2 == "" ||
        nm2 %in% uninf ||
        startsWith(nm2, "uncultured_")                  # e.g. "uncultured_bacterium"
      if (!bad) return(paste0(ranks[i], "__", nm, "_NA"))
    }
    NA_character_                                           # all were uninformative
  }
  
  vapply(x, pick_one, FUN.VALUE = character(1))
}

# ======= Useful Function to include =====
# Function: remove ASV marker to get the "base" method name
get_base_name <- function(x) {
  x %>%
    # remove ASV surrounded by underscores/dots/dashes
    str_replace(regex("([._-])?ASV([._-])?", ignore_case = TRUE), "\\1") %>%
    # remove multiple trailing separators
    str_replace(regex("[._-]+$", ignore_case = TRUE), "") %>%
    # remove multiple leading separators
    str_replace(regex("^[._-]+", ignore_case = TRUE), "") %>%
    trimws()
}



# Helper: detect ASV id in 'feature' (customizable: you can tune regex if your features have another pattern)
detect_asv <- function(x) {
  # returns the ASV string if it looks like an ASV id; otherwise NA
  is_asv <- grepl("^ASV[0-9A-Za-z_-]*$", x) | grepl("(^|[_\\.-])ASV([_\\.-]|$)", x, ignore.case = TRUE)
  ifelse(is_asv, x, NA_character_)
}



# Compare Genus per base 
compare_base <- function(base_name) {
  # matched dfs
  n_asv   <- base_asv  %>% filter(base == base_name) %>% pull(name)
  n_noasv <- base_noasv %>% filter(base == base_name) %>% pull(name)
  
  # if multiple per base, take union of genus
  genus_asv <- unique(unlist(
    lapply(dfs_asv[n_asv],   function(df) extract_clean_genus(df$bacteria))
  ))
  genus_no  <- unique(unlist(
    lapply(dfs_noasv[n_noasv], function(df) extract_clean_genus(df$bacteria))
  ))
  
  # metrics
  inters <- length(intersect(genus_asv, genus_no))
  jacc   <- if ((length(genus_asv) + length(genus_no) - inters) == 0) 1 else
    inters / (length(genus_asv) + length(genus_no) - inters)
  identical_sets <- setequal(genus_asv, genus_no)
  
  # detailed table
  detail <- tibble(
    base     = base_name,
    Genus    = union(genus_asv, genus_no),
    in_noASV = Genus %in% genus_no,
    in_ASV   = Genus %in% genus_asv
  )
  
  summary <- tibble(
    base              = base_name,
    n_genus_noASV     = length(genus_no),
    n_genus_ASV       = length(genus_asv),
    intersection      = inters,
    jaccard           = jacc,
    identical_sets    = identical_sets
  )
  
  list(summary = summary, detail = detail)
}

get_scale_alpha_values <- function(qval_vector) {
  # Customizable
  #  q < 0.05 → alpha = 1; else → alpha = 0.3
  alphas <- ifelse(qval_vector < 0.05, 1, 0.3)
  names(alphas) <- ifelse(qval_vector < 0.05, "*", "")
  return(alphas)
}

# ======== Plotting  ========

sets_building <- function(df, tools_sub, tipo = "noBH") {
  lapply(tools_sub, function(herr) {
    patron <- if (tipo == "BH") paste0(herr, ".*[^a-zA-Z]BH([^a-zA-Z]|$)")
    else                paste0(herr, ".*noBH")
    cols_herr <- grep(patron, colnames(df), value = TRUE)
    if (length(cols_herr) == 0) return(character(0))
    seleccionadas <- rowSums(df[, cols_herr, drop = FALSE]) > 0
    df$bacteria[seleccionadas]
  }) |> `names<-`(tools_sub)
}

#  2–5 sets Venn Diagram with ggVennDiagram
plot_venn <- function(sets_named, titulo = NULL) {
  sets_filtrados <- sets_named[sapply(sets_named, length) > 0]
  k <- length(sets_filtrados)
  if (k < 2) {
    warning("No hay suficientes conjuntos para Venn.")
    return(invisible(NULL))
  }
  if (k > 5) {
    stop("Para más de 5 conjuntos, usa UpSet (ver opción B).")
  }
  pal <- if (k <= 3) brewer.pal(3, "Set2")[seq_len(k)] else brewer.pal(max(3, k), "Set3")[seq_len(k)]
  p <- ggVennDiagram(sets_filtrados, label_alpha = 0) +
    scale_fill_gradient(low = "#F5F5F5", high = "#9ECAE1") +
    scale_color_manual(values = pal) +
    labs(title = titulo)
  print(p)
  invisible(p)
}

# All combinations
venn_all_combinations <- function(df, tools, tipo = "noBH", k = 4,
                                     save = FALSE, prefijo = "venn") {
  stopifnot(k >= 2, k <= 5)
  sets_all <- sets_building(df, tools, tipo)
  combos <- combn(names(sets_all), k, simplify = FALSE)
  res_plots <- vector("list", length(combos))
  for (i in seq_along(combos)) {
    sub_names <- combos[[i]]
    sub_sets  <- sets_all[sub_names]
    titulo    <- paste0("Venn (", tipo, ") - ", paste(sub_names, collapse = " + "))
    p <- try(plot_venn(sub_sets, titulo), silent = TRUE)
    res_plots[[i]] <- p
    if (save && inherits(p, "gg")) {
      fn <- paste0(prefijo, "_", tipo, "_", paste(sub_names, collapse = "_"), ".png")
      ggsave(fn, plot = p, width = 8, height = 6, dpi = 300)
    }
  }
  invisible(res_plots)
}


# ======== Plotting UpSet =======


# Helper function for creating the upset-compatible data
# Converts sets into a binary presence/absence data frame for UpSet
create_upset_df <- function(df, methods, type = "noBH") {
  sets <- sets_building(df, methods, type)
  all_elements <- unique(unlist(sets))
  
  binary_df <- tibble(bacteria = all_elements)
  for (method in methods) {
    binary_df[[method]] <- binary_df$bacteria %in% sets[[method]]
  }
  binary_df
}

# Plot 4 upset plots: prevalence vs no prevalence, BH vs noBH
compare_upset_4 <- function(df_prev, df_noprev, methods, save = FALSE, filename = "upset_comparison_4.png") {
  types <- c("noBH", "BH")
  titles <- c("No Prevalence (noBH)", "No Prevalence (BH)",
              "Prevalence (noBH)", "Prevalence (BH)")
  
  plots <- vector("list", 4)
  
  for (i in seq_along(titles)) {
    type <- types[(i - 1) %% 2 + 1]
    df <- if (i <= 2) df_noprev else df_prev
    df_upset <- create_upset_df(df, methods, type)
    plots[[i]] <- ComplexUpset::upset(df_upset,
                        intersect = methods,
                        name = "Methods",
                        base_annotations = list('Intersections' = ComplexUpset::intersection_size())) +
      ggtitle(titles[i])
  }
  
  combined_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  
  if (save) {
    ggsave(filename, plot = combined_plot, width = 16, height = 12, dpi = 300)
  }
  
  print(combined_plot)
  invisible(combined_plot)
}

# --- Cleaner for tax labels in your "genus" column ---
# one string
clean_tax_label_one <- function(x) {
  if (is.null(x)) return(NA_character_)
  x <- as.character(x)
  if (length(x) == 0) return(NA_character_)
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  
  x <- str_squish(x)
  
  # split in tokens, drop exact duplicates
  toks <- unlist(str_split(x, "\\s+"))
  toks <- unique(toks)
  
  # pick the longest token per rank (often most informative)
  pick_longest <- function(v) if (length(v)) v[which.max(nchar(v))] else NA_character_
  
  s_tok <- pick_longest(toks[str_detect(toks, "^s__")])
  g_tok <- pick_longest(toks[str_detect(toks, "^g__")])
  f_tok <- pick_longest(toks[str_detect(toks, "^f__")])
  o_tok <- pick_longest(toks[str_detect(toks, "^o__")])
  
  # species considered uninformative
  is_uninformative_s <- function(s) {
    if (is.na(s)) return(FALSE)
    str_detect(s, regex(
      "uncultured|unknown|unclassified|metagenome|ambiguous|organism|^s__sp(_|\\b)|\\bsp\\b",
      ignore_case = TRUE
    ))
  }
  
  strip_na_number <- function(tok) {
    if (is.na(tok) || !nzchar(tok)) return(tok)
    sub("(?<=_NA)_\\d+$", "", tok, perl = TRUE)
  }
  
  out <- NA_character_
  
  if (!is.na(s_tok)) {
    if (is_uninformative_s(s_tok)) {
      # keep g__ + s__ (if genus exists), with inline cleanup
      g_tok <- strip_na_number(g_tok)
      s_tok <- strip_na_number(s_tok)
      out <- str_squish(paste(g_tok, s_tok))
    } else {
      # informative species -> keep species only
      out <- strip_na_number(s_tok)
    }
  } else if (!is.na(g_tok)) {
    out <- strip_na_number(g_tok)
  } else if (!is.na(f_tok)) {
    out <- strip_na_number(f_tok)
  } else if (!is.na(o_tok)) {
    out <- strip_na_number(o_tok)
  } else if (length(toks)) {
    out <- strip_na_number(toks[1])
  }
  
  if (!nzchar(out)) NA_character_ else out
}

# vectorizada (base R)
clean_tax_label <- function(x) {
  vapply(x, clean_tax_label_one, FUN.VALUE = character(1), USE.NAMES = FALSE)
}


# Vectorize option:
clean_tax_label <- function(x) purrr::map_chr(x, clean_tax_label_one)



# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}



## plotting
plot_da_bars <- function(df,
                         alpha = 0.05,
                         remove_untrusted = FALSE,
                         ytitle = "Log2 fold change") {
  
  if (remove_untrusted)
    df <- df[!grepl("Incertae|uncultured", df$taxon_id), ]
  
  df$qval.txt <- cut(df$qval,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "-"),
                     right  = FALSE)
  df$qval.txt <- factor(df$qval.txt,
                        levels = c("***", "**", "*", "-"))
  
  ggplot2::ggplot(df,
                  ggplot2::aes(x = taxon_id, y = lfc,
                               fill = direction, alpha = qval.txt)) +
    ggplot2::geom_bar(stat = "identity",
                      position = ggplot2::position_dodge2(width = 0.9,
                                                          preserve = "single")) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lfc - se,
                                        ymax = lfc + se),
                           width = 0.5,
                           position = ggplot2::position_dodge2(width = 0.9,
                                                               preserve = "single"),
                           colour = "grey36") +
    ggplot2::geom_text(ggplot2::aes(label = qval.txt,
                                    y = lfc + orientation * (abs(lfc) + 0.22)),
                       vjust = 0.7, colour = "grey36") +
    ggplot2::geom_hline(yintercept = 0, alpha = 0.5) +
    ggplot2::facet_wrap(~target_level, scales = "free_y", ncol = 1) +
    ggplot2::scale_alpha_manual(
      values = c("***" = 1, "**" = 0.9, "*" = 0.7, "-" = 0.3)) +
    ## Predefined colours: “Low”vs“Sev”
    ggplot2::scale_fill_manual(values = c("Low" = "#22A884FF",
                                          "Sev" = "#440154FF")) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = ytitle)
}


# Plot ALDEx2
plot_aldex_barplot <- function(df) {
  df$qval.txt <- cut(df$qval,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "-"))
  
  df$qval.txt <- factor(df$qval.txt, levels = c("***", "**", "*", "-"))
  
  ggplot(data = df, aes(x = taxon_id, y = lfc, fill = direction, alpha = qval.txt)) +
    scale_alpha_manual(values = c("***" = 1, "**" = 0.9, "*" = 0.7, "-" = 0.3)) +
    geom_bar(stat = "identity",
             position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.5,
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  color = "gray36") +
    geom_text(aes(label = qval.txt, y = lfc + orientation * (abs(lfc) + 0.2)),
              vjust = 0.7, color = "gray36",
              position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    facet_wrap(~target_level, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = "Log2 fold change") +
    scale_fill_manual(values = c("Low" = "#22A884FF", "Sev" = "#440154FF")) +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10, face = "italic"),
      legend.text = element_text(size = 8)
    )
}



# Plot Linda

linda_plot <- function(df_fig_noBH) {
  # Requires: ggplot2
  # Expects columns: Genus, logFC, direction, alpha_grp, se, qval_txt, orientation, target_level
  
  ggplot(df_fig_noBH,
         aes(x = Genus,
             y = logFC,
             fill = direction,
             alpha = alpha_grp)) +
    
    scale_alpha_manual(values = c("1" = 1, "0.6" = 0.6, "0.3" = 0.3)) +
    
    geom_bar(stat = "identity",
             position = position_dodge2(width = 0.9, preserve = "single")) +
    
    geom_errorbar(aes(ymin = logFC - se, ymax = logFC + se),
                  width = 0.5,
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  color = "gray36") +
    
    geom_text(aes(label = qval_txt,
                  y = logFC + orientation * (se + 0.22)),
              vjust = 0.7, color = "gray36",
              position = position_dodge2(width = 0.9, preserve = "single"),
              size = 3) +
    
    geom_hline(yintercept = 0, alpha = 0.5) +
    
    facet_wrap(~target_level, scales = "free_y", ncol = 1) +
    
    coord_flip() +
    labs(x = NULL, y = "Log2 Fold-Change") +
    theme_bw() +
    theme(
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 10),
      axis.title.y = element_text(size = 10, face = "italic"),
      legend.text  = element_text(size = 8)
    ) +
    scale_color_manual(values = c("#22A884FF", "#440154FF")) +
    scale_fill_manual(values = c("#22A884FF", "#440154FF"))
}



# ZicoSeq Plotting Auxiliar function
plot_logFC_bar <- function(df, title_plot) {
  df <- df %>%
    arrange(logFC) %>%
    mutate(bacteria = factor(bacteria, levels = unique(bacteria)))
  
  ggplot(df, aes(x = bacteria, y = logFC, fill = direction, alpha = qval.txt)) +
    geom_bar(stat = "identity",
             position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_text(aes(label = qval.txt, y = logFC + 0.2 * orientation), 
              vjust = 0.7, color = "gray36",
              position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_hline(yintercept = 0, alpha = 0.5) +
    coord_flip() +
    scale_fill_manual(values = c("Positive LFC" = "#440154FF", 
                                 "Negative LFC" = "#22A884FF")) +
    scale_alpha_manual(values = c("***" = 1, "**" = 0.8, "*" = 0.6, "-" = 0.3)) +
    labs(x = NULL, y = "Log fold change", title = title_plot, 
         fill = "Direction", alpha = "q-value") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10, face = "italic"),
      legend.text = element_text(size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}


# ======== DAA METHODS ========



## ANCOM-BC
### The following code regarding the function run_ancombc developed by Pablo Aja-Macaya.Additional features and modifications were added by Elsa Martin-De Arribas in function run_ancombc_1vs1

run_ancombc_1vs1 <- function(physeq_obj, target_level, target_variables, prv_cut=0, qval_min=0.05, p_adj_method="BH"){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=8,
                 p_adj_method=p_adj_method, alpha = 0.05, conserve = TRUE)
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  pval_long <- res$p_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "pval")
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, qval_long, pval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"
  
  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1
  
  # Add target level
  df_fig$target_level <- target_level
  
  return(df_fig)
}


## ancombc
### ancombc_noprev_noBH_function
set.seed(123)


library(ANCOMBC)
run_ancombc_nPnB <- function(physeq_obj, target_level, target_variables, prv_cut=0, qval_min=0.05){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=8,
                 p_adj_method="none", alpha = 0.05, conserve = TRUE)
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  pval_long <- res$p_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "pval")
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, qval_long, pval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"
  
  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  #merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1
  
  # Add target level
  df_fig$target_level <- target_level
  
  return(df_fig)
}

# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}


### ancombc_noprev_BH_function

library(ANCOMBC)
run_ancombc_nPB <- function(physeq_obj, target_level, target_variables, prv_cut=0, qval_min=0.05){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=8,
                 p_adj_method="BH", alpha = 0.05, conserve = TRUE)
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  pval_long <- res$p_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "pval")
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, qval_long, pval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"
  
  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig_nopre_noBH <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig_nopre_noBH$direction <- factor(df_fig_nopre_noBH$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig_nopre_noBH$qval.txt <- cut(df_fig_nopre_noBH$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig_nopre_noBH$orientation <- 0
  df_fig_nopre_noBH$orientation[df_fig_nopre_noBH$lfc > 0] <- 1
  df_fig_nopre_noBH$orientation[df_fig_nopre_noBH$lfc < 0] <- -1
  
  # Add target level
  df_fig_nopre_noBH$target_level <- target_level
  
  return(df_fig_nopre_noBH)
}

# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}



### ancombc_prev_noBH_function
run_ancombc_PnB_toolprev <- function(physeq_obj, target_level, target_variables, prv_cut=0, qval_min=0.05){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=8,
                 p_adj_method="none", alpha = 0.05, conserve = TRUE)
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  pval_long <- res$p_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "pval")
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, qval_long, pval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"
  
  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1
  
  # Add target level
  df_fig$target_level <- target_level
  
  return(df_fig)
}

# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}



### ancombc_prev_BH_function


run_ancombc_PB <- function(physeq_obj, target_level, target_variables, prv_cut=0, qval_min=0.05){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }
  
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=8,
                 p_adj_method="BH", alpha = 0.05, conserve = TRUE) # conserve = TRUE due to low sample size
  # Get results
  res = out$res
  
  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")
  pval_long <- res$p_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "pval")
  
  
  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, pval_long, qval_long))
  
  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged
  
  # Change column name
  names(merged)[1] <- "taxon_id"
  
  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }
  
  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)
  
  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))
  
  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))
  
  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1
  
  # Add target level
  df_fig$target_level <- target_level
  
  return(df_fig)
}

# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}



# ALDEx2



## DeSeq2





## LefSe




## LinDA





## ZicoSeq


























