#' get name of directory which package is in
#'
#' @return a string
get_package_dir <- function() {
  "~/git_repos/deltaomicron/"
}

#' read plaque assay data
#'
#' @import tidyr
#' @import dplyr
#' @return a tibble
read_plaque <- function() {
  filename <- paste0(get_package_dir(), "data/raw_data/Copy of Plaque assays.xlsx")
  a <- xlsx::read.xlsx(filename, 1)
  a <- a[5:11, 6:ncol(a)]
  a <- t(a)
  colnames(a) <- a[1,]
  colnames(a)[1:2] <- c("drug", "strain")
  a <- a[-1,]
  for(i in seq_len(nrow(a))) {
    for(j in seq_len(2)) {
      if(is.na(a[i,j])) {
        a [i,j] <- a[i-1, j]
      }
    }
  }

  supernatant_vol <- .2
  a <- as_tibble(a) %>%
    mutate(across(!c("drug", "strain"), as.numeric)) %>%
    mutate(replicate = rep(seq_len(3), n() / 3)) %>%
    pivot_longer(-c("drug", "strain", "replicate"), names_to = "t", values_to = "V") %>%
    filter(!is.na(V)) %>%
    mutate(t = gsub("inoc", "-1", t),
           t = as.numeric(t),
           drug = gsub("Amphotericin", "AmphoB", drug))
  a
}

#' read qPCR data
#'
#' @import tidyr
#' @import dplyr
#' @return a tibble
read_PCR <- function() {
  filename <- paste0(get_package_dir(), "data/raw_data/Omicron replication data compiled.csv")
  data_tibble <- readr::read_csv(filename)
  data_tibble <- data_tibble[18:nrow(data_tibble),2:ncol(data_tibble)]
  n_replicates <- 3
  strains <- c("delta", "omicron")
  drugs <- c("no_drug", "Camostat", "AmphoB")
  t1 <- data_tibble[,1] %>% pull() %>% as.numeric
  extract_V_tot <- function(strain, drug, col_idx) {
    V_tot <- data_tibble[,seq(col_idx, col_idx + n_replicates - 1)] %>%
      mutate(across(where(is.character), as.numeric))
    colnames(V_tot) <- seq_len(n_replicates)

    V_tot <- V_tot %>%
      mutate(t = t1) %>%
      pivot_longer(-t, names_to = "replicate", values_to = "V_tot") %>%
      mutate(drug = drug, strain = strain)
    V_tot
  }

  par_grid <- tibble(strain = rep(strains, each = length(drugs)),
                     drug = rep(drugs, length(strains))) %>%
    mutate(col_idx = c(2, 5, 8, 13, 16, 19))
  g <- apply_named_args(par_grid, 1, extract_V_tot) %>%
    bind_rows()
  g

}

#' read plaque assay data from Peacock et al. (2021) Nat Microbiol 10.1038/s41564-021-00908-w
read_plaque_peacock <- function() {
  filename <- paste0(get_package_dir(), "data/raw_data/peacock_2021_nat_microbiol.csv")
  data_tibble <- readr::read_csv(filename)
  data_tibble <- data_tibble[,2:ncol(data_tibble)]
  data_tibble <- t(data_tibble)
  colnames(data_tibble) <- data_tibble[1,]
  colnames(data_tibble)[1] <- "strain_drug"
  data_tibble <- data_tibble[-1,]
  for(i in seq(2, nrow(data_tibble))) {
    if(is.na(data_tibble[i, 1])) {
      data_tibble[i,1] <- data_tibble[i-1, 1]
    }
  }

  get_strain_drug <- function(strain_drug, strain) {
    strain_drug <- strsplit(strain_drug,  " + ", fixed = TRUE)
    if(strain) {
      vcapply(strain_drug, function(x) x[1])
    } else {
      vcapply(strain_drug, function(x) x[2])
    }
  }
  data_tibble %>%
    as_tibble %>%
    mutate(strain = get_strain_drug(strain_drug, strain = TRUE),
           strain = recode(strain, ΔCS = "DeltaCS"),
           drug = get_strain_drug(strain_drug, strain = FALSE),
           drug = recode(drug, vehicle = "no drug", Camo = "Camostat")) %>%
    select(-strain_drug) %>%
    pivot_longer(-c("strain", "drug"), names_to = "t", values_to = "V") %>%
    mutate(V = recode(V, `<3` = "1"),
           t = as.numeric(t),
           V = as.numeric(V))

}



#' read Calu 3 data with different drugs
#'
read_Calu3_data <- function(strain, from_raw = FALSE, PCR = TRUE) {
  if(PCR) {
    filename <- paste0(get_package_dir(), "data/PCR_", strain, "_Calu3.csv")
  } else {
    filename <- paste0(get_package_dir(), "data/plaque_", strain, "_Calu3.csv")
  }
  if(!from_raw) {
    data_df <- readr::read_csv(filename) %>%
      as.data.frame()
    colnames(data_df) <- gsub("_no_drug", "", colnames(data_df))
    return(data_df)
  }

  n_replicates <- 3
  inoculum_vol <- .4 # to check
  # supernatant_vol <- 1

  # two-drug combination is in a different file
  read_inner_double <- function(strain, PCR) {
    filename <- paste0(get_package_dir(), "data/raw_data/RG viruses Calu3 E gene copies for Ada.csv")
    if(strain == "omicron") {
      col_select <- c(1, seq(2, 4))
    } else {
      col_select <- c(1, seq(8, 10))
    }
    if(PCR) {
      data_df <- readr::read_csv(filename, skip = 1, n_max = 6, col_select = col_select, col_types = "d")
      values_to <- "V_tot_Camostat_AmphoB"
    } else {
      data_df <- readr::read_csv(filename, skip = 19, col_select = col_select, col_types = "d")
      data_df <- data_df[c(1,3),]
      data_df[2,1] <- 24
      values_to <- "V_Camostat_AmphoB"
    }
    data_df[1,1] <- -1
    data_df[1,2] <- data_df[1,2] * inoculum_vol
    data_df[1,seq(3,4)] <- data_df[1,2]
    colnames(data_df) <- c("t", seq_len(n_replicates))
    data_df <- data_df %>%
      pivot_longer(-t, names_to = "replicate", values_to = values_to) %>%
      mutate(replicate = as.numeric(replicate))
    data_df
  }

  read_inner <- function(strain, PCR) {
    assay_str <- ifelse(PCR, "qPCR", "Plaque")
    filename <- paste0(get_package_dir(), "data/raw_data/", assay_str, "_Calu3_", strain, ".csv")
    V_str <- ifelse(PCR, "V_tot_", "V_")
    data_df <-readr::read_csv(filename) %>%
      select(where(is.numeric)) %>%
      rename(t = "...1") %>%
      pivot_longer(-t, names_to = "drug", values_to = "V") %>%
      mutate(replicate = rep(seq_len(n_replicates), n() / n_replicates),
             drug = vcapply(drug, function(x) strsplit(x, "...", fixed = TRUE)[[1]][1])) %>%
      pivot_wider(names_from = drug, names_prefix = V_str, values_from = V)
    if(!PCR) {
      # hand coded inoculum

      inoculum <- ifelse(strain == "omicron", 1.5e4, 2.5e4) * inoculum_vol

      data_df_inoculum <- data_df %>%
        slice(seq_len(n_replicates)) %>%
        mutate(across(starts_with("V_"), function(x) inoculum),
               t = -1)
      data_df <- rbind(data_df, data_df_inoculum)
    }
    data_df
  }

  if(PCR) {
    data_df <- lapply(c(FALSE, TRUE), read_inner, strain = strain) %>%
      do.call(full_join, .)
    data_df_double <- lapply(c(FALSE, TRUE), read_inner_double, strain = strain) %>%
      do.call(full_join, .)
  } else {
    data_df <- read_inner(strain, PCR = FALSE)
    data_df_double <- read_inner_double(strain, PCR = FALSE)
  }
  data_df <- data_df %>%
    full_join(data_df_double) %>%
    rename_with(function(x) gsub("_No drug", "", x), everything())
  readr::write_csv(data_df, filename)
  data_df
}
#' read data with and without camostat
#'
#' @param strain WT, DeltaCS, delta or omicron
#' @param from_raw logical, whether to read from raw data
#' @param PCR logical, whetehr to include PCR data
#' @import tidyr
#' @import dplyr
#' @return a data frame
read_camostat_data <- function(strain1, from_raw = FALSE, PCR = TRUE) {
  if(!from_raw) {
    if(PCR) {
      filename <- paste0(get_package_dir(), "data/PCR_", strain1, ".csv")
    } else {
      filename <- paste0(get_package_dir(), "data/plaque_", strain1, ".csv")
    }

    data_df <- readr::read_csv(filename) %>%
      as.data.frame()
    colnames(data_df) <- gsub("_no_drug", "", colnames(data_df))
    return(data_df)
  }

  if(strain1 %in% c("WT", "DeltaCS") && PCR) {
    stop("PCr data not available")
  }

  supernatant_vol <- .2

  if(strain1 %in% c("delta", "omicron")) {
    data_df <- read_plaque()
  } else {
    data_df <- read_plaque_peacock()
  }

  n_replicates <- 3
  n_timepoints <- length(unique(data_df$t))
  drugs <- c("no drug", "Camostat", "AmphoB")

  strain2 <- switch(strain1,
                    delta = "Delta1",
                    omicron = "Omicron1",
                    strain1)
  data_df <- data_df %>%
    filter(drug %in% drugs,
           strain == strain2) %>%
    select(-strain) %>%
    mutate(drug = gsub(" ", "_", drug),
           V = V * supernatant_vol)# convert from pfu/mL to pfu
  if(!("replicate" %in% colnames(data_df))) {
    data_df$replicate <- rep(rep(seq_len(n_replicates), each = n_timepoints), length(drugs))
  }
  # return(data_df)
  data_df <- data_df %>%
    pivot_wider(names_from = drug, values_from = V, names_prefix = "V_") %>%
    as.data.frame()

  # specify 0 hour timepoints to be below threshold, which we know is likely from PCR --
  # won't need this hacky solution after adding PCR data

  if(PCR) {
    data_df_pcr <- read_PCR() %>%
      filter(strain == strain1) %>%
      select(-strain) %>%
      pivot_wider(names_from = drug,
                  names_prefix = "V_tot_",
                  values_from = V_tot) %>%
      mutate(replicate = as.numeric(replicate))
    data_df <- full_join(data_df, data_df_pcr)
    readr::write_csv(data_df, paste0(get_package_dir(), "data/PCR_", strain1, ".csv"))
  } else {
    extend_df <- data.frame(replicate = seq_len(n_replicates),
                            t = 0,
                            V_no_drug = 1,
                            V_Camostat = 1)

    data_df <- rbind(data_df, extend_df)
    # hardcode peacock et al. inoculum for now

    if(!(strain1 %in% c("Delta1", "Omicron1"))) {
      extend_df <- data.frame(replicate = seq_len(n_replicates),
                              t = -1,
                              V_no_drug = 2.5e5*.1,
                              V_Camostat = 2.5e5*.1)
      data_df <- rbind(data_df, extend_df)
    }
  }

  colnames(data_df) <- gsub("_no_drug", "", colnames(data_df))

  data_df

}
