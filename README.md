# deltaomicron1
Model simulations for *The altered entry pathway and antigenic distance of the SARS-CoV-2 Omicron variant map to separate domains of spike protein*


## Installation

In R, type

```r
devtools::install_github("ada-w-yan/deltaomicron1")
```

## Usage

As an example,

```r
run_exp_camostat_amphoB("results/",
                 strain = "omicron",
                 camostat_vec = c(FALSE, TRUE),
                 amphoB_vec = c(FALSE, FALSE),
                 Calu3 = TRUE,
                 PCR = TRUE,
                 mvr = FALSE,
                 length_run = 2,
                 run_flag = TRUE)

```

fits the model to the Omicron hNEC data without drugs and with Camostat, and saves the results in results/.  To reproduce the full study, see vignettes.