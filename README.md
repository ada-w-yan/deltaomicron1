# deltaomicron1
Model simulations for *Insert paper name here*


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
                 camostat_vec = c(FALSE, TRUE, FALSE),
                 amphoB_vec = c(FALSE, FALSE, TRUE),
                 Calu3 = TRUE,
                 PCR = TRUE,
                 mvr = FALSE,
                 length_run = 2,
                 run_flag = TRUE)

```

fits the model to the Omicron Calu-3 data without drugs and with one drug, and saves the results in results/.  To reproduce the full study, see vignettes.