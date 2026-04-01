# iLBA

**iLBA** is an R package for confidentially disseminating aggregated frequency tables from microdata with hierarchical variables. It implements the **Information-Loss-Bounded Aggregation (iLBA)** algorithm of *[Park et al. (2024)](https://link.springer.com/article/10.1007/s42952-023-00248-x)*, together with **Small Cell Adjustment (SCA)** at the finest table level.

The package is designed for settings in which users may request many frequency tables across different combinations of hierarchical key variables and non-hierarchical key variables. In such settings, masking small cells in each table separately may still leave disclosure risks through differencing across released tables. The iLBA algorithm addresses this problem by masking aggregated counts while bounding information loss.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("SLTLab-SNU/iLBA_package")
```

## Features

- Build the **finest-level frequency table** from microdata
- Apply **SCA masking** to small cells at the finest table level
- Generate **confidential aggregated tables** at coarser hierarchical levels using iLBA
- Export **information loss summaries**
- Return a **masked frequency for a single queried cell**

## Main Functions

### `save_full_tb()`

Constructs the finest-level frequency table from microdata and applies SCA masking to small cells.

```r
save_full_tb(
  data,
  hkey,
  key = NULL,
  mask_thr = 5,
  hkey_rank = NULL,
  key_thr = 100,
  output_path = "full_tb.rds"
)
```

### `save_agg_tb()`

Generates a masked aggregated table at a specified hierarchical level using iLBA.

```r
save_agg_tb(
  hkey_level,
  key,
  input_path = "full_tb.rds",
  output_tb_path = "agg_tb.csv",
  output_iL_path = "info_loss.csv"
)
```

### `get_agg_freq()`

Returns a masked frequency for a single user-specified aggregated cell.

```r
get_agg_freq(
  hkey_level,
  key,
  hkey_value,
  key_value,
  input_path = "full_tb.rds"
)
```


## Reference

Park, M.-J., Kim, H. J., and Kwon, S. (2024). *[Disseminating massive frequency tables by masking aggregated cell frequencies](https://doi.org/10.1007/s42952-023-00248-x)*. *Journal of the Korean Statistical Society*, 53, 328–348.