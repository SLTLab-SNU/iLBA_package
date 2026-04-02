# iLBA

**iLBA** is an R package for the confidential dissemination of aggregated frequency tables from microdata with hierarchical variables. It implements the **Information-Loss-Bounded Aggregation (iLBA)** algorithm described in *[Disseminating massive frequency tables by masking aggregated cell frequencies](https://link.springer.com/article/10.1007/s42952-023-00248-x)*, together with **Small Cell Adjustment (SCA)** at the finest table level.

The package is primarily intended for **statistical agencies** and other producers of official statistics. It focuses on the **confidential dissemination of aggregated frequency tables** from microdata with hierarchical variables. In this setting, two naive approaches are unsatisfactory. Applying **SCA directly to aggregated cell counts** often results in **disclosure risk** through differencing-based inference. On the other hand, **summing SCA-masked counts from the finest level** can lead to substantial **information loss**. The iLBA algorithm is designed to address this problem by protecting aggregated counts while bounding information loss.

## Installation

You can install the package from GitHub:

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

## Main Functions and Example Usage

The **Census dataset** is a synthetic dataset designed to mimic the statistical distribution of the **2010 Korean Census microdata**. It is provided to demonstrate the package's capabilities in handling large-scale administrative and demographic data.

### Dataset Overview
The dataset contains **1,000,000 records** with a structured combination of geographic and demographic attributes:

* **Hierarchical Administrative Variables (4):**
    * `LA1`, `LA2`, `LA3`, and `OA` (Output Area).
* **Key Demographic Variables (5):**
    * `gender`: Sex of the individual.
    * `age`: Age or age group.
    * `edu`: Educational attainment.
    * `mar`: Marital status.
    * `htype`: Household type.

### Metadata
Detailed information regarding variable definitions, code categories, and dataset attributes can be found in the metadata file:

[**census_metadata.xlsx**](https://github.com/SLTLab-SNU/iLBA_package/blob/main/census_metadata.xlsx)

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

Example:

```r
save_full_tb(
  data = census_microdata,
  hkey = c("LA1", "LA2", "LA3", "OA"),
  key = c("gender", "age", "edu", "mar", "htype"),
  mask_thr = 5,
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

Example:

```r
save_agg_tb(
  hkey_level = 3,
  key = c("gender", "age", "htype"),
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

Example:

```r
get_agg_freq(
  hkey_level = 3,
  key = c("gender", "age", "htype"),
  hkey_value = c('01','0104','010407'),
  key_value = c(2, 4, 6),
  input_path = "full_tb.rds"
)
```


## Reference

Park, M.-J., Kim, H. J., and Kwon, S. (2024). *[Disseminating massive frequency tables by masking aggregated cell frequencies](https://doi.org/10.1007/s42952-023-00248-x)*. *Journal of the Korean Statistical Society*, 53, 328–348.
