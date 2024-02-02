# Micasense RedEdge-MX DUAL processing
Simon Oiry

This workflow is an adaptation of the micasense workflow to process
manually images coming from the micasense RedEdge-MX Dual camera. The
original workflow can be found
here\[https://github.com/micasense/imageprocessing\] and has originally
been written in Python.

## Packages

``` r
require(tidyverse)
```

    Le chargement a nécessité le package : tidyverse

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
