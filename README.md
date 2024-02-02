# Micasense RedEdge-MX DUAL processing
Simon Oiry

This workflow is an adaptation of the micasense workflow to process
manually images coming from the micasense RedEdge-MX Dual camera. The
original workflow can be found
[here](https://github.com/micasense/imageprocessing) and has originally
been written in Python.

## Packages

The first thing to do is to ensure that all the packages are ready to be
used. The [exiftoolr](https://github.com/JoshOBrien/exiftoolr) packages
is used to read Exif of tiff files. After installing the package you
will need to run this line of code : `exiftoolr::install_exiftool()`

``` r
exiftoolr::install_exiftool()
```

<details>
<summary>Code</summary>

``` r
require(tidyverse)
require(exiftoolr)
require(terra)
require(rstudioapi)
```

</details>
