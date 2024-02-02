# Micasense RedEdge-MX DUAL processing
Simon Oiry

This workflow adapts the Micasense workflow for manual processing of
images from the Micasense RedEdge-MX Dual camera. he original workflow,
written in Python, is available
[here](https://github.com/micasense/imageprocessing). This repository
aims to translate the Python workflow into an R workflow. The python
workflow is no longer maintained and requires specific versions of each
package to function correctly. Don’t enjoy spending hours just setting
up your Python environment? This repo is made for you!

The original aims of micasense when they created this processing
workflow was to help researchers and developers to do their own image
processing. While a number of commercial tools fully support processing
MicaSense data into reflectance maps, there are a number of reasons to
process your own data, including controlling the entire radiometric
workflow (for academic or publication reasons), pre-processing images to
be used in a non-radiometric photogrammetry suite, or processing single
sets of 10 images without building a larger map.

I will personally use this workflow to process single static images,
where the usual structure-from-motion photogrammetry technique cannot be
used.

## Dual-MX Sensor <img src="img/Bandswv.png" width="50%" align="right" style="padding-left:10px;background-color:white;" />

The dual-MX camera have a spectral resolution of 10 bands, ranging from
the blue (444nm) to the NIR (840nm).

## Packages

The first thing to do is to ensure that all the packages are ready to be
used. The [exiftoolr](https://github.com/JoshOBrien/exiftoolr) packages
is used to read Exif of tiff files. After installing the package you
will need to run this line of code : `exiftoolr::install_exiftool()`.

<details>
<summary>Code</summary>

``` r
require(tidyverse)
require(exiftoolr)
require(terra)
```

</details>

## Locate images and Reading metadata used to orthorectify, calibrate and align images

The code is used to find the path for each individual image, identify
each band, and extract all the necessary metadata.

<details>
<summary>Code</summary>

``` r
image_df<-"Dual_MX_Images" %>% 
  list.files(recursive = T, full.names = T) %>% 
  as.data.frame() %>% 
  rename(path = ".") %>% 
  mutate(image_name = gsub(".*/","",path),
         image_ID = substr(image_name,5,8),
         Band = paste0("B",gsub(".*_","",image_name) %>% gsub(".tif","",.))) 

meta <-data.frame(
  Image_name = image_df$image_name,
  Image_Path = image_df$path,
  Make = NA,
  Model = NA,
  Exposure_Time = NA,
  Gain = NA,
  Resolution = NA,
  Band_Name = NA,
  Central_Wavelength = NA,
  Band_Width = NA,
  Capture_ID = NA,
  Flight_ID = NA,
  Focal_Length = NA,
  Black_Level = NA,
  Radiometric_Calibration_a1 = NA,
  Radiometric_Calibration_a2 = NA,
  Radiometric_Calibration_a3 = NA
)

for (i in 1:nrow(image_df)){
  exif<-exif_read(image_df$path[i])
  
  meta$Make[i]<-exif$Make
  meta$Model[i]<-exif$Model
  meta$Exposure_Time[i]<-exif$ExposureTime
  meta$Gain[i]<-exif$ISOSpeed
  meta$Resolution[i]<-paste0(exif$ImageWidth,"x",exif$ImageHeight)
  meta$Band_Name[i]<-exif$BandName
  meta$Central_Wavelength[i]<-exif$CentralWavelength
  meta$Band_Width[i]<-exif$WavelengthFWHM
  meta$Capture_ID[i]<-exif$CaptureId
  meta$Flight_ID[i]<-exif$FlightId
  meta$Focal_Length[i]<-exif$FocalLength
  meta$Black_Level[i]<-mean(as.numeric(str_split(exif$BlackLevel," " )[[1]]))
  meta$Radiometric_Calibration_a1[i]<-as.numeric(exif$RadiometricCalibration[[1]][1])
  meta$Radiometric_Calibration_a2[i]<-as.numeric(exif$RadiometricCalibration[[1]][2])
  meta$Radiometric_Calibration_a3[i]<-as.numeric(exif$RadiometricCalibration[[1]][3])
}
```

</details>

|    Name of metadata     |                      Description                       |    Unit    |                                                                                                                                                 Comments                                                                                                                                                 |
|:-----------------------:|:------------------------------------------------------:|:----------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|      Exposure_Time      |              Exposure time of the picture              |   Second   |                                                                                                                                                                                                                                                                                                          |
|          Gain           |               Gain applied to the image                |  Unitless  |                                                                                                                                                                                                                                                                                                          |
|       Resolution        |                Resolution of the image                 |   pixels   |                                                                                                                                                                                                                                                                                                          |
|   Central_Wavelength    |             Central wavelength of the band             | Nanometer  |                                                                                                                                                                                                                                                                                                          |
|       Band_Width        |       Full Width Half Maximum (FWHM) of the band       | Nanometer  |                                                                                                                                                                                                                                                                                                          |
|       Capture_ID        |         Unique ID of all the band of the image         |            |                                                                                                                                                                                                                                                                                                          |
|        Flight_ID        |        Unique ID of all the image of the flight        |            |                                                                                                                                                                                                                                                                                                          |
|      Focal_Length       |           Focal length of the optical system           | millimeter |                                                                                                                                                                                                                                                                                                          |
|       Black_Level       |                    darkPixel value                     |     DN     | Average of 4 values. These values come from optically-covered pixels on the imager which are exposed at the same time as the image pixels. They measure the small amount of random charge generation in each pixel, independent of incoming light, which is common to all semiconductor imaging devices. |
| Radiometric_Calibration | Optical parameter used for the row gradient correction |            |                                                                                                                                       imager-specific calibrations                                                                                                                                       |

Description of the metadata extracted from micasense images

## Vignetting correction
