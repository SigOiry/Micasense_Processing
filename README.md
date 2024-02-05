# Micasense RedEdge-MX DUAL processing
Simon Oiry

This workflow adapts the Micasense workflow for manual processing of
images from the Micasense RedEdge-MX Dual camera. he original workflow,
written in Python, is available
[here](https://github.com/micasense/imageprocessing). This repository
aims to translate the Python workflow into an R workflow. I’m not really
used to code in Python and I thought it can be a good exercise to try to
translate this repository in R.

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
require(tidyterra)
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
  Unique_ID = NA,
  Image_Path = image_df$path,
  Make = NA,
  Model = NA,
  Exposure_Time = NA,
  Gain = NA,
  Width = NA,
  Height = NA,
  Band_Name = NA,
  Central_Wavelength = NA,
  Band_Width = NA,
  Capture_ID = NA,
  Flight_ID = NA,
  Focal_Length = NA,
  Black_Level = NA,
  Radiometric_Calibration_a1 = NA,
  Radiometric_Calibration_a2 = NA,
  Radiometric_Calibration_a3 = NA,
  Vignetting_Center_X = NA,
  Vignetting_Center_Y = NA,
  Vignetting_Polynomial = NA
)

for (i in 1:nrow(image_df)){
  exif<-exif_read(image_df$path[i])
  
  meta$Make[i]<-exif$Make
  meta$Model[i]<-exif$Model
  meta$Exposure_Time[i]<-exif$ExposureTime
  meta$Gain[i]<-exif$ISOSpeed
  meta$Width[i]<-exif$ImageWidth
  meta$Height[i]<-exif$ImageHeight
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
  meta$Vignetting_Center_Y[i]<-as.numeric(exif$VignettingCenter[[1]][1])
  meta$Vignetting_Center_y[i]<-as.numeric(exif$VignettingCenter[[1]][2])
  meta$Vignetting_Polynomial[i]<-c(exif$VignettingPolynomial)
  meta$Unique_ID[i]<-paste(sep = "_",meta$Capture_ID[i],meta$Central_Wavelength[i])
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

Vignetting refers to the reduction of image brightness toward the
periphery compared to the image center, a phenomenon due to the lens
properties of the camera. Before converting each image to radiance, the
digital numbers need to be corrected for this vignetting effect.

The function `vignette_map()` takes a `spatRaster` object or the path to
a TIFF file as input and outputs the vignetting map of this image :

<details>
<summary>Code</summary>

``` r
vignette_map<- function(img){
  
  if(typeof(img) == "S4"){
    img<-gsub(paste0(getwd(),"/"),"",sources(img))
  }
  
  metadata<-exif_read(img)
  x_dim <-metadata$ImageWidth
  y_dim <-metadata$ImageHeight
  # get vignette center
  x_vignette <- metadata$VignettingCenter[[1]][1]
  y_vignette <- metadata$VignettingCenter[[1]][2]

  # get vignette polynomial
  nvignette_poly <- length(metadata$VignettingPolynomial[[1]])
  vignette_poly_coef <- ((metadata$VignettingPolynomial[[1]]))
  
  vignette_poly<-function(dist,vignette_poly_coef = vignette_poly_coef){
    x <- vignette_poly_coef[6]*dist^6+vignette_poly_coef[5]*dist^5+vignette_poly_coef[4]*dist^4+vignette_poly_coef[3]*dist^3+vignette_poly_coef[2]*dist^2+vignette_poly_coef[1]*dist+1
    return(x)
  }
  
  
  # perform vignette correction
  # get coordinate grid across image
  
  x = matrix(rep(seq(1:x_dim),y_dim), nrow = y_dim, byrow = TRUE)
  y = matrix(rep(seq(1:y_dim),x_dim),ncol = x_dim)
  
  
  # compute matrix of distances from image center
  
  dist_list<-sqrt(((as.vector(x)-x_vignette)^2)+((as.vector(y)-y_vignette)^2))
  dist<-matrix(dist_list, ncol = x_dim,nrow = y_dim)
  
  
  vignette_list <-1/vignette_poly(dist,vignette_poly_coef)
  
  vignette<-rast(matrix(vignette_list, ncol = x_dim, nrow = y_dim))

  return(vignette)
}

### Plot an example of vignetting map : 

img_example<-meta$Image_Path[1]
vignette_map_exemple<-vignette_map(img_example)

plot_exemple<-ggplot() +
  geom_spatraster(data = vignette_map_exemple, aes(fill = lyr.1))+
   scale_fill_viridis_c()+
  labs(fill = "Correction factor")+
  theme_bw()+
  theme(legend.position = "top")

ggsave("export/plot/exemple_vignetting.png", plot_exemple, width = 10, height = 10)
```

</details>

The output of `vignette_map()` should give something looking like that :

<img src="export/plot/exemple_vignetting.png" width="50%" style="padding-left:25%;background-color:white;" />

To correct the original image, we simply need to multiply the image by
the vignetting map.
