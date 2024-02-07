# Micasense RedEdge-MX DUAL processing
Simon Oiry

**WORK IN PROGRESS (last update : 2024-02-07 16:24:06.172526)**

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

I will personally use this workflow to **process single static images**,
where the usual structure-from-motion photogrammetry technique cannot be
used.

## Dual-MX Sensor <img src="img/Bandswv.png" width="50%" align="right" style="padding-left:10px;background-color:white;"/>

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
require(patchwork)
require(opencv) # install.packages("opencv", repos = "https://ropensci.r-universe.dev")
require(sf)
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

## Vignetting correction <img src="Output/plot/exemple_vignetting.png" width="50%" align="left" style="padding-left:10px;background-color:white;"/>

Vignetting refers to the reduction of image brightness toward the
periphery compared to the image center, a phenomenon due to the lens
properties of the camera. Before converting each image to radiance, the
digital numbers need to be corrected for this vignetting effect.

The function `vignette_map()` takes a `spatRaster` object or the path to
a TIFF file as input and outputs the vignetting map of this image :

The output of `vignette_map()` should give something looking the plot on
the left :

To correct the original image, we simply need to multiply the image by
the vignetting map. While it may not be immediately obvious, in the plot
below, the corners of the corrected images (right) are brighter than
those of the raw image (left).

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
    x <- vignette_poly_coef[6]*dist^6
      +vignette_poly_coef[5]*dist^5
      +vignette_poly_coef[4]*dist^4
      +vignette_poly_coef[3]*dist^3
      +vignette_poly_coef[2]*dist^2
      +vignette_poly_coef[1]*dist+1
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

# img_example<-meta$Image_Path[1]
# vignette_map_exemple<-vignette_map(img_example)
# 
# plot_exemple<-ggplot() +
#   geom_spatraster(data = vignette_map_exemple, aes(fill = lyr.1))+
#    scale_fill_viridis_c()+
#   labs(fill = "Correction factor")+
#   theme_bw()+
#   theme(legend.position = "top")
# 
# ggsave("export/plot/exemple_vignetting.png", plot_exemple, width = 10, height = 10)
```

</details>

<img src="Output/plot/micasense_compare.png" width="100%" align="left" style="padding-left:10px;background-color:white;"/>

The following code is a loop designed to correct all images present in
the Dual_MX_Images folder.

<details>
<summary>Code</summary>

``` r
img_list<-list.files("Dual_MX_Images", pattern = ".tif",recursive = T,full.names = T)

for(i in 1:length(img_list)){
  img_raw<-rast(img_list[i])
  img_map<-vignette_map(img_raw)
  img_corrected<-img_raw*img_map

  writeRaster(img_corrected, paste0("Output/RAW/Vignetting/Vign_",gsub(".*/","",img_list[i])),overwrite = T)
}
```

</details>

## Row Gradient correction <img src="Output/plot/exemple_RowGradient.png" width="50%" align="right" style="padding-left:10px;background-color:white;"/>

The next step involves correcting what MicaSense refers to as the Row
Gradient I haven’t been able to find any resources explaining what it
is, which led me to consult ChatGPT.

It appears there’s a disparity in the amount of light captured at the
top of the sensor versus what’s recorded at the bottom.

The row gradient correction applied by MicaSense on raw images before
processing is a calibration step aimed at compensating for any
non-uniformities and artifacts that may be present across the rows of
the sensor in the captured images.

Within the `XMP:RadiometricCalibration` tag of MicaSense images’
metadata, one can find all the necessary information to correct for this
row gradient effect.

<details>
<summary>Code</summary>

``` r
row_gradient_map<-function(img){
  
  if(typeof(img) == "S4"){
    img<-gsub(paste0(getwd(),"/"),"",sources(img))
  }
  
  exposure_time<-exif_read(img)$ExposureTime
  RadiometricCalibration <- as.numeric(exif_read(img)$RadiometricCalibration[[1]])
  
  x_dim<-exif_read(img)$ImageWidth 
  y_dim<-exif_read(img)$ImageHeight
  
  y = as.vector(matrix(rep(seq(1:y_dim),x_dim),ncol = x_dim))
  
  
  R<- 1 / (1 + RadiometricCalibration[2] * y / exposure_time - RadiometricCalibration[3] * y)
  
  R<-rast(matrix(R, nrow = y_dim, ncol = x_dim))
  
  return(R)
}

### Plot an example of vignetting map : 
# 
# img_example<-meta$Image_Path[1]
# gradient_map_exemple<-row_gradient_map(img_example)
# 
# plot_exemple<-ggplot() +
#   geom_spatraster(data = gradient_map_exemple, aes(fill = lyr.1))+
#    scale_fill_viridis_c()+
#   labs(fill = "Correction factor")+
#   theme_bw()+
#   theme(legend.position = "top")
# 
# ggsave("Output/plot/exemple_RowGradient.png", plot_exemple, width = 10, height = 10)
```

</details>

The following code is a loop designed to correct all images present in
the Dual_MX_Images folder.

<details>
<summary>Code</summary>

``` r
img_list<-list.files("Dual_MX_Images", pattern = ".tif",recursive = T,full.names = T)

for(i in 1:length(img_list)){
  img_raw<-rast(img_list[i])
  img_map<-row_gradient_map(img_raw)
  img_corrected<-img_raw*img_map

  writeRaster(img_corrected, paste0("Output/RAW/Row_Gradient/Row_Grad_",gsub(".*/","",img_list[i])),overwrite = T)
}
```

</details>

## Subtract the dark level and adjust for vignette and row gradient

At this stage, we will simultaneously apply vignetting, row gradient,
and dark level corrections. The purpose of dark level correction is to
mitigate the camera sensor’s inherent noise and ensure that the baseline
level of the image data is accurately established. This enhances the
image quality and accuracy for analysis.

<details>
<summary>Code</summary>

``` r
img_correction <- function(img){
  
  if(typeof(img) == "S4"){
    img<-gsub(paste0(getwd(),"/"),"",sources(img))
  }
  
  img_raw<-rast(img)
  
  img_vignette<-vignette_map(img_raw)
  img_rowGradient<-row_gradient_map(img_raw)
  
  DarkLevel <- mean(as.numeric(str_split(exif$BlackLevel," " )[[1]]))
  
  L = img_vignette*img_rowGradient*(img_raw - DarkLevel) 
  
  describe(img_raw)
  
  return(L)
}

img_list<-list.files("Dual_MX_Images", pattern = ".tif", recursive = T, full.names = T)

for (i in 1:length(img_list)) {
  
  corrected_img<-img_correction(img_list[i])
  
  writeRaster(corrected_img, paste0("Output/RAW/Radiometric_Calibration/All_Corr_",gsub(".*/","",img_list[i])),overwrite = T)

}
```

</details>

## DN to Radiance

After adjusting the digital numbers for sensor and lens uncertainties,
we can convert them to radiance to ensure that each image is expressed
in the same unit (W/m^2/nm/sr). To achieve this, it is necessary to
retrieve the exposure time and the gain applied to each image. It’s
important to note that the ISO value stored in the EXIF data must be
divided by 100. This is because the gain is represented in the
photographic parameter ISO, with a base ISO of 100. Dividing the ISO
value by 100 allows us to obtain a numeric gain.

Note also that during this conversion, it’s essential to normalize by
the image’s bit depth (2^16 for 16-bit images, 2^12 for 12-bit images)
because the calibration coefficients are designed to work with
normalized input values. This normalization ensures the coefficients are
applied correctly across different image bit depths. It’s important to
note that the `terra` package in R does not support the management of
TIFF file metadata, meaning there is no EXIF data in files saved using
`writeRaster()`. To address this issue, I am extracting metadata from
the raw files and associating it with the corrected files.

<details>
<summary>Code</summary>

``` r
DN_to_Radiance<-function(img){
  
  if(typeof(img) == "S4"){
    img_RAW<-img
    img<-gsub(paste0(getwd(),"/"),"",sources(img))
  }else{
    img_RAW<-rast(img)
  }
 
  img_name_current<-gsub(".*/","",img) %>% substr(.,10,24)
  img_path_RAW<-list.files("Dual_MX_Images",pattern = ".tif", recursive = T, full.names = T) %>% 
    as.data.frame() %>% 
    rename(path = ".") %>% 
    mutate(image_name = gsub(".*/","",path)) %>% 
    filter(image_name == img_name_current) %>% 
    pull(path)
    
   
  exposure_time<-exif_read(img_path_RAW)$ExposureTime
  Gain<-exif_read(img_path_RAW)$ISOSpeed/100
  bitsPerPixel<-exif_read(img_path_RAW)$BitsPerSample
  dnMax <- 2**bitsPerPixel
  a1<- as.numeric(exif_read(img_path_RAW)$RadiometricCalibration[[1]])[1]
  
  radianceImage <- img_RAW/(Gain*exposure_time)*a1/dnMax
  
  return(radianceImage)
}

img_path_cal<-list.files("Output/RAW/Radiometric_Calibration",pattern = ".tif", full.names = T)

for(i in 1 : length(img_path_cal)){
  
  Image_Radiance<-DN_to_Radiance(img_path_cal[i])
  
   writeRaster(Image_Radiance, paste0("Output/Radiance/Radiance_",gsub(".*/","",img_path_cal[i]) %>% gsub("All_Corr_","",.)),overwrite = T)
}
```

</details>

## Radiance to Reflectance <img src="Output/plot/exemple_QR_detection.png" width="50%" align="left" style="padding-left:10px;background-color:white;"/>

Now that we have a flat and calibrated radiance image, we can convert
into reflectance. To do this, we will use the radiance values of the
panel image of known reflectance to determine a scale factor between
radiance and reflectance.

Initially, it’s essential to detect the calibration panel in images. To
achieve this, we must locate the QR code of the calibration panel on
each image by utilizing the `opencv` library. The function
`ocv_qr_detect()` is used to find the coordinates of the QR code’s
corners as shown on this image.

Note that the Y-axis references used by `opencv` and `terra` differ,
necessitating a correction to ensure that points are plotted correctly
on the image. Please note also that this detection process is performed
on the raw image because the `OpenCV` library cannot open 32-bit data;
it is only capable of handling 16-bit images.

<details>
<summary>Code</summary>

``` r
Qr_detection<-function(img){
  
  if(typeof(img) == "S4"){
    img<-gsub(paste0(getwd(),"/"),"",sources(img))
  }
  
  img_openCV<-opencv::ocv_read(img) ### Open the image with opencv
  y_dim<-exif_read(img)$ImageHeight ### Retrieve the number of rows of the image
  
  qr_coordinate<-attr(opencv::ocv_qr_detect(img_openCV),which = "points") %>% 
  as.data.frame() %>% 
  rename(x = "V1",
         y = "V2") %>% 
  mutate(y= y_dim - y ,# Convert Y coordinate to match terra coordinate
         names = paste0("A_",c(1:4)))
  
  return(qr_coordinate)
}


df<-Qr_detection("Dual_MX_Images/Red/IMG_0002_1.tif")

img_rast<-rast("Dual_MX_Images/Red/IMG_0002_1.tif")
names(img_rast)<-"value"


### PLOT
plot_qr<-ggplot()+
  geom_spatraster(data = img_rast, aes(fill = value))+
  labs(fill = "DN")+
  scale_fill_gradientn(
    colours = grey(0:100 / 100),
    na.value = "transparent",
    trans = "sqrt"
  )+
  geom_point(data = df, aes(x =x , y = y), color = "red", size = 3)+
  geom_text(data = df, aes(x =x+30 , y = y+30, label = names), fontface = "bold",color = "red", size = 6)+
  theme_void()+
  theme(legend.position = "none")+
  coord_equal()

ggsave("Output/plot/exemple_QR_detection.png", plot_qr, width = 10, height = 10)
```

</details>
