###### custom function to add filename when reading data from excel spreadsheets

required <- list("tidyverse", 
                 "readxl", 
                 "colorspace",
                 "ggh4x",
                "glmmTMB",
                "emmeans",
                "ggpubr",
                "ggtext",
                "patchwork",
                "DHARMa",
                "scales",
                "sjPlot",
                "showtext")

lapply(required, library, character.only = T)


## Functions

read_plus <- function(flnm, sheet, skip) {
  read_excel(flnm ,sheet=sheet, skip=skip, na=c("na","NA", "-")) %>% 
    mutate(filename = flnm)
  
}


read_data <- function(sheet){
  list.files(path="./data/", pattern = "*.xlsx", full.names = T) %>% 
    map_df(~read_plus(., sheet=sheet, skip=1)) %>%  # read and input all files ending .xlsx
    mutate(filename=str_match(filename, "/data//Crossing data summary A_s*(.*?)\\s*B_1590")[,2])
}

#####

##### force bind, bind rows of data when column names do not match #####

force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  bind_rows(df1, df2)
}

#####


##### DHARMa_check function to simulate residuals and plot - for mixed models

DHARMa_check <- function(model){
  sim <- DHARMa::simulateResiduals(model) 
 plot(sim, asFactor=T)
}


