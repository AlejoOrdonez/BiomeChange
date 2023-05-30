# Download CHALSA 2.1 [1980-2018]
rm(list=ls());gc()
require(snowfall)
setwd("/Users/alejandroordonez/Library/CloudStorage/Dropbox/Aarhus Assistant Professor/Projects/5. BiomeChange (BIOCHANGE)/BiomeChange//Data/CHELSA V21 Climatology/Data")
dir()
WDStrg <- getwd()
VARS <- dir()
sfInit( parallel=TRUE, cpus=5)
sfExport("WDStrg")
sfExport("VARS")

sfLapply (VARS,
          function(Vars){#(Vars <- VARS[2])
            setwd(paste0(WDStrg,"/",Vars))
            UrlList <- read.table("envidatS3paths.txt")
            for(DownLoadFile in UrlList$V1){#(DownLoadFile <- UrlList$V1[1])
              if(!strsplit(DownLoadFile,"CHELSA_")[[1]][2] %in% dir(paste0(WDStrg,"/",Vars))){
                download.file(url = DownLoadFile,
                              destfile = paste0(paste0(WDStrg,"/",Vars),"/",
                                                strsplit(DownLoadFile,"CHELSA_")[[1]][2]),
                              method = "wget",
                              quiet = TRUE)
              }		
            }
          })
sfStop()

