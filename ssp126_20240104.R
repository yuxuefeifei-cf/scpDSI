#########################################################
originwd<-"/media/fwq/8TB/cf/R"
setwd(originwd)

library(pdsi)
library(dplR)
library(ncdf4)

getwd()

prdata1<-nc_open("./ssp126/pr/ACCESS_CM2_pr_resampled.bc.nc")
prdata2<-nc_open("CanESM5.pr.bc.nc")
prdata3<-nc_open("CESM2_WACCM.pr.bc.nc")
prdata4<-nc_open("FGOALS_g3_time.pr.bc.nc")
prdata5<-nc_open("FIO_ESM_2_0.pr.bc.nc")
prdata6<-nc_open("MIROC6.pr.bc.nc")
prdata7<-nc_open("MRI_ESM2_0.pr.bc.nc")

tasdata1<-nc_open("./ssp126/tas/ACCESS_CM2_tas_resampled.bc.nc")
tasdata2<-nc_open("CanESM5.tas.bc.nc")
tasdata3<-nc_open("CESM2_WACCM.tas.bc.nc")
tasdata4<-nc_open("FGOALS_g3_time.tas.bc.nc")
tasdata5<-nc_open("FIO_ESM_2_0.tas.bc.nc")
tasdata6<-nc_open("MIROC6.tas.bc.nc")
tasdata7<-nc_open("MRI_ESM2_0.tas.bc.nc")

awcdata<-nc_open("awc_resampled.nc")
awc<-ncvar_get(awcdata,'AWC_CLASS')
lon <- ncvar_get(awcdata,'lon') 
lat <- ncvar_get(awcdata,'lat')
lon_lat <- expand.grid(lon=lon, lat=lat)
awc_vec <- as.vector(awc)
awc_df <- data.frame(lon_lat, awc=awc_vec)

pr1<-ncvar_get(prdata1,'pr')
pr2<-ncvar_get(prdata2,'pr')
pr3<-ncvar_get(prdata3,'pr')
pr4<-ncvar_get(prdata4,'pr')
pr5<-ncvar_get(prdata5,'pr')
pr6<-ncvar_get(prdata6,'pr')
pr7<-ncvar_get(prdata7,'pr')

tas1<-ncvar_get(tasdata1,'tas')
tas2<-ncvar_get(tasdata2,'tas')
tas3<-ncvar_get(tasdata3,'tas')
tas4<-ncvar_get(tasdata4,'tas')
tas5<-ncvar_get(tasdata5,'tas')
tas6<-ncvar_get(tasdata6,'tas')
tas7<-ncvar_get(tasdata7,'tas')

pr<-pr1
tas<-tas1

year<-rep(1901:2100,each=12)
month<-rep(1:12,time = 200)
data<-cbind(year,month)
dim(data)

lon <- ncvar_get(prdata1,'lon') 
lat <- ncvar_get(prdata1,'lat')
# lon_lat <- expand.grid(lon=lon, lat=lat)
# pr_data<-pr[,,1]
# pr_vec <- as.vector(pr)
# pr_df <- data.frame(lon_lat, pr=pr_vec)

lat<-as.data.frame(lat)
lon<-as.data.frame(lon)

#' Calculation of (sc)PDSI
library(bootRes)
# data(muc.clim)
library(pdsi)
library(dplR)

scpdsi<-pdsi(120, 50, muc.clim, 1960, 2000,mode = "scpdsi",verbose = TRUE)
#scpdsi2<-pdsi(125, 82.46, data_subset, 2016, 2100,mode = "scpdsi",verbose = TRUE)
# verbose <- TRUE

library(treeclim)
library(devtools)
library(digest)
# install_github("cszang/pdsi")
# install_github("https://github.com/cszang/pdsi.git")
# data_subset$tas_subset<-data_subset$tas_subset*1.8+32
# data_subset$pr_subset<-data_subset$pr_subset/25.4
start<-1901
end<-2100
df4<-na.omit(awc_df)
dir.create("./ssp126/temp")
setwd("./ssp126/temp")
exec_path <- file.path(system.file(package = "pdsi"), "scpdsi.o")
m <- array(NA, dim = c(240,120,2400))

# m2<-m
###########先运行pmat
pmat <- function(x, start = 1, end = 12, vnames = NULL) {
  years <- unique(x[, 1])
  n <- length(years)
  no.vars <- dim(x)[2] - 2 # number of variables
  months <- paste(c(rep("prev.", 12), rep("curr.", 12)), rep(c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), 2), sep = "")
  month.ids <- c(-1:-12, 1:12)
  used.months <- months[which(month.ids == start):which(month.ids == end)]
  no.months <- length(used.months)
  
  ## check for specified variable names, else default to colnames of
  ## x, if not present to V1, V2 etc.
  if (is.null(vnames) || length(unique(vnames)) != no.vars) {
    if (!is.null(colnames(x))) {
      vnames <- colnames(x)[3:(no.vars+2)]
    } else {
      vnames <- paste(rep("V", no.vars), 1:no.vars, sep = "")
    }
  }
  
  ## create unique names for variables
  vnames.mat <- matrix(NA, nrow = no.months, ncol = no.vars)
  for (i in 1:no.vars) {
    vnames.mat[, i] <- paste(vnames[i], ".", used.months, sep = "")
  }
  
  m <- matrix(NA, nrow = no.months*no.vars, ncol = n )
  colnames(m) <- years
  rownames(m) <- as.vector(vnames.mat)
  
  for (i in 1:n) {
    if (start < 0) {
      start.with <- which(x[, 1] == years[i - 1])[abs(start)] # start month in previous year
    } else {
      start.with <- which(x[, 1] == years[i])[start] # start month in current year
    }
    for (k in 1:no.vars) { # loop through variables
      for (j in 1:no.months) { # loop through months
        m[(j + (no.months*(k-1))), (i)] <- x[(start.with + j - 1), 2+k]
      }
    }
  }
  
  pmatrix <- as.data.frame(t(m))
  attributes(pmatrix)$npar <- no.vars
  
  pmatrix
}

for (i in 1:nrow(df4)) {
  index_x<- which(df4[i,1] == lon)
  index_y<- which(df4[i,2] == lat)
  latitude<-df4[i,2]
  awc<-df4[i,3]#/25.4
  pr_subset<-pr[index_x,index_y,]
  tas_subset<-tas[index_x,index_y,]
  data_subset<-as.data.frame(cbind(year,month,tas_subset,pr_subset))
  ## create temp dir
  tdir <- paste("./", index_x," ",index_y,sep = "")##################路径修改
  dir.create(tdir)
  climate<-data_subset
  ## truncate and reformat climate data
  climate_start <- which(climate[,1] == start)[1]
  climate_end <- which(climate[,1] == end)[12]
  climate <- climate[climate_start:climate_end,]
  climate_reform <- pmat(climate, start = 1, end = 12)
  climate_reform[is.na(climate_reform)]<-0
  ## split in temp and prec
  pmat_temp <- climate_reform[,1:12]
  pmat_prec <- climate_reform[,13:24]
  ## write to files
  temp_path <- file.path(tdir, "monthly_T")
  prec_path <- file.path(tdir, "monthly_P")
  write.table(pmat_temp, temp_path, col.names = F, quote = F)
  write.table(pmat_prec, prec_path, col.names = F, quote = F)
  ## calculate mean values and write to files
  normal_temp <- round(t(as.vector(colMeans(pmat_temp))), 3)
  normal_prec <- round(t(as.vector(colMeans(pmat_prec))), 3)
  normal_temp_path <- file.path(tdir, "mon_T_normal")
  normal_prec_path <- file.path(tdir, "mon_P_normal")
  write.table(normal_temp, normal_temp_path, col.names = F, quote = F,
              row.names = F)
  write.table(normal_prec, normal_prec_path, col.names = F, quote = F,
              row.names = F)
  ## write parameter files to tempdir
  params <- t(c(awc, latitude))
  param_path <- file.path(tdir, "parameter")
  write.table(params, param_path, col.names = F, quote = F,
              row.names = F)
  oldwd <- getwd()
  setwd(tdir)
  cmd <- paste(exec_path, " -m -i", shQuote(tdir), start, end)
  system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  setwd(oldwd)
  scpdsi_path <- file.path(tdir, "monthly", "self_cal", "PDSI.tbl")
  scPDSI <- read.fwf(scpdsi_path, c(5, rep(7, 12)))
  colnames(scPDSI) <- c("YEAR", toupper(month.abb))
  scPDSI<-scPDSI[,-1]
  scPDSI<- as.matrix(t(scPDSI))
  dim(scPDSI)<-c(200*12)
  m[index_x,index_y,] <- scPDSI
}

##################################################################################################
lon <- ncvar_get(awcdata,'lon') 
lat <- ncvar_get(awcdata,'lat')
time<-ncvar_get(prdata1,"time")##############修改pr

lon <- ncdim_def( "lon", "degrees_east", vals=lon)
lat <- ncdim_def( "lat", "degrees_north", vals=lat)
time <- ncdim_def( "time", "days since 1901-01-16", time)
setwd(originwd)
dir.create("./scpdsi")
setwd("./scpdsi")

# mv <- -9999 # missing value to use
scpdsi<- ncvar_def( name = 'scpdsi',units = 'scpdsi', dim = list(lon,lat,time), missval = NA, prec = 'float' )

#创建文档
# ncnew <- nc_create( filename = 'MRI_ESM2_0_scpdsi.nc', vars =list(scpdsi) )######修改模型名称
ncnew <- nc_create(filename = output_filename, vars = list(scpdsi))

#写入数据
ncvar_put( nc = ncnew, varid = scpdsi, vals = m)

#写入属性
nc_close(ncnew)






