rm(list = ls())

library(dplyr)
library(xts)
library(raster)
library(ED2.Congo)

################################################################################
# Extract

slat = 1.
slon = 24.5
in.path = "/data/gent/vo/000/gvo00074/ED_common_data/met/CB/"
start_date = "1960-01-01"
end_date = "2019-12-31"
outfolder = "/data/gent/vo/000/gvo00074/ED_common_data/met/CB/extracted"
in.prefix = "ERA5_"
newsite = "YGB"
vars = NULL
overwrite = TRUE


years <- seq(lubridate::year(start_date),
             lubridate::year(end_date),
             1
)
ensemblesN <- seq(1, 1)


tryCatch({
  #for each ensemble
  one.year.out <- years %>%
    purrr::map(function(year) {

      # for each year
      point.data <-  ensemblesN %>%
        purrr::map(function(ens) {


          ncfile <- file.path(in.path, paste0(in.prefix, year, ".nc"))

          PEcAn.logger::logger.info(paste0("Trying to open :", ncfile, " "))

          if (!file.exists(ncfile)){PEcAn.logger::logger.severe("The nc file was not found.")}

          #msg
          PEcAn.logger::logger.info(paste0(year, " is being processed ", "for ensemble #", ens, " "))
          #open the file
          nc_data <- ncdf4::nc_open(ncfile)
          # time stamp

          t <- ncdf4::ncvar_get(nc_data, "time")
          tunits <- ncdf4::ncatt_get(nc_data, 'time')
          tustr <- strsplit(tunits$units, " ")
          timestamp <-
            as.POSIXct(t * 3600, tz = "UTC", origin = tustr[[1]][3])
          try(ncdf4::nc_close(nc_data))


          # set the vars
          if (is.null(vars))
            vars <- names(nc_data$var)
          # for the variables extract the data

          all.data.point <- vars %>%
            purrr::map_dfc(function(vname) {
              PEcAn.logger::logger.info(paste0(" \t ",vname, "is being extracted ! "))

              brick.tmp <-
                raster::brick(ncfile, varname = vname, level = ens)
              nn <-
                raster::extract(brick.tmp,
                                sp::SpatialPoints(cbind(slon, slat)),
                                method = 'simple')

              if (!is.numeric(nn)) {
                PEcAn.logger::logger.severe(paste0(
                  "Expected raster object to be numeric, but it has type `",
                  paste0(typeof(nn), collapse = " "),
                  "`"
                ))
              }


              # replacing the missing/filled values with NA
              nn[nn == nc_data$var[[vname]]$missval] <- NA
              # send out the extracted var as a new col
              t(nn)

            }) %>%
            `colnames<-`(vars)
          #close the connection

          # send out as xts object
          xts::xts(all.data.point, order.by = timestamp)
        }) %>%
        setNames(paste0("ERA_ensemble_", ensemblesN))

      #Merge mean and the speard
      return(point.data)

    }) %>%
    setNames(years)


  # The order of one.year.out is year and then Ens - Mainly because of the spead  / I wanted to touch each file just once.
  # This now changes the order to ens - year
  point.data <- ensemblesN %>%
    purrr::map(function(Ensn) {
      one.year.out %>%
        purrr::map( ~ .x [[Ensn]]) %>%
        do.call("rbind.xts", .)
    })


  # Calling the met2CF inside extract bc in met process met2CF comes before extract !
  out <- met2CF.ERA5(
    slat,
    slon,
    start_date,
    end_date,
    sitename=newsite,
    outfolder,
    point.data,
    overwrite = TRUE,
    verbose = TRUE
  )

}, error = function(e) {
  PEcAn.logger::logger.severe(paste0(conditionMessage(e)))
})

saveRDS(point.data,"TS_YGB.RDS")

PEcAn.ED2::met2model.ED2(in.path =  file.path(outfolder,"ERA5_YGB_1"),
                         in.prefix = "ERA5.1",
                         outfolder = file.path(outfolder, "ERA5_YGB_1","ED2"),
                         start_date = start_date,
                         end_date = end_date,
                         lat = slat,
                         lon = slon,
                         overwrite = TRUE)

# scp /home/femeunier/Documents/projects/YGB/scripts/process.ERA5.files.R hpc:/data/gent/vo/000/gvo00074/felicien/R
