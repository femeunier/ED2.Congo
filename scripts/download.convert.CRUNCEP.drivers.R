rm(list = ls())

# Workflow ED2
# Libraries

library(ncdf4)
library(stringr)
library(dplyr)
library(udunits2)

################################################################################

census.file <- "./data/2_data_by_censusinterval.csv"
data.census.interval <- read.csv(census.file)
sites <- data.census.interval$PlotCode
sites <- sites[grepl("YGB",sites)]

################################################################################
method = "ncss"
maxErrors = 10
sleep = 2
verbose = FALSE
overwrite = TRUE

###############################################################################
# CO2
fileCO2 <- "./data/CO2_1700_2019_TRENDYv2020.txt"
dataC02 <- read.table(fileCO2,stringsAsFactors = FALSE)

dataCO2.n <- dataC02 %>% mutate(years = str_sub(V1,7,10),
                                CO2 = as.numeric(str_sub(V1,12,17))) %>% dplyr::select(years,CO2)



for (isite in seq(1,length(sites))){

  csite <- sites[isite]
  print(csite)

  cdf  <- data.census.interval %>%
    filter(PlotCode == csite)

  # Years of drivers
  years <- min(floor(cdf$int_ini)):max(ceiling(cdf$int_fin))

  # Your site coordinates
  lon <- mean(cdf$Lon)
  lat <- mean(cdf$Lat)

  # Output folders
  outfolder <- file.path("./inputs/drivers",csite)
  dir.create(outfolder, showWarnings = FALSE)

  download.CRUNCEP(outfolder,
                   start_date = paste0(min(years),"-01-01"),
                   end_date = paste0(max(years),"-12-31"), lat, lon,
                   overwrite = overwrite, verbose = verbose)

  # Step 2) Convert drivers into a ED2-readable format

  directory <- outfolder

  in.path = directory
  in.prefix = "CRUNCEP"
  outfolder = file.path(directory,"ED2")
  start_date = paste0(min(years),"/01/01")
  end_date = paste0(max(years),"/12/31")
  lst = 0
  lat = NA
  lon = NA
  overwrite = FALSE
  verbose = FALSE
  leap_year = TRUE

  ###############################################################################

  overwrite <- as.logical(overwrite)
  start_date <- as.POSIXlt(start_date, tz = "UTC")
  end_date <- as.POSIXlt(end_date, tz = "UTC")
  met_folder <- outfolder
  met_header_file <- file.path(met_folder, "ED_MET_DRIVER_HEADER")
  results <- data.frame(file = met_header_file, host = PEcAn.remote::fqdn(),
                        mimetype = "text/plain", formatname = "ed.met_driver_header files format",
                        startdate = start_date, enddate = end_date, dbfile.name = "ED_MET_DRIVER_HEADER",
                        stringsAsFactors = FALSE)
  dir.create(met_folder, recursive = TRUE, showWarnings = FALSE)
  dm <- c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305,
          335, 366)
  dl <- c(0, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306,
          336, 367)
  month <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
             "AUG", "SEP", "OCT", "NOV", "DEC")
  mon_num <- c("01", "02", "03", "04", "05", "06", "07", "08",
               "09", "10", "11", "12")
  day2mo <- function(year, day, leap_year) {
    mo <- rep(NA, length(day))
    if (!leap_year) {
      mo <- findInterval(day, dm)
      return(mo)
    }
    else {
      leap <- lubridate::leap_year(year)
      mo[leap] <- findInterval(day[leap], dl)
      mo[!leap] <- findInterval(day[!leap], dm)
      return(mo)
    }
  }
  start_year <- lubridate::year(start_date)
  end_year <- lubridate::year(end_date)
  year_seq <- seq(start_year, end_year)
  day_secs <- udunits2::ud.convert(1, "day", "seconds")
  need_input_files <- file.path(in.path, paste(in.prefix,
                                               year_seq, "nc", sep = "."))
  have_input_files <- file.exists(need_input_files)
  if (!all(have_input_files)) {
    PEcAn.logger::logger.severe("Missing the following required input files: ",
                                paste(sprintf("'%s'", need_input_files[!have_input_files]),
                                      collapse = ", "))
  }
  month_seq <- seq(lubridate::floor_date(start_date, "month"),
                   lubridate::floor_date(end_date, "month"), by = "1 month")
  target_fnames <- paste0(toupper(strftime(month_seq, "%Y%b",
                                           tz = "UTC")), ".h5")
  target_out_files <- file.path(met_folder, target_fnames)
  have_target_out_files <- file.exists(target_out_files)
  if (any(have_target_out_files)) {
    if (overwrite) {
      PEcAn.logger::logger.warn("The following existing target output files will be overwritten:",
                                paste(sprintf("'%s'", target_out_files[have_target_out_files]),
                                      collapse = ", "))
    }
    else {
      have_output_byyear <- split(have_target_out_files,
                                  lubridate::year(month_seq))
      complete_years <- vapply(have_output_byyear, all,
                               logical(1))
      skip_years <- tryCatch(as.numeric(names(complete_years[complete_years])),
                             warning = function(e) PEcAn.logger::logger.severe(e))
      PEcAn.logger::logger.warn("The following output files already exist:",
                                paste(target_out_files[have_target_out_files]),
                                ". This means the following complete years will be skipped: ",
                                skip_years)
      year_seq <- setdiff(year_seq, skip_years)
    }
  }
  for (year in year_seq) {
    ncfile <- file.path(in.path, paste(in.prefix, year,
                                       "nc", sep = "."))
    nc <- ncdf4::nc_open(ncfile)
    flat <- try(ncdf4::ncvar_get(nc, "latitude"), silent = TRUE)
    if (!is.numeric(flat)) {
      flat <- nc$dim[[1]]$vals[1]
    }
    if (is.na(lat)) {
      lat <- flat
    }
    else if (lat != flat) {
      PEcAn.logger::logger.warn("Latitude does not match that of file",
                                lat, "!=", flat)
    }
    flon <- try(ncdf4::ncvar_get(nc, "longitude"), silent = TRUE)
    if (!is.numeric(flon)) {
      flat <- nc$dim[[2]]$vals[1]
    }
    if (is.na(lon)) {
      lon <- flon
    }
    else if (lon != flon) {
      PEcAn.logger::logger.warn("Longitude does not match that of file",
                                lon, "!=", flon)
    }

    # lat <- eval(parse(text = lat))
    # lon <- eval(parse(text = lon))

    sec <- nc$dim$time$vals
    Tair <- ncdf4::ncvar_get(nc, "air_temperature")
    Qair <- ncdf4::ncvar_get(nc, "specific_humidity")
    U <- try(ncdf4::ncvar_get(nc, "eastward_wind"), silent = TRUE)
    V <- try(ncdf4::ncvar_get(nc, "northward_wind"), silent = TRUE)
    Rain <- ncdf4::ncvar_get(nc, "precipitation_flux")
    pres <- ncdf4::ncvar_get(nc, "air_pressure")
    SW <- ncdf4::ncvar_get(nc, "surface_downwelling_shortwave_flux_in_air")
    LW <- ncdf4::ncvar_get(nc, "surface_downwelling_longwave_flux_in_air")

    CO2.actual_level <- dataCO2.n %>% filter(years == year) %>% pull(CO2)
    CO2 <- CO2.actual_level*(LW**0)


    use_UV <- is.numeric(U) & is.numeric(V)
    if (!use_UV) {
      U <- try(ncdf4::ncvar_get(nc, "wind_speed"), silent = TRUE)
      if (is.numeric(U)) {
        PEcAn.logger::logger.info("eastward_wind and northward_wind are absent, using wind_speed to approximate eastward_wind")
        V <- rep(0, length(U))
      }
      else {
        PEcAn.logger::logger.severe("No eastward_wind and northward_wind or wind_speed in the met data")
      }
    }
    useCO2 <- is.numeric(CO2)
    sec <- udunits2::ud.convert(sec, unlist(strsplit(nc$dim$time$units,
                                                     " "))[1], "seconds")
    ncdf4::nc_close(nc)
    dt <- PEcAn.utils::seconds_in_year(year, leap_year)/length(sec)
    toff <- -as.numeric(lst) * 3600/dt
    # slen <- seq_along(SW)
    # Tair <- c(rep(Tair[1], toff), Tair)[slen]
    # Qair <- c(rep(Qair[1], toff), Qair)[slen]
    # U <- c(rep(U[1], toff), U)[slen]
    # V <- c(rep(V[1], toff), V)[slen]
    # Rain <- c(rep(Rain[1], toff), Rain)[slen]
    # pres <- c(rep(pres[1], toff), pres)[slen]
    # SW <- c(rep(SW[1], toff), SW)[slen]
    # LW <- c(rep(LW[1], toff), LW)[slen]
    # if (useCO2) {
    #   CO2 <- c(rep(CO2[1], toff), CO2)[slen]
    # }
    skip <- FALSE
    nyr <- floor(length(sec) * dt/86400/365)
    yr <- NULL
    doy <- NULL
    hr <- NULL
    asec <- sec
    for (y in seq(year, year + nyr - 1)) {
      diy <- PEcAn.utils::days_in_year(y, leap_year)
      ytmp <- rep(y, udunits2::ud.convert(diy/dt, "days",
                                          "seconds"))
      dtmp <- rep(seq_len(diy), each = day_secs/dt)
      if (is.null(yr)) {
        yr <- ytmp
        doy <- dtmp
        hr <- rep(NA, length(dtmp))
      }
      else {
        yr <- c(yr, ytmp)
        doy <- c(doy, dtmp)
        hr <- c(hr, rep(NA, length(dtmp)))
      }
      rng <- length(doy) - length(ytmp):1 + 1
      if (!all(rng >= 0)) {
        skip <- TRUE
        PEcAn.logger::logger.warn(year, " is not a complete year and will not be included")
        break
      }
      asec[rng] <- asec[rng] - asec[rng[1]]
      hr[rng] <- (asec[rng] - (dtmp - 1) * day_secs)/day_secs *
        24
    }
    mo <- day2mo(yr, doy, leap_year)
    if (length(yr) < length(sec)) {
      rng <- (length(yr) + 1):length(sec)
      if (!all(rng >= 0)) {
        skip <- TRUE
        PEcAn.logger::logger.warn(paste(year, "is not a complete year and will not be included"))
        break
      }
      yr[rng] <- rep(y + 1, length(rng))
      doy[rng] <- rep(1:366, each = day_secs/dt)[1:length(rng)]
      hr[rng] <- rep(seq(0, length = day_secs/dt, by = dt/day_secs *
                           24), 366)[1:length(rng)]
    }
    if (skip) {
      print("Skipping to next year")
      next
    }


    ugrdA <- U
    vgrdA <- V
    shA <- Qair
    tmpA <- Tair
    dlwrfA <- LW
    presA <- pres
    prateA <- Rain

    if (useCO2) {
      co2A <- CO2
    }

    nbdsfA <- nddsfA <- vbdsfA <- vddsfA <- NA*pres

    clat <- lat
    clon <- lon

    cosz <- PEcAn.data.atmosphere::cos_solar_zenith_angle(doy,
                                                          clat, clon, dt, hr)
    rpot <- 1366 * cosz

    if(all(is.na(SW))) next()
    SW[rpot < SW] <- rpot[rpot < SW]

    frac <- SW/rpot
    frac[frac > 0.9] <- 0.9
    frac[frac < 0] <- 0
    frac[is.na(frac)] <- 0
    frac[is.nan(frac)] <- 0

    SWd <- SW * (1 - frac)
    nbdsfA <- (SW - SWd) * 0.57
    nddsfA <- SWd* 0.48
    vbdsfA <- (SW - SWd) * 0.43
    vddsfA <- SWd * 0.52

    hgtA <- 50*(Rain**0)


    for (y in year + 1:nyr - 1) {
      sely <- which(yr == y)
      for (m in unique(mo[sely])) {
        selm <- sely[which(mo[sely] == m)]
        mout <- paste(met_folder, "/", y, month[m],
                      ".h5", sep = "")
        if (file.exists(mout)) {
          if (overwrite) {
            file.remove(mout)
            ed_met_h5 <- hdf5r::H5File$new(mout)
          }
          else {
            PEcAn.logger::logger.warn("The file already exists! Moving to next month!")
            next
          }
        }
        else {
          ed_met_h5 <- hdf5r::H5File$new(mout)
        }


        nbdsf <- nddsf <- vddsf <- vbdsf <- prate <-
          dlwrf <- pres <- hgt <- ugrd <- vgrd <- sh <- tmp <- co2 <- array(data = NA,dim = c(length(selm),1,1))

        nbdsf[,1,1] <- nbdsfA[selm]
        nddsf[,1,1] <- nddsfA[selm]
        vbdsf[,1,1] <- vbdsfA[selm]
        vddsf[,1,1] <- vddsfA[selm]
        prate[,1,1] <- prateA[selm]
        dlwrf[,1,1] <- dlwrfA[selm]
        pres[,1,1] <- presA[selm]
        hgt[,1,1] <- hgtA[selm]
        ugrd[,1,1] <- ugrdA[selm]
        vgrd[,1,1]<- vgrdA[selm]
        sh[,1,1]<- shA[selm]
        tmp[,1,1]<- tmpA[selm]
        if (useCO2) {
          co2[,1,1]<- co2A[selm]
        }

        Grid.lat <- lat
        Grid.lon <- lon

        ed_met_h5[["lat"]] <- Grid.lat
        ed_met_h5[["lon"]] <- Grid.lon
        ed_met_h5[["nbdsf"]] <- nbdsf
        ed_met_h5[["nddsf"]] <- nddsf
        ed_met_h5[["vbdsf"]] <- vbdsf
        ed_met_h5[["vddsf"]] <- vddsf
        ed_met_h5[["prate"]] <- prate
        ed_met_h5[["dlwrf"]] <- dlwrf
        ed_met_h5[["pres"]] <- pres
        ed_met_h5[["hgt"]] <- hgt
        ed_met_h5[["ugrd"]] <- ugrd
        ed_met_h5[["vgrd"]] <- vgrd
        ed_met_h5[["sh"]] <- sh
        ed_met_h5[["tmp"]] <- tmp
        if (useCO2) {
          ed_met_h5[["co2"]] <- co2
        }
        ed_met_h5$close_all()
      }
    }
    metvar <- c("nbdsf", "nddsf", "vbdsf", "vddsf", "prate",
                "dlwrf", "pres", "hgt", "ugrd", "vgrd", "sh", "tmp",
                "co2")
    metvar_table <- data.frame(variable = metvar, update_frequency = dt,
                               flag = 1)
    if (!useCO2) {
      metvar_table_vars <- metvar_table[metvar_table$variable !=
                                          "co2", ]
    }
    else {
      metvar_table_vars <- metvar_table
    }
    ed_metheader <- list(list(path_prefix = met_folder,
                              nlon = length(lon), nlat = length(lat), dx = 0.5, dy = 0.5, xmin = min(lon),
                              ymin = min(lat), variables = rbind(metvar_table_vars)))

    # ed_metheader <- list(list(path_prefix = met_folder,
    #                           nlon = length(lon), nlat = length(lat), dx = 0.5, dy = 0.5, xmin = min(lon),
    #                           ymin = min(lat), variables = metvar_table_vars))

    check_ed_metheader(ed_metheader)
    write_ed_metheader(ed_metheader, met_header_file, header_line = shQuote("Made_by_PEcAn_met2model.ED2"))
  }
  PEcAn.logger::logger.info("Done with met2model.ED2")

  # Step 3) run the ED2 model (cluster)

  # Modify ED_MED_DRIVER!

  system2("rsync",
          paste("-avz",
                file.path("/home/femeunier/Documents/projects/ED2.Congo/inputs/drivers/",csite),
                "hpc:/data/gent/vo/000/gvo00074/ED_common_data/met/"))


}

