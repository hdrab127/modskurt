#' Tuangi (*Austrovenus structburyi*, NZ Cockle) estuarine abundance
#'
#' Discrete counts of tuangi in intertidal estuaries of New Zealand
#'
#' @format ## `tuangi`
#' A data frame with 764 rows and 11 columns:
#' \describe{
#'   \item{year}{Year of survey}
#'   \item{lon, lat}{Longitude and latitude of the sample}
#'   \item{count}{Total tuangi counted in a 250 gram sediment core}
#'   \item{mud_pct}{\% of sediment that was muddy (<63microns?) at the time of sampling}
#'   \item{tn_conc}{total nitrogen concentration (mg/L) at the time of sampling}
#'   \item{toc_pct}{total organic carbon content (\%) at the time of sampling}
#'   \item{tp_conc}{total phosphorus concentration (mg/L) at the time of sampling}
#'   \item{sst_min, sst_avg, sst_max}{Min, mean, and max sea surface temperature (deg C) at the location of sampling over the 12 months preceding}
#' }
#' @source <https://www.saltecology.co.nz/>
#' @usage data(tuangi)
"tuangi"

#' *Aonides trifida* estuarine abundance
#'
#' Discrete counts of aonides trifida in intertidal estuaries of New Zealand
#'
#' @format ## `aonides`
#' A data frame with 251 rows and 7 columns:
#' \describe{
#'   \item{estu}{Estuary where the count was observed}
#'   \item{year}{Year of survey}
#'   \item{count}{Total counted in a 250 gram sediment core}
#'   \item{mud_pct}{\% of sediment that was muddy (<63microns?) at the time of sampling}
#'   \item{tn_conc}{total nitrogen concentration (mg/L) at the time of sampling}
#'   \item{toc_pct}{total organic carbon content (\%) at the time of sampling}
#'   \item{tp_conc}{total phosphorus concentration (mg/L) at the time of sampling}
#' }
#' @source <https://www.saltecology.co.nz/>
#' @usage data(aonides)
"aonides"
