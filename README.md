# Pesticide ecological risk screening tool

An R Shiny app that screens pesticide risk to waterways, pollinators and soil fauna at a chosen site, with human health information for context. Site soil, hydrology, rainfall and ground cover are fetched for any location; Australian users can narrow compounds by crop, pest and product label rate.

Version 1 (18 November 2026).

Authors: A. D. Canning, S. Lewis and D. J. Irvine

## Run

```r
install.packages(c("shiny", "dplyr", "ggplot2", "forcats", "ggpubr", "DT", "ggrepel", "jsonlite", "shinyjs",
                   "leaflet", "leaflet.extras", "httr", "jose", "tidyr", "gridExtra", "patchwork"))
shiny::runApp()
```

## Files

| File | Purpose |
|---|---|
|`app.R`|The app|
|`Pesticides_APVMA.csv`|Compound properties, toxicity and reference rates|
|`Pesticides_r.csv`|Human health classifications|
|`Pesticides_APVMA_uses.rds`|Registered crops and pests by product (PubCRIS)|
|`Pesticides_APVMA_use_rates.rds`|Label rates by product, crop and pest|

The four data files are not included in this repository because they contain values from the Pesticide Properties Database (PPDB), which cannot be redistributed without a licence from the University of Hertfordshire. Contact the authors for access. The underlying values can be viewed openly at their sources: compound properties and toxicity in the [PPDB](https://sitem.herts.ac.uk/aeru/ppdb/), and product registrations and labels in [APVMA PubCRIS](https://portal.apvma.gov.au/pubcris).

## Site data

Live site data need a Google Earth Engine service-account key. Save it as `gee-service-account.json` beside `app.R`, or point to it:

```r
Sys.setenv(PESTICIDE_GEE_KEY = "path/to/key.json")
```

Never commit the key. Without it the app still runs, using default site values that can be overwritten by hand.

## Data sources

PPDB, APVMA PubCRIS and product labels, ANZG default guideline values, HydroATLAS, CHIRPS, CSIRO SLGA, SoilGrids, OpenLandMap, MERIT Hydro and MERIT DEM, Landsat NDVI. Full references are listed in the app.

## Licence and citation

Code: [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Free to use, share and adapt for non-commercial purposes with credit. Commercial use needs written permission from the authors.

If you use the tool, please cite: Canning, A. D., Lewis, S. and Irvine, D. J. (2026). A globally applicable pesticide ecological risk screening tool with site-specific environmental context. Manuscript in preparation.

## Use

A screening tool for comparing compounds and practices at a site. It does not replace product labels or regulatory advice.
