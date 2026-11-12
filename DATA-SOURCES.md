# Data sources

The current application uses the Jura topsoil data distributed with `gstat`:
259 original training sites and 100 reserved validation sites. The response is
natural-log cadmium concentration; original cadmium is measured in mg/kg.
Geographic coordinates are transformed to EPSG:2056, with the analysis using
kilometres after projection and centring.

The original columns, original training/validation roles and site identifiers are
retained in `revision/data/jura-sites.csv` and `jura-data.rds`.
`revision/R/prepare_jura_data.R` regenerates the export from `gstat::jura`.
The saved export was prepared with gstat 2.1.1.

- [gstat Jura documentation](https://r-spatial.github.io/gstat/reference/jura.html)
- Atteia, O., Dubois, J.-P. and Webster, R. (1994). Geostatistical analysis of soil
  contamination in the Swiss Jura. *Environmental Pollution*, 86, 315–327.
- Goovaerts, P. (1997). *Geostatistics for Natural Resources Evaluation*.
  Oxford University Press, Appendix C.

The data retain their upstream attribution and terms; the project-code licence
does not relicense third-party data or software. Simulation responses are generated
by the included scripts with explicit seeds. No PM2.5 or bird data are used by the
current reproduction pipeline.
