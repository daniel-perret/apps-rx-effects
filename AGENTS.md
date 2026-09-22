# apps-rx-effects — Project Memory

## Project Overview
Modeling canopy severity outcomes of prescribed burns in VA/WV forests. The pipeline runs from raw RAWS weather + burn history data through fire weather summarization, terrain processing, and a random forest classifier predicting post-fire canopy severity class.

## Code Pipeline
| Script | Role |
|---|---|
| `00_setup.R` | Libraries and global ggplot theme |
| `0.2.1_raws_processing.R` | Standardizes RAWS xlsx files → `RAWS_all_stations_092126.csv` |
| `0.2.2_raws_timeseries.R` | Builds per-fire weather summaries from RAWS timeseries |
| `0.2.4_topo_terrain.R` | Builds `topo_stack.tif` from CHILI, TPI, slope, DEM |
| `0.2.3_terrain_roughness.R` | Multi-scale rugosity/TRI metrics from topo_stack |
| `0.3.1_exploratory_topo.R` | EDA: topographic predictors vs. severity class |
| `0.4_exploratory_model.R` | Random forest classifier (main model script) |
| `0.5_prelim_prediction.R` | Applies rf1 to generate pixel-level severity predictions per burn unit |
| `0.6_prelim_leaflet.R` | Interactive leaflet map of predictions |

## Key Data Files

### Burn outcome polygons
- **`data/from_jean/080726_update/Rxb units and canopy layers/rxbmodel/1burn_cnpy.shp`**
- **Critical structure**: multiple sub-polygons per burn unit — one per canopy severity class (0 = Closed, 1 = Open, 2 = Early Seral). `Severity_1` is the outcome variable.
- Key columns: `burn_name`, `Severity_1`, `Burn_Yea_1`
- **When doing raster extractions per burn unit**, always dissolve to one polygon per unit first — do NOT extract against raw `rx.shp`:
  ```r
  rx_units <- rx.shp %>%
    group_by(burn_name) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
  ```
- Then use `terra::extract(..., fun = sd/mean/etc.)` to keep summarization inside terra and avoid materializing large pixel vectors in R.

### Burn history table
- **`data/Rxb_unit_history_092126.csv`** → loaded as `unit.hist`
- Key columns: `Fire_name_fixed` (matches `burn_name` in spatial data), `date_fixed` (format `%m/%d/%Y`), `keep.`, `burn_entry`, `RAWS_station`, `RAWS1styr`, `ACRES`
- Always filter to `keep. == "yes"` and `burn_entry == "1"` for analysis
- Date range: 2002–2026

### RAWS weather data
- **`data/RAWS/RAWS_all_stations_092126.csv`** → loaded as `all.raws`
- 11 stations. Two modern xlsx schemas (A/B) plus one legacy schema (Wise).
- **Wise station** uses a legacy 4-row header and a single combined `YYMMDDhhmm` datetime column (no separate date/time). Handled in `0.2.1`.
- **Wise and Flatwoods are treated as one continuous station** in `0.2.2` — Wise is renamed to Flatwoods on load (`ifelse(station == "Wise", "Flatwoods", station)`). Their date ranges are perfectly complementary: Wise 1998-05-18 → 2016-05-24, Flatwoods 2016-05-24 → 2026-04-24.
- RAWS station column is uppercased in `all.raws` for joining; `UPPER TRACT  WEST` is recoded to `UPPER TRACT`.

### Fire danger variables
- **`data/RAWS/fire_danger_vars.csv`** → loaded as `fire.vars`
- Date range: 2005-01-01 to 2026-06-29. Note: burn history starts 2002, so pre-2005 burns have no matching `fire.vars` data.
- Station column is named `StationName` (not `station`).

### Fire weather summary
- **`data/fire_weather_summaries/fireweather_summary_Rx1b_082126.csv`** → loaded as `final.weather`
- Filter to `date_confirm == "y"`. Joins to spatial data by `fire_name` ↔ `burn_name`.

### Terrain rasters
- **`data/topography/topo_stack.tif`** — 4-layer SpatRaster, layer names: `"chili"`, `"tpi"`, `"slp.dg"`, `"elev"`. All layers resampled to CHILI resolution/extent.
- **CHILI**: Continuous Heat-Load Index, topographic moisture/insolation proxy.
- **TPI**: Topographic Position Index — signed relative-elevation index. Can be near zero or negative; CV (sd/mean) is unreliable. Use SD for within-unit roughness: `terra::extract(tpi_layer, vect(rx_units), fun = sd, na.rm = TRUE)`.
- Slope source: LANDFIRE LF2020_SlpD; DEM source: LANDFIRE LF2020_Elev. Both projected and resampled to CHILI grid.

### Winter NDVI
- **`data/winter_NDVI/NDVI_bins_mosaic.tif`** → `eric_stack`
- One layer per year, named `y_YYYY`. Match to burn year − 1 to get prior-year NDVI.

### Canopy cover (RAP)
- **`data/RAP/TreeCover_mosaic.tif`** → `cc_stack`
- One layer per year, named `YYYY_TRE_YYYY`. Match to burn year − 1.

### RAWS station metadata
- **`data/RAWS/fw13/RAWSfw13list.csv`** → `rawsMeta`
- Columns include `Name`, `LatDegrees`, `LonDegrees`, `Elevation`. Trailing whitespace in `Name` — trim on load.

### AOI
- **`data/aoi.shp`** — project area of interest, used for cropping/masking rasters.

## Random Forest Model (0.4)
- **Outcome**: `Severity_1` (factor: Closed / Open / Early Seral)
- **Spatially balanced sample**: 1000 random points per severity class from `rx.shp`, filtered to burns in `final.weather`.
- **Point-level predictors** (extracted at sample points): `chili`, `tpi`, `slp.dg`, `elev`, `ndvi_winter_prior_year`, `canopy_cover_prior_year`
- **Unit-level predictors** (joined by `burn_name`): `ACRES`, `tpi_sd`, all fire weather variables (`daily_max_temp`, `daily_minRH`, `daily_mean_temp`, `daily_meanRH`, `days_since_precip`, `days_since_precip_025`, `ts_MAP`, `daily_mean_temp_aft`, `daily_minRH_aft`, `X1HrFM`, `X10HrFM`, `X100HrFM`, `X1000HrFM`, `HerbFM`, `WoodyFM`, `KBDI`, `ERC`, `rh_trend`, `temp_trend`, `KBDI_trend`, `ERC_trend`), `burn_doy`
- `rf_data` is `model_data |> na.omit() |> as.data.frame()` — must be a plain data.frame (not tibble) for `partialPlot()` compatibility.
- `tpi_sd`: SD of TPI pixels across each dissolved burn unit polygon. Computed in section 1b of 0.4 by dissolving `rx.shp` by `burn_name` and using `terra::extract(..., fun = sd)`.
