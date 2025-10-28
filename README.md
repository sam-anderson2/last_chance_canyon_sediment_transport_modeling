# last_chance_canyon_sediment_transport_modeling
Modeling sediment transport in Last Chance Canyon, New Mexico, using the landlab component "riverbeddynamics"

Author: Sam Anderson
Affiliation: Tulane University / USDA–ARS Southwest Watershed Research Center
Last Updated: October 2025

This repository contains a reproducible modeling framework for coupling **Landlab’s OverlandFlow** (hydraulics) with **RiverBedDynamics** (mixed-grain sediment transport). It simulates how rainfall events mobilize sediment across variable grain-size distributions (GSDs) within the Last-Chance Creek watershed. The workflow is designed for high reproducibility, low RAM usage, and detailed post-event diagnostics.

---

## 📂 Repository Structure

### Core Scripts
| File | Description |
|------|--------------|
| **driver_expt1a.py** | Main single-storm driver. Couples OverlandFlow and RiverBedDynamics, runs one storm per DEM, and writes complete diagnostics to `run_<tag>/`. |
| **run_many_storms.py** | Batch manager for running all rainfall files in `working_climate/` across CPU cores. Handles multi-storm automation. |
| **rainfall_manager.py** | Reads rainfall hyetographs and updates rainfall intensity during the simulation loop. |
| **update_time.py** | Manages adaptive time stepping to maintain stable and efficient hydraulic solutions. |
| **save_data_to_file.py** | Utility for structured export of numeric data (CSV or TXT). |
| **save_raster.py** | Writes Landlab raster fields to ESRI ASCII grids for GIS visualization. |
| **clean_old_output.py** | Removes or resets prior `run_<tag>/` folders to avoid conflicts between experiments. |

---

### Input Data
| File/Folder | Description |
|--------------|-------------|
| **filled_lc1_dem.asc**, **filled_lc3_dem.asc** | Filled digital elevation models (ESRI ASCII format) for LC-1 and LC-3 domains. |
| **lc1_gsd_locations.asc**, **lc3_gsd_locations.asc** | Raster maps assigning nodes to upstream, intermediate, or downstream grain-size sections (coded 0/1/2). |
| **GSDs/** | Directory containing text tables of grain-size distributions (mm and cumulative %, by section). |
| **working_climate/** | Folder containing rainfall hyetographs (two-column text: `time_s`, `intensity_m/s`) used to drive model runs. |

---

### Utilities and Metadata
| File | Description |
|------|--------------|
| **.gitattributes** | Git configuration for consistent line endings and text encoding. |
| **README.md** | This documentation file describing the model structure and usage. |

---


