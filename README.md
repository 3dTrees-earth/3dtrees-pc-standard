# 3DTrees LAS/LAZ Validation and Standardization Tool

This repository provides a command‑line tool (with a Docker image) to:

- **Validate and standardize single LAS/LAZ point clouds**, and  
- **Check consistency and basic quality of LAS/LAZ collections** (tile sets),

with a focus on LiDAR / TLS datasets for tree and forest applications.

The core logic is written in R version 4.5 using [`lidR`](https://github.com/r-lidar/lidR),  [`sf`](https://r-spatial.github.io/sf/), and [`jsonlite`](https://github.com/jeroen/jsonlite).  
A ready‑to‑use Docker image is defined in the `Dockerfile`.

---

## Features

### Single‑file mode (`--method single_file`)

For a single LAS/LAZ file, the tool:

- Loads the file via `lidR::readLAS()`
- Runs `lidR::las_check()` and logs errors, warnings, and messages
- Fixes some header issues:
  - If points lie outside the header bounding box, `las_update()` is applied
- Checks basic geometry / plausibility:
  - Reports extent in X, Y, and Z
  - Warns if Z range is implausible (`< -100` or `> 5500`)
- Computes point density:
  - Via a convex hull of XY coordinates
  - Warns if density is ≤ `--min-density` (points / m²)
- Checks and fixes attributes:
  - Computes per‑column statistics (`attribute_stats()`):  
    counts, NA counts, distinct values, min / max / mean / median / SD, sums, etc.
  - Validates LAS “Extra Bytes” definitions against actual R data types and hot‑fixes mismatches in the header
  - Optionally **removes attributes** by name (`--removeable_attributes`)
- Checks and handles CRS:
  - Verifies that a CRS is defined
  - Tries to transform a random sample to **EPSG:4978**
  - Builds a convex hull polygon, transforms it to **EPSG:4326**, and checks:
    - Longitude/latitude plausibility (inside world bounds)
    - Whether the bounding box centroid lies on land  
      (using `/src/geodata/ne_10m_land.shp` and `ne_10m_minor_islands.shp`)
  - If CRS is not transformable or implausible:
    - Logs an error
    - Re‑offsets X, Y, Z so the cloud starts at (0, 0, 0)
    - Assigns a **fake CRS EPSG:4978** (so that later processes still have a defined CRS)
- Harmonizes point data format:
  - If `Point Data Format ID > 7`, demotes it to **7**
  - If a `NIR` channel exists, moves it to an Extra Byte (`NIRe`) and drops the original `NIR` column
- Writes outputs:
  - A **standardized `.laz`** file
  - A **GeoJSON** convex hull in WGS84 (`*_convex_hull_wgs84.GeoJSON`) when CRS is usable
  - A **JSON log file** with LAS metadata pre‑ and post‑standardization

---

### Collection mode (`--method collection`)

For a directory containing multiple LAS/LAZ files, the tool:

1. **Scans the folder** for `*.las` and `*.laz` files.

2. For each file:
   - Reads the LAS/LAZ
   - Stores:
     - File ID and name
     - CRS (`sf::st_crs`)
     - Attribute names and R data types
     - A convex hull of XY coordinates (sampled if extremely large)
     - Number of points
     - Per‑attribute statistics via `attribute_stats()` (as above)

3. **Collection‑level checks**:

   - **CRS homogeneity**
     - Extracts CRS for each file and compares EPSG codes
     - Reports whether all tiles share the same CRS

   - **Attribute name homogeneity**
     - Checks if all tiles have the same set and order of columns
     - Computes:
       - All attribute names (union)
       - Common attribute names (intersection)
       - Non‑overlapping attribute names (present only in some tiles)

   - **Attribute content screening**
     - Identifies attributes that are:
       - **Only NA** in all tiles
       - **Always zero** (mean = 0 and SD = 0 in every tile where they are present)
     - Derives a set of **“removable” attributes**:  
       attributes that are either all‑NA or constant zero across the collection

   - **Attribute type homogeneity**
     - Compares R classes of attributes against the first tile
     - Reports whether types are consistent across the collection

   - **Tile geometry / tiling sanity checks**  
     (using convex hull polygons of each tile):
     - Builds per‑tile `sf` polygons from the convex hulls
     - Uses buffered geometry to:
       - Detect overlaps (using `st_overlaps` on slightly **shrunk** polygons)
       - Check whether tiles are (approximately) disjoint (using `st_disjoint` on slightly **expanded** polygons)

   - **Concave hull / coverage map**
     - Combines per‑tile polygons into a single `MULTIPOLYGON`
       using `polygons_to_multipolygon()`
     - Assigns the CRS of the first tile
     - If the CRS is defined and transformable, transforms the multipolygon to **EPSG:4326**
     - Exports the **multipolygon as WKT** (text), suitable for downstream use

   - **Global attribute statistics**
     - For each attribute that has **non‑NA and non‑constant** values across the collection:
       - Aggregates per‑tile statistics (`n`, sum, sum of squares)
       - Computes global mean, variance, and standard deviation  
         using `multi_tile_variance()`

4. **Writes one JSON summary** per collection:
   - `collection_summary.json` in `--output-dir`
   - Structure:
     - `collection`: collection‑level checks and stats
       - `n_tiles`, CRS homogeneity, example CRS, EPSG codes
       - Attribute name & type homogeneity
       - All / common / non‑overlapping attribute names
       - Lists of attributes that are only NA, only zeros, and suggested removable attributes
       - Tile overlap/disjoint flags
       - Global attribute statistics
       - Multipolygon WKT of the coverage
     - `files`: per‑file entries
       - `crs` (as a simple list), attributes, attribute types
       - Per‑attribute stats, convex hull coordinates, point counts, etc.

This mode does **not** modify data; it is purely for validation and summary.

---

## Installation

### 1. Using Docker (recommended)

Build the image:

```bash
docker build -t 3dtrees-lidar .