# Ordinary Kriging Raster Production

## Purpose

Generate complete Pb, Zn and Cu continuous concentration/anomaly surfaces and corresponding kriging variance maps.

## Fixed inputs

- Sample data: `Project/Data/Suppl_material.csv`
- Reference grid: `Project/Raster/reference_grid_25m.tif`
- CRS: EPSG:32635
- Cell size: 25 m

## Selected models

- Pb: Raw concentrations, Exponential variogram
- Zn: Raw concentrations, Exponential variogram
- Cu: Raw concentrations, Gaussian variogram

## Neighborhood

- Search radius: no hard maximum distance
- Local neighborhood: nearest available samples
- Minimum neighbors: 8
- Maximum neighbors: 20

## Outputs

- Three Ordinary Kriging concentration GeoTIFFs
- Three kriging variance GeoTIFFs
- Six base-R quick-look PNG maps
- One production summary CSV

## Interpretation note

The OK concentration rasters are continuous geochemical surfaces. Formal anomaly classes or thresholds have not yet been applied. Those will be produced in the subsequent anomaly-extraction step.

## Quality-control note

A hard maximum search distance is not imposed because it created large internal NoData areas. The nearest 8-20 samples are used at each valid reference-grid cell; uncertainty in weakly supported areas is represented by the kriging-variance raster. Negative predictions, if present, are retained and counted rather than silently truncated.
