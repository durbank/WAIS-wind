# Script to extract REMA surface topography values at 
# all points of interest and save results

# Import modules required for script
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from myModule import *

# Import PAIPR-generated data
PAIPR_dir = ROOT_DIR.joinpath('data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)

# Format accumulation data
accum_long = format_PAIPR(
    data_0, start_yr=1979, end_yr=2009).drop('elev', axis=1)
traces = accum_long.groupby('trace_ID')

# Create df for mean annual accumulation
accum_trace = traces.aggregate(np.mean).drop('Year', axis=1)
accum_trace = gpd.GeoDataFrame(
    accum_trace, geometry=gpd.points_from_xy(
        accum_trace.Lon, accum_trace.Lat), 
    crs="EPSG:4326").drop(['Lat', 'Lon'], axis=1)

# Import Antarctic outline shapefile
ant_path = ROOT_DIR.joinpath(
    'data/Ant_basemap/Coastline_medium_res_polygon.shp')
ant_outline = gpd.read_file(ant_path)

# Convert accum crs to same as Antarctic outline
accum_trace = accum_trace.to_crs(ant_outline.crs)



# Define path to index shapefile
index_file = ROOT_DIR.joinpath("data/REMA", 
    "REMA_Tile_Index_Rel1.1", "REMA_Tile_Index_Rel1.1.shp")

# Import shapefile of DEM tile locations
dem_index = gpd.read_file(index_file)

# Keep only DEMs that contain accum traces
dem_index = (
    gpd.sjoin(dem_index, accum_trace, op='contains')
    .iloc[:,0:dem_index.shape[1]]).drop_duplicates()

# Find and download missing REMA DSM tiles
REMA_outdir = ROOT_DIR.joinpath("data/REMA/tiles_8m_v1.1")
tiles_list = pd.DataFrame(
    dem_index.drop(columns='geometry'))
get_REMA(tiles_list, REMA_outdir) 
# get_REMA(tiles_list[:2], REMA_outdir) 


# Generate list of paths to downloaded DEMs for topo
#  calculations
REMA_dir = ROOT_DIR.joinpath("data/REMA/tiles_8m_v1.1")
dem_list = [path for path in REMA_dir.glob('**/*dem.tif')]

# Calculate slope and aspect for each DEM
[calc_topo(dem) for dem in dem_list]
# [calc_topo(dem) for dem in dem_list[:1]]

# Extract elevation, slope, and aspect values for each trace 
# location
raster_vals = accum_trace.drop(['accum', 'std'], axis=1)
tif_dirs = [path.parent for path in dem_list[:1]]
for path in tif_dirs:
    raster_vals = topo_vals(path, raster_vals)

# Save raster_vals to GeoJSON
file_out = ROOT_DIR.joinpath('data/REMA/raster_vals.geojson')
raster_vals.to_file(file_out, driver='GeoJSON')