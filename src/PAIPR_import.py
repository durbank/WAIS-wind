# Script to import and format PAIPR .csv files to Python

# Import required modules
import os
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import richdem as rd
from pathlib import Path
import rasterio as rio

# Set current working directory to root of project repo
dirname=os.path.dirname
root_path = dirname(dirname(__file__))
os.chdir(root_path)

# Function to import and concatenate PAIPR .csv files
def import_PAIPR(input_dir):
    """
    Function to import PAIPR-derived accumulation data.
    This concatenates all files within the given directory into a single pandas dataframe.

    Parameters:
    input_dir (str): Path of directory containing .csv files of PAIPR results (gamma-distributions fitted to accumulation curves). This directory can contain multiple files, which will be concatenated to a single dataframe.
    """
    data = pd.DataFrame()
    for file in glob.glob(os.path.join(input_dir, "*.csv")):
        f_path = os.path.join(input_dir, file)
        data_f = pd.read_csv(f_path)
        data = data.append(data_f)
    return data

# Import PAIPR-generated data
PAIPR_dir = os.path.join(os.getcwd(), 'data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)

# Remove time series with data missing from period
# of interest (and clip to period of interest)
traces = data_0.groupby(['Lat', 'Lon', 'elev'])
data = data_0.assign(trace_ID = traces.ngroup())
traces = data.groupby('trace_ID')
data = traces.filter(
    lambda x: min(x['Year']) <= 1979 
    and max(x['Year']) >= 2009)
data = data.query("Year >= 1979 & Year < 2009")

# Ensure each trace has only one time series 
# (if not, take the mean of all time series)
data = data.groupby(['trace_ID', 'Year']).mean()

# Generate descriptive statistics based on imported 
# gamma-fitted parameters
alpha = data['gamma_shape']
alpha.loc[alpha<1] = 1
beta = 1/data['gamma_scale']
mode_accum = (alpha-1)/beta
var_accum = alpha/beta**2

# New df (in long format) with accum data assigned
accum_long = (
    data.filter(['trace_ID', 'Lat', 'Lon', 'elev', 'Year'])
    .assign(accum = mode_accum, std = np.sqrt(var_accum))
    .reset_index()
)
traces = accum_long.groupby('trace_ID')

# New accum and std dfs in wide format
accum = accum_long.pivot(
    index='trace_ID', columns='Year', values='accum')
accum_std = accum_long.pivot(
    index='trace_ID', columns='Year', values='std')

# Create df for mean annual accumulation
accum_mu = traces.aggregate(np.mean).drop('Year', axis=1)
accum_mu = gpd.GeoDataFrame(
    accum_mu, geometry=gpd.points_from_xy(
        accum_mu.Lon, accum_mu.Lat)
).drop(['Lat', 'Lon', 'elev'], axis=1)
accum_mu.crs = "EPSG:4326"

# Create a gdf from accum df
accum_gdf = gpd.GeoDataFrame(
    accum, geometry=gpd.points_from_xy(
        traces.aggregate(np.mean)['Lon'], 
        traces.aggregate(np.mean)['Lat'])
)
accum_gdf.crs = "EPSG:4326"

# Import Antarctic outline shapefile
ant_path = os.path.join(
    os.getcwd(), 
    'data/Ant_basemap/Coastline_medium_res_polygon.shp')
ant_outline = gpd.read_file(ant_path)

# Convert accum crs to same as Antarctic outline
accum_gdf = accum_gdf.to_crs(ant_outline.crs)
accum_mu = accum_mu.to_crs(ant_outline.crs)

# Plot inset map
base = ant_outline.plot(color='grey', edgecolor='black')
accum_gdf.sample(n=1000).plot(ax=base, color='red')

# Plot mean accumulation spatial map
accum_mu.sample(n=2500).plot(column='accum')

## Import required DSM data


# Define REMA DEM directory and path to index shapefile
REMA_dir = Path("data/REMA/")
index_file = (
    REMA_dir / "REMA_Tile_Index_Rel1.1" / 
    "REMA_Tile_Index_Rel1.1.shp")

# Import shapefile of DEM tile locations
dem_index = gpd.read_file(index_file)

# Keep only DEMs that contain accum traces
dem_index = (
    gpd.sjoin(dem_index, accum_mu, op='contains')
    .iloc[:,0:dem_index.shape[1]]).drop_duplicates()


import requests
import shutil

def get_REMA(tile_idx, output_dir):
    """

    """

    for idx, row in tile_idx.iterrows():
        f_dir = output_dir.joinpath(row.tile)

        if not f_dir.exists():
            f_dir.mkdir(parents=True)
            zip_path = f_dir.joinpath('tmp.tar.gz')
            r = requests.get(row.fileurl, stream=True)
            print(f"Downloading tile {f_dir.name}")
            with open(zip_path, 'wb') as zFile:
                for chunk in r.iter_content(
                        chunk_size=1024*1024):
                    if chunk:
                        zFile.write(chunk)
            # open(zip_path, 'wb').write(tile.content)
            print(f"Unzipping tile {f_dir.name}")
            shutil.unpack_archive(zip_path, f_dir)
            os.remove(zip_path)
        else:
            print(f"REMA tile {f_dir.name} already exists locally, moving to next download")


output_dir = Path(
    '/home/durbank/Downloads/REMA/tiles_8m_v1.1')
tile_index = pd.DataFrame(
    dem_index.drop(columns='geometry')).iloc[:2]
get_REMA(tile_index, output_dir)



def calc_topo(dem_path):
    """
    Calculates slope and aspect from given DEM and saves output.
    The function checks to see whether a slope/aspect file has already been created so as to avoid needless processing.
    
    Parameters:
    dem_path (pathlib.PosixPath): The relative or absolute path to an input DEM file.

    Dependencies: 
    richdem module
    GDAL binaries
    pathlib module
    """
    slope_path = Path(
        str(dem_path).replace("dem", "slope"))
    aspect_path = Path(
        str(dem_path).replace("dem", "aspect"))

    if ((not slope_path.is_file()) or 
            (not aspect_path.is_file())):
        
        dem = rd.LoadGDAL(str(dem_path))

    if not slope_path.is_file():
        slope = rd.TerrainAttribute(
            dem, attrib='slope_riserun')
        rd.SaveGDAL(str(slope_path), slope)
    
    if not aspect_path.is_file():
        aspect = rd.TerrainAttribute(dem, attrib='aspect')
        rd.SaveGDAL(str(aspect_path), aspect)
        
    




def topo_vals(tile_dir, trace_locs):
    """
    Extracts elevation, slope, and aspect values at given locations.
    
    Parameters:
    tile_dir (pathlib.PosixPath): The relative or absolute path to a directory containing REMA tile DSM, slope and aspect geotiffs.
    trace_locs (geopandas.geodataframe.GeoDataFrame): A geodataframe containing the locations at which to extract raster data. These data should have the same coordinate reference system as the raster data, with the geometries stored in a column named "geometry".

    Dependencies: Requires the rasterio (as rio) module and, by extension, GDAL binaries.
    Requires the geopandas module.
    """
    trace_locs = (
        trace_locs.assign(elev=None)
        .assign(slope=None).assign(aspect=None))

    coords = (
        [(x,y) for x, y in zip(
            trace_locs.geometry.x, trace_locs.geometry.y)]
    )
    
    # Extract elevation values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*dem.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.elev[tile_mask] = tile_vals[tile_mask]
    src.close()

    # Extract slope values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*slope.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.slope[tile_mask] = tile_vals[tile_mask]
    src.close()

    # Extract aspect values for all points within tile
    tile_path = [
        file for file in tile_dir.glob("*aspect.tif")][0]
    src = rio.open(tile_path)
    tile_vals = np.asarray(
        [x[0] for x in src.sample(coords, masked=True)])
    tile_mask = ~np.isnan(tile_vals)
    trace_locs.aspect[tile_mask] = tile_vals[tile_mask]
    src.close()

    return trace_locs




test_dem = Path(
    'data/REMA/tiles_8m_v1.1/26_20/26_20_8m_dem.tif')
test = topo_vals(test_dem.parent, accum_mu.sample(n=1000))



