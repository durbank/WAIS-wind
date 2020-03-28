# Script to import REMA surface elevation, slope, and aspect
# for each accumulation trace location

# Import all required modules and variables
from myModule import *

# Import PAIPR-generated data
PAIPR_dir = ROOT_DIR.joinpath('data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)


accum_long = format_PAIPR(data_0, start_yr=1979, end_yr=2009)
traces = accum_long.groupby('trace_ID')

# Create gdf of trace locations and mean accum
accum_mu = gpd.GeoDataFrame(
    traces.mean().drop(['Year', 'Lat', 'Lon'], axis=1),
    geometry=gpd.points_from_xy(
        traces.aggregate(np.mean)['Lon'], 
        traces.aggregate(np.mean)['Lat']), 
    crs="EPSG:4326")
accum_mu = accum_mu.to_crs("EPSG:3031")

# Define path to index shapefile
index_file = ROOT_DIR.joinpath("data/REMA", 
    "REMA_Tile_Index_Rel1.1", "REMA_Tile_Index_Rel1.1.shp")

# Import shapefile of DEM tile locations
dem_index = gpd.read_file(index_file)

# Keep only DEMs that contain accum traces
dem_index = (
    gpd.sjoin(dem_index, accum_mu, op='contains')
    .iloc[:,0:dem_index.shape[1]]).drop_duplicates()

# Find and download missing REMA DSM tiles
REMA_outdir = ROOT_DIR.joinpath("data/REMA/tiles_8m_v1.1")
tiles_list = pd.DataFrame(
    dem_index.drop(columns='geometry'))
# For development purposes, only using first 2 tiles
# get_REMA(tiles_list, REMA_outdir) 
get_REMA(tiles_list[:2], REMA_outdir) 