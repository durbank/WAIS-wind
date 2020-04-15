# Script to generate visualizations of WAIS-wind results

# Import modules
import pandas as pd
import geopandas as gpd
import seaborn as sns
from pathlib import Path
import geoviews as gv
import holoviews as hv
from cartopy import crs as ccrs
from bokeh.io import output_notebook
output_notebook()
gv.extension('bokeh')

## Intialize script variables

# Define root of the project directory
ROOT_DIR = Path(__file__).parent.parent

# Load preprocessed results as gdf
accum_trace = gpd.read_file()

# Create non-geo df of results
accum_df = pd.DataFrame(accum_trace.drop(columns='geometry'))
accum_df['East'] = accum_trace.geometry.x
accum_df['North'] = accum_trace.geometry.y

# Define plotting projection to use
ANT_proj = ccrs.SouthPolarStereo(true_scale_latitude=-71)

# Define Antarctic boundary file
shp = str(ROOT_DIR.joinpath('data/Ant_basemap/Coastline_medium_res_polygon.shp'))

## Plot data inset map
Ant_bnds = gv.Shape.from_shapefile(shp, crs=ANT_proj).opts(
    projection=ANT_proj, width=500, height=500)
trace_plt = gv.Points(accum_trace, crs=ANT_proj).opts(
    projection=ANT_proj, color='red')
Ant_bnds * trace_plt

# Plot mean accumulation across study region
accum_plt = gv.Points(accum_trace,vdims=['accum', 'sd'], 
    crs=ANT_proj).opts(projection=ANT_proj, color='accum', 
    cmap='viridis', colorbar=True, tools=['hover'], width=800, height=500)
accum_plt

# Plot linear temporal accumulation trends
trends_insig = gv.Points(
    accum_trace[accum_trace['p.val_yr']>0.05], 
                       vdims=['coeff_perc', 'err_perc'], crs=ANT_proj).opts(
    alpha=0.05, projection=ANT_proj, color='coeff_perc', cmap='coolwarm_r', 
    symmetric=True, colorbar=True, tools=['hover'], width=800, height=500)
trends_sig = gv.Points(
    accum_trace[accum_trace['p.val_yr']<=0.05], 
                       vdims=['coeff_perc', 'err_perc'], crs=ANT_proj).opts(
    projection=ANT_proj, color='coeff_perc', cmap='coolwarm_r', symmetric=True, 
    colorbar=True, tools=['hover'], width=800, height=500)
trends_insig * trends_sig

## Display/plot summary statistics
accum_trace.describe()

sns.pairplot(
    accum_df[['accum','East','North','trnd','trnd_perc']], 
    diag_kind='kde')