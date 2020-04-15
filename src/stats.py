# Script to run multilinear and spatial regression on accumulation results

# Import modules
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
import time
from myModule import *
import statsmodels.api as SM
import statsmodels.formula.api as sm

# Import PAIPR-generated data
PAIPR_dir = ROOT_DIR.joinpath('data/gamma_20111109')
data_0 = import_PAIPR(PAIPR_dir)

# Format accumulation data
accum_long = format_PAIPR(
    data_0, start_yr=1979, end_yr=2009).drop('elev', axis=1)
traces = accum_long.groupby('trace_ID')

# New accum and std dfs in wide format
accum = accum_long.pivot(
    index='Year', columns='trace_ID', values='accum')
accum_std = accum_long.pivot(
    index='Year', columns='trace_ID', values='std')

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


##### Estimate time series regressions

# Preallocate arrays for linear regression
lm_data = accum.transpose()
std_data = accum_std.transpose()
coeff = np.zeros(lm_data.shape[0])
std_err = np.zeros(lm_data.shape[0])
p_val = np.zeros(lm_data.shape[0])
R2 = np.zeros(lm_data.shape[0])

# Robust OLS model
tic = time.perf_counter()
for i in range(lm_data.shape[0]):
    X = SM.add_constant(lm_data.columns)
    y = lm_data.iloc[i]
    model = SM.RLM(y, X, M=SM.robust.norms.HuberT())
    results = model.fit()
    coeff[i] = results.params[1]
    std_err[i] = results.bse[1]
    p_val[i] = results.pvalues[1]
    # R2_r[i] = results.rsquared
toc = time.perf_counter()
print(f"Execution time of RLS: {toc-tic} s")

# Add regression results to gdf
accum_trace['trnd'] = coeff
accum_trace['p_val'] = p_val

# Add %change data (relative to long-term mean accumulation)
accum_trace['trnd_perc'] = accum_trace.trnd/accum_trace.accum
accum_trace['std_err'] = std_err / accum_trace.accum

# Create normal df of results
accum_df = pd.DataFrame(accum_trace.drop(columns='geometry'))
accum_df['East'] = accum_trace.geometry.x
accum_df['North'] = accum_trace.geometry.y

##### Large-scale multi-linear regression models
 # NOTE: Will want to incoporate surface topography 
 # and wind vector data when available

# Purely positional predictor (add elevation data later)
accum_lm = sm.ols(
    formula="accum ~ East + North", 
    data=accum_df).fit()
print(accum_lm.summary())

# Positional predictor for abs. trend data
trend_lm_abs = sm.ols(
    formula='trnd ~ East + North + accum', 
    data=accum_df).fit()
print(trend_lm_abs.summary())

# Positional predictor for rel. trend data
trend_lm_rel = sm.ols(
    formula='trnd_perc ~ East + North + accum', 
    data=accum_df).fit()
print(trend_lm_rel.summary())

##### Geographically weighted regression section

from mgwr.sel_bw import Sel_BW
from mgwr.gwr import GWR, MGWR
import libpysal as ps

