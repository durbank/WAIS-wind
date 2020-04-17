# Python script to generate and analyze WAIS wind results

# Import modules required for script
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

# # Create a gdf from accum df
# accum_gdf = gpd.GeoDataFrame(
#     accum.transpose(), geometry=gpd.points_from_xy(
#         traces.aggregate(np.mean)['Lon'], 
#         traces.aggregate(np.mean)['Lat']), 
#         crs="EPSG:4326")

# Import Antarctic outline shapefile
ant_path = ROOT_DIR.joinpath(
    'data/Ant_basemap/Coastline_medium_res_polygon.shp')
ant_outline = gpd.read_file(ant_path)

# Convert accum crs to same as Antarctic outline
accum_trace = accum_trace.to_crs(ant_outline.crs)
# accum_gdf = accum_gdf.to_crs(ant_outline.crs)



# Load ERA wind vectors
ERA_path = ROOT_DIR.joinpath(
    'data/downscaling_data_qc.npy')
ERA_raw = pd.DataFrame(
    np.load(str(ERA_path), allow_pickle=True))

# df for mean wind vectors
wind_mu = pd.DataFrame(
    {'trace_ID': ERA_raw.iloc[:,0], 
    'wind_u': ERA_raw.iloc[:,7].apply(lambda x: np.mean(x)), 
    'wind_v': ERA_raw.iloc[:,8].apply(lambda x: np.mean(x))})

# Add ERA data to main df
accum_trace = accum_trace.merge(wind_mu, on='trace_ID')




## Load R-generated slope data and add to dataframe
R_long = pd.read_csv(
    ROOT_DIR.joinpath('data/R-trends/accum.csv')).drop(
        ['accum_mu', 'elev'], axis='columns')
R_groups = R_long.groupby('trace_ID')
R_sites = R_groups.mean()[[
    'Easting', 'Northing', 'elev.REMA', 'slope', 'aspect']
    ].rename(columns={'elev.REMA': 'elev_REMA'})
R_gdf = gpd.GeoDataFrame(
    R_sites, 
    geometry=gpd.points_from_xy(R_sites.Easting, 
    R_sites.Northing), crs='EPSG:3031').drop(
        ['Easting', 'Northing'], axis='columns')

# Join R-topo data to main gdf
R_gdf['geometry'] = R_gdf.buffer(1)
accum_trace = gpd.sjoin(
    accum_trace, R_gdf, how='inner', op='within').drop(
        'index_right', axis='columns')

## Estimate time series regressions

# Simple polyfit model
lin_coeffs = np.polyfit(accum.index, accum, 1)
beta_yr = pd.Series(lin_coeffs[0], index=accum.columns)

# Preallocate arrays for linear regression
lm_data = accum.transpose()
std_data = accum_std.transpose()
coeff = np.zeros(lm_data.shape[0])
std_err = np.zeros(lm_data.shape[0])
p_val = np.zeros(lm_data.shape[0])
R2 = np.zeros(lm_data.shape[0])

# Full stats (with diagnostics) OLS model
tic = time.perf_counter()

for i in range(lm_data.shape[0]):
    X = SM.add_constant(lm_data.columns)
    y = lm_data.iloc[i]
    model = SM.OLS(y, X)
    results = model.fit()
    coeff[i] = results.params[1]
    std_err[i] = results.bse[1]
    p_val[i] = results.pvalues[1]
    R2[i] = results.rsquared
toc = time.perf_counter()
print(f"Execution time of OLS: {toc-tic}s")



# Full stats (with diagnostics) WLS model
tic = time.perf_counter()
for i in range(lm_data.shape[0]):
    X = SM.add_constant(lm_data.columns)
    y = lm_data.iloc[i]
    w = 1/(std_data.iloc[i] ** 2)
    model = SM.WLS(y, X, weights=w)
    results = model.fit()
    coeff[i] = results.params[1]
    std_err[i] = results.bse[1]
    p_val[i] = results.pvalues[1]
    R2[i] = results.rsquared
toc = time.perf_counter()
print(f"Execution time of WLS: {toc-tic}s")



# Full stats (with diagnostics) Robust OLS model
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
accum_trace['trnd_perc'] = accum_trace.trnd/accum_trace.accum
accum_trace['std_err'] = std_err / accum_trace.accum

## Large-scale spatial multi-linear regression

# Create normal df of results
accum_df = pd.DataFrame(accum_trace.drop(columns='geometry'))
accum_df['East'] = accum_trace.geometry.x
accum_df['North'] = accum_trace.geometry.y

# Diagnostic summary statistics
accum_df.describe()

# import seaborn as sns
# sns.pairplot(
#     accum_df[['accum', 'East', 'North', 'trnd', 'trnd_perc']]
#     .sample(100), diag_kind='kde')

# Large scale linear models
accum_lm = sm.ols(
    formula="accum ~ East + North + slope + wind_u + wind_v", 
    data=accum_df).fit()
accum_lm.summary()

trend_lm_abs = sm.ols(
    formula='trnd ~ East + North + accum + slope + aspect + wind_u + wind_v', 
    data=accum_df).fit()
trend_lm_abs.summary()

trend_lm_rel = sm.ols(
    formula='trnd_perc ~ East + North + accum + slope + aspect + wind_u + wind_v', 
    data=accum_df).fit()
trend_lm_rel.summary()



## Geographically weighted regression
from mgwr.sel_bw import Sel_BW
from mgwr.gwr import GWR, MGWR
import libpysal as ps

