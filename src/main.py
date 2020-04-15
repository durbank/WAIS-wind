# Python script to generate and analyze WAIS wind results

# Import modules required for script
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
import time
from myModule import *

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








## Estimate time series regressions

import statsmodels.api as SM
import statsmodels.formula.api as sm

# # Simple polyfit model
# lin_coeffs = np.polyfit(accum.index, accum, 1)
# beta_yr = pd.Series(lin_coeffs[0], index=accum.columns)

# Preallocate arrays for linear regression
lm_data = accum.transpose()
std_data = accum_std.transpose()
coeff = np.zeros(lm_data.shape[0])
std_err = np.zeros(lm_data.shape[0])
p_val = np.zeros(lm_data.shape[0])
R2 = np.zeros(lm_data.shape[0])

# # Full stats (with diagnostics) OLS model
# tic = time.perf_counter()

# for i in range(lm_data.shape[0]):
#     X = SM.add_constant(lm_data.columns)
#     y = lm_data.iloc[i]
#     model = SM.OLS(y, X)
#     results = model.fit()
#     coeff[i] = results.params[1]
#     std_err[i] = results.bse[1]
#     p_val[i] = results.pvalues[1]
#     R2[i] = results.rsquared
# toc = time.perf_counter()
# print(f"Execution time of OLS: {toc-tic}s")



# # Full stats (with diagnostics) WLS model
# tic = time.perf_counter()
# for i in range(lm_data.shape[0]):
#     X = SM.add_constant(lm_data.columns)
#     y = lm_data.iloc[i]
#     w = 1/(std_data.iloc[i] ** 2)
#     model = SM.WLS(y, X, weights=w)
#     results = model.fit()
#     coeff[i] = results.params[1]
#     std_err[i] = results.bse[1]
#     p_val[i] = results.pvalues[1]
#     R2[i] = results.rsquared
# toc = time.perf_counter()
# print(f"Execution time of WLS: {toc-tic}s")



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

# import matplotlib.pyplot as plt
# idx = np.arange(0, len(coeff), 100)
# plt.scatter(idx, coeff[idx], color='blue', label='OLS')
# plt.scatter(idx, coeff_w[idx], color='red', label='WLS')
# plt.scatter(idx, coeff_r[idx], color='purple', label='RLS')
# plt.legend(loc='best')
# plt.show()

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
import seaborn as sns

# accum_df.describe()
# sns.pairplot(
#     accum_df[['accum', 'East', 'North', 'trnd', 'trnd_perc']]
#     .sample(100), diag_kind='kde')

# Large scale linear models
accum_lm = sm.ols(
    formula="accum ~ East + North", 
    data=accum_df).fit()
print(accum_lm.summary())

trend_lm_abs = sm.ols(
    formula='trnd ~ East + North + accum', 
    data=accum_df).fit()
print(trend_lm_abs.summary())

trend_lm_rel = sm.ols(
    formula='trnd_perc ~ East + North + accum', 
    data=accum_df).fit()
print(trend_lm_rel.summary())

## Geographically weighted regression
from mgwr.sel_bw import Sel_BW
from mgwr.gwr import GWR, MGWR
import libpysal as ps

