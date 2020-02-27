# Script to import and format PAIPR .csv files to Python

# Import required modules
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set current working directory to root of project repo
os.chdir("../")

# Function to import and concatenate PAIPR .csv files
def import_PAIPR(input_dir):
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
# data_old = data_0
# data_0 = data_0[0:1000]
traces = data_0.groupby(['Lat', 'Lon', 'elev'])
data = data_0.assign(ID = traces.ngroup())
traces = data.groupby('ID')
data = traces.filter(
    lambda x: min(x['Year']) <= 1979 
    and max(x['Year']) >= 2009)
data = data.query("Year >= 1979 & Year < 2009")

# Ensure each trace has only one time series 
# (if not, take the mean of all time series)
data = data.groupby(['ID', 'Year']).mean()

# Generate descriptive statistics based on imported 
# gamma-fitted parameters
alpha = data['gamma_shape']
alpha.loc[alpha<1] = 1
beta = 1/data['gamma_scale']
mode_accum = (alpha-1)/beta
var_accum = alpha/beta**2

# New df (in long format) with accum data assigned
accum = (
    data.filter(['ID', 'Lat', 'Lon', 'elev', 'Year'])
    .assign(accum = mode_accum, std = np.sqrt(var_accum))
    .reset_index()
)