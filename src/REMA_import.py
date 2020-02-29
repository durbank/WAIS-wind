# Script to import REMA surface elevation, slope, and aspect
# for each accumulation trace location

# Import required modules
import geopandas as gpd
from ftplib import FTP
import xarray as xr