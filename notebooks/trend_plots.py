# An exploratory notebook for recreating plots and 
# graphics for the spatiotemporal trends in accumulation

# Import modules
import numpy as np
import pandas as pd
import geopandas as gpd
import seaborn as sns
from pathlib import Path

# Define project root directory
ROOT_DIR = Path(__file__).parent.parent

# Import R-generated accumulation and core data
accum_long = pd.read_csv(
    ROOT_DIR.joinpath('data/R-trends/accum.csv')).drop(
        ['accum_mu', 'elev'], axis='columns')
cores_long = pd.read_csv(
    ROOT_DIR.joinpath('data/R-trends/cores.csv')).drop(
        'accum_mu', axis='columns')


core_groups = cores_long.groupby('Site')
core_sites = core_groups.mean().drop('Year', axis=1)

trace_groups = accum_long.groupby('trace_ID')
accum_sites = trace_groups.mean().drop('Year', axis=1)