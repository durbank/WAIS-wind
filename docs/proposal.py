# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # DATA SCIENCE PROJET PROPOSAL
# 
# **Project Title:** Correlation of surface mass balance, surface topography, and wind vectors in the West Antarctic Ice Sheet 
# 
# **Team members**
# 
# - Durban Keeler (u1046484): durban.keeler@utah.edu
# - Eric Johnson (u0929154@utah.edu)
# %% [markdown]
# ## Background and motivation
# 
# Current and future changes in climate have important implications for the polar ice sheets. 
# Such changes are significant in part due to the outsized influence these ice sheets play in global climate mechanisms and sea levels. 
# In particular, the West Antarctic Ice Sheet (WAIS) has the potential to add an additional 4.3 meters to sea level in the coming centuries, with some sectors of the ice sheet possibly already in the early stages of collapse. 
# Improving our understanding of changes in the various mass components of the WAIS in the recent past and near future, therefore, represents an essential target for both the scientific and policymaking communities.
# 
# A key aspect of the overall ice sheet mass budget is the surface mass balance (SMB) i.e. the net inputs and outputs of mass occurring at the surface during a given period of time (traditionally one year). 
# For Antarctica, this term is largely synonymous with integrated precipitation (as snow), dominantly controlled by large synoptic-scale weather systems. 
# SMB, however, can vary drastically over small spatial and temporal scales. 
# A number of studies have linked this small-scale variability to surface topography and/or post-depositional redistribution by near-surface winds. 
# The remoteness of West Antarctica combined with the extreme environments present on the ice sheet, however, has limited the ability to collect in-situ data with sufficient coverage and density to fully investigate this proposed relationship.
# 
# %% [markdown]
# ## Project objectives
# 
# Here we propose to investigate the spatiotemporal relationships between annual SMB, surface topography, and near-surface wind vectors for a region of central West Antarctica over recent decades (~1980-2010). 
# To do this, we plan to use several new remote sensing and climate datasets. 
# We will use these datasets to determine the statistical relationships between annual SMB and surface topography and wind vectors, and how those relationships vary over the spatial and temporal extents of our data. 
# This will improve our understanding of the link between post-depositional redistribution of precipitation and spatial patterns in ice sheet mass balance.
# %% [markdown]
# ## Data
# 
# ### Antarctic surface mass balance
# 
# Mass balance data for the West Antarctic Ice Sheet (WAIS) derived from airborne ice-penetrating radar surveys (using [NASA Operation Ice Bridge Snow radar](https://nsidc.org/data/IRSNO1B/versions/2)) using the [PAIPR algorithm](https://github.com/UofU-Cryosphere/PAIPR).
# These data are already processed, so all we will need to do is importing, formatting, and cleaning.
# Additionally, permissions and use of these data (the PAIPR-derived SMB data) are not an issue, as these data are produced by Durban Keeler.
# The data covers the interior of central West Antarctica within the following bounding box: Longitude (-115.86289, -94.9989) Latitude (-79.9461, -75.58234).
# 
# ### Climate reanalysis variables
# 
# The European Centre for Medium-Range Weather Forecasts (ECMRWF) produces the [ERA5 climate reanalysis](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) (C3S, 2017). 
# This is a global 30x30 km sub-daily resolution dataset consisting of multiple atmospheric, land, and oceanic climate variables. 
# We plan to utilize the near-surface wind vectors and precipitation outputs, averaged over yearly (and possibly monthly) timescales, to investigate relationships to annual SMB.
# 
# ### Digital surface model
# 
# For surface topography data, we plan to use the [Reference Elevation Model for Antarctica (REMA)](https://www.pgc.umn.edu/data/rema/). 
# These data are 8-m resolution elevations derived from several remotely sensed datasets including, including [Cryosat-2](https://earth.esa.int/web/eoportal/satellite-missions/c-missions/cryosat-2), [ICESat](https://icesat.gsfc.nasa.gov/icesat/), and [DigitalGlobe](https://www.digitalglobe.com/products/satellite-imagery) satellite imagery (Howat et al., 2019). 
# 
# %% [markdown]
# ## Ethical considerations
# 
# Requirements to use each dataset utilized in this project are presented here:
# ERA5 data are readily available for download via API. Their conditions for service state that they agree to provide “...free and unrestricted access to meteorological data and products for educational and research purposes…”, if you meet their conditions, which we do (see the full conditions [here](https://www.ecmwf.int/sites/default/files/meteorological_data_research_education.pdf)).
# REMA data are available for download via FTP. Use of REMA data requires you to “agree to cite PGC and its sponsorship by the NSF”, by including an acknowledgement (provided in Acknowledgements section below) and citation (also provided below). It is otherwise freely available for research purposes.
# Antarctic mass balance data was derived from *in-situ* airborne radar data collected in 2010-2011 and are stored in the University of Utah Snow, Water, and Ice Research Lab. These data can be freely used and shared without issue.
# 
# %% [markdown]
# ## Data processing
# 
# These data should be fairly clean datasets overall to begin with. Data extraction will be performed through API/FTP for the DEM and climate data, and the mass balance data is already housed in-house. We will likely need to calculate slope and aspect across the DEM, and may further need to convert these values to cartesian vectors. The wind data from ERA5 will need to be downscaled to a finer resolution, likely through a simple linear interpolation, though possibly through some other statistical downscaling approach. The Antarctic mass balance data will require some filtering and formatting to remove co-located and data lacking a sufficient temporal extent, and to coerce the data to coherent standards. Additionally, special attention will need to be paid to ensure that each layer of data is correctly georeferenced so as to align correctly with the other layers. We anticipate using GeoPandas to accomplish this, though each dataset already contains spatial reference data already.
# %% [markdown]
# ## Exploratory analysis
# 
# We will initially explore the data with standard summary and descriptive statistics. 
# We will also look at georeferenced data maps of monthly averaged climatologies for precipitation (the primary driver of ice sheet mass balance), wind speed, and wind direction. 
# Development of a way to characterize surface roughness may also be important at this stage.
# %% [markdown]
# ## Data methods and modeling
# 
# In addition to standard multilinear regression, also plan to use a selection of spatial modelings (e.g. spatial lag models, geographically-weighted regression) and time-series analysis to investigate how these relationships vary over time and space.
# We may also try a variogram or Mantel test.
# %% [markdown]
# ## Project schedule
# 
# |  **Date**  | **Deadline**                                                                   |
# |------------|--------------------------------------------------------------------------------|
# |  3/08/2020 | Have all data extracted from APIs/FTPs                                         |
# | 3/15/2020  | Spring break                                                                   |
# | 3/22/2020  | Ensure data geographic overlap, data exploration, downscaling                  |
# | 3/29/2020  | Project Milestone: Functional project prototype - Initial statistical analyses |
# | 4/05/2020  | Optimization of statistical models                                             |
# | 4/12/2020  | Compile results, data visualization, project write-up                          |
# | 4/19/2020  | Project completely submitted                                                   |
# %% [markdown]
# ## Acknowledgements
# 
# DEMs provided by the Byrd Polar and Climate Research Center and the Polar Geospatial Center under NSF-OPP awards 1543501, 1810976, 1542736, 1559691, 1043681, 1541332, 0753663, 1548562, 1238993 and NASA award NNX10AN61G. Computer time provided through a Blue Waters Innovation Initiative. DEMs produced using data from DigitalGlobe, Inc. 
# %% [markdown]
# ## References
# 
# Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate . Copernicus Climate Change Service Climate Data Store (CDS), date of access. https://cds.climate.copernicus.eu/cdsapp#!/home
# 
# Howat, I. M., Porter, C., Smith, B. E., Noh, M.-J., and Morin, P.: The Reference Elevation Model of Antarctica, The Cryosphere, 13, 665-674, https://doi.org/10.5194/tc-13-665-2019, 2019.
# 

# %%


