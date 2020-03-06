# -*- coding: utf-8 -*-

"""
Created on Mon Mar  2 14:30:08 2020

@author: Eric Johnson
"""
# xarray for netcdfs

def Extract_ERA5():

    # ERA5 Login Credentials:
    # username: u0929154@utah.edu
    # password: DataScience2020
    
    # Extent: Longitude (-115.86289, -94.9989) Latitude (-79.9461, -75.58234) 
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This section modified from https://stackoverflow.com/questions/44210656/how-to-check-if-a-module-is-installed-in-python-and-if-not-install-it-within-t/44210735
    # Check to insure all necessary packages are installed, install missing packages
    # import sys
    # import subprocess
    # import pkg_resources
    
    # required = {'cdsapi', 'netCDF4'}
    # installed = {pkg.key for pkg in pkg_resources.working_set}
    # missing = required - installed
    
    # if missing:
    #     python = sys.executable
    #     subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
    ###############################################################################
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Imports
    
    import cdsapi
    import netCDF4
    import numpy as np
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Initial conditions
    
    # Extents
    iLat_min = -75.58234
    iLat_max = -79.9461
    iLon_min = -115.86289
    iLon_max = -75.58234
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Download ERA5 variables into files
    
    # url = 'https://cds.climate.copernicus.eu/api/v2'
    # key = '38293:414bf99c-171e-48e6-b13d-75661652a8de'
    # login_file = '.cdsapirc'
    
    
    # 10-m u-component of wind
    c = cdsapi.Client()
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_u_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_u_component_of_wind.nc')
    
    
    # 10-m v component of wind
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : '10m_v_component_of_wind',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('10m_v_component_of_wind.nc')
    
    
    # Total precipitation
    r = c.retrieve(
        'reanalysis-era5-land-monthly-means', {
                'format'      : 'netcdf',
                'product_type': 'monthly_averaged_reanalysis',
                'variable'    : 'total_precipitation',
                'stream'      : 'moda',
                'year'        : [
                    '1980','1981','1982','1983','1984','1985','1986','1987','1988',
                    '1989','1990','1991','1992','1993','1994','1995','1996','1997',
                    '1998','1999','2000','2001','2002','2003','2004','2005','2006',
                    '2007','2008','2009','2010','2011','2012','2013','2014','2015'
                    ],
                'month'       : [
                    '01','02','03','04','05','06','07','08','09','10','11','12'
                    ],
                'time'        : '00:00'
        })
    r.download('total_precipitation.nc')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extract necessary portions of files
    
    # Extract lat/lon data
    data = netCDF4.Dataset('10m_u_component_of_wind.nc')
    vLat = data.variables['latitude']
    vLat = np.asarray(vLat[:])
    vLon = data.variables['longitude']
    vLon = np.asarray(vLon[:])
    vLon[vLon > 180] = -(360 - vLon[vLon > 180]) # Is this right? What's with the 360 deg latitudes?
    
    # Find indices for required extents
    iLat_min_idx = np.argmin(np.abs(vLat-iLat_min))
    iLat_max_idx = np.argmin(np.abs(vLat-iLat_max))
    iLon_min_idx = np.argmin(np.abs(vLon-iLon_min))
    iLon_max_idx = np.argmin(np.abs(vLon-iLon_max))
    
    # Read NetCDF files
    data = netCDF4.Dataset('10m_u_component_of_wind.nc')
    wind_u_10m = data.variables['u10']
    data = netCDF4.Dataset('10m_v_component_of_wind.nc')
    wind_v_10m = data.variables['v10']
    data = netCDF4.Dataset('total_precipitation.nc')
    precip = data.variables['prcp']
    
    # Extract required portions of each variable
    wind_u_10m = wind_u_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    wind_v_10m = wind_v_10m[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]
    precip = precip[:,iLat_min_idx:iLat_max_idx,iLon_min_idx:iLon_max_idx]














    return wind_u_10m, wind_v_10m, precip
























































