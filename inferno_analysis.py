from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_INFERNO = Dataset('../../model_data/Inferno_S1_fFirepft.nc', 
                'r', format = 'NETCDF4')
BA_INFERNO = Dataset('../../model_data/Inferno_S1_burntArea.nc', 
                'r', format = 'NETCDF4')              
grid_INFERNO = Dataset('../../model_data/Inferno_grid.nc', 
                    'r', format = 'NETCDF4')
landmask_INFERNO = Dataset('../../model_data/CRU-NCEP-LandMask.nc',
                    'r', format = 'NETCDF4')
landCover_INFERNO = Dataset('../../model_data/Inferno_S1_LandCoverFrac.nc',
                    'r', format = 'NETCDF4')
                

 
#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, 
                    landmask, landCover_data, keep_time=False):
    np.seterr(over='ignore')
    time = int(year*12)
    sec_per_month = []
    for i in range(month_period):
        if time+i == 0:
            sec_per_month.append(BA_data["time"][0])
        else:
            sec_per_month.append(BA_data["time"][time+i]-BA_data["time"][time+i-1])
  
    sec_per_month= np.array(sec_per_month)
    
    # Remove fill values.
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.] = 0
    
    landCover = landCover_data["LandCoverFrac"][time:time+month_period, :9]
    landCover= np.multiply(landCover,landmask)
    
    # Possibly convert to decimals and multiply by landCover (?)
    BA = BA_data["burntArea"][time:time+month_period]
    BA[BA<0.]=0
    BA = np.multiply(BA, landCover)
    BA = np.multiply(BA, grid_data["cell_area"])
    BA = np.sum(BA, axis=1)
    BA = np.multiply(BA, sec_per_month[:, np.newaxis, np.newaxis])
    
    if keep_time:
        return BA
    BA = np.sum(BA, axis=0)
    return BA
    
def get_global_BA_yearly(year, BA_data, grid_data, landmask, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data,
                                        landmask, landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#
                    
def get_grid_emissions(year, month_period, emis_data, grid_data, 
                    landmask, landCover_data, keep_time=False):
    time = int(year*12)
    sec_per_month = []
    for i in range(month_period):
        if time+i == 0:
            sec_per_month.append(emis_data["time"][0])
        else:
            sec_per_month.append(emis_data["time"][time+i]-emis_data["time"][time+i-1])      
    sec_per_month= np.array(sec_per_month)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove fill values.    
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.] = 0
    
    landCover = landCover_data["LandCoverFrac"][time:time+month_period, :9]
    landCover= np.multiply(landCover,landmask)
    
    emis = emis_data["fFirepft"][time:time + month_period]
    emis = np.multiply(emis, landCover)
    emis = np.multiply(emis, grid_data["cell_area"])
    emis = np.sum(emis, axis = 1)
    emis = np.multiply(emis, sec_per_month[:, np.newaxis, np.newaxis])
    
    if keep_time:
        return emis
    emis = np.sum(emis, axis = 0)
    return emis
    
def get_global_emissions_yearly(year, emis_data, grid_data, landmask, landCover_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landmask, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year, month_period, emis_data, BA_data, 
                            landmask, landCover_data, monthly=False):
    time = int(year*12)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    # Remove fill values.
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.] = 0
    
    landCover = landCover_data["LandCoverFrac"][time:time+month_period, :9]
    landCover= np.multiply(landCover,landmask)
    
    emis = emis_data["fFirepft"][time:time+month_period]
    emis = np.multiply(emis, landCover)
    emis = np.sum(emis,axis=1)
    if not monthly:
        emis = np.sum(emis,axis=0)
    
    BA = BA_data["burntArea"][time:time+month_period]
    BA[BA<0.]=0
    BA = np.multiply(BA, landCover)
    BA = np.sum(BA,axis=1)
    if not monthly:
        BA = np.sum(BA,axis=0)

    inv_BA = 1./BA
    # Remove infinities due to division by 0.
    inv_BA[inv_BA == np.inf] = 0
    
    fuel_consumption = np.multiply(emis, inv_BA)
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    return fuel_consumption
    
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, 
                             grid_data, landmask, landCover_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, 
                                        landmask, landCover_data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    

def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, landmask, landCover_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,landmask,landCover_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,landmask,landCover_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC
   
