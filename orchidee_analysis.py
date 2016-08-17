from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm


emis_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')             
BA_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_burntArea.nc', 
                'r', format = 'NETCDF4')                              
grid_ORCHIDEE = Dataset('../../model_data/HalfDegree-gridarea-8975-inverted.nc', 
                    'r', format = 'NETCDF4')
landCover_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_landCoverFrac.nc', 
                'r', format = 'NETCDF4')                                  
                    
# ORCHIDEE has no useful separate time data, must use JSBACH data which matches it.             
time_data = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')


#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, 
                        landCover_data, keep_time=False):
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year:int(year+month_period/12)]
    landCover = np.array(landCover)
    landCover = np.divide(landCover,100.)
    landCover = np.repeat(landCover, 12, axis=0)    
    
    BA = BA_data["burntArea"][year*12:year*12+month_period]
    BA = np.multiply(landCover, BA)
    
    
    BA = np.sum(BA, axis=1)
    BA[BA>1e5]=0.
    BA = np.multiply(BA, grid_data["cell_area"])
    
    if keep_time:
        return BA
    BA = np.sum(BA, axis=0)
    return BA
    
def get_global_BA_yearly(year, BA_data, grid_data, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data, landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#

                
def get_grid_emissions(year, month_period, emis_data, grid_data, 
                    landCover_data, time_data, keep_time=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(time_data["time"]):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
        else:
            days_per_month.append(31)
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year:int(year+month_period/12.)]
    landCover = np.array(landCover)
    landCover = np.divide(landCover,100.)
    landCover = np.repeat(landCover, 12, axis=0)

    emis = emis_data["fFirepft"][time:time+month_period]
    emis = np.array(emis)
    emis = np.multiply(emis, landCover)
    emis = np.sum(emis, axis=1)
    emis = np.multiply(emis, grid_data["cell_area"])
    emis = np.multiply(emis, sec_per_month[:, np.newaxis, np.newaxis])
    emis[emis==np.inf]=0.
    
    if keep_time:
        return emis
    emis = np.sum(emis, axis = 0)
    return emis

def get_global_emissions_yearly(year, emis_data, grid_data, landCover_data, time_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, 
                                        landCover_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year, month_period, emis_data, BA_data, 
                            landCover_data, time_data, monthly=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(time_data["time"]):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
        else:
            days_per_month.append(31)
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore division by zero warning. Returns inf or NaN.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year:int(year+month_period/12.)]
    landCover = np.array(landCover)
    #landCover[landCover>100] = 0
    landCover = np.divide(landCover,100.)
    landCover = np.repeat(landCover, 12, axis=0)

    BA = np.array(BA_data["burntArea"][time:time+month_period])
    BA = np.multiply(landCover, BA)

    # Assume fractional, add up pft dependency.
    BA = np.sum(BA, axis=1)
    if not monthly:
        BA=np.sum(BA,axis=0)
    
    BA[BA>1e5]=0.
    inv_BA = 1./BA
    # Remove infinities.
    inv_BA[inv_BA==np.inf] = 0.
    
    emis = emis_data["fFirepft"][time:time+month_period]
    emis = np.array(emis)
    emis= np.multiply(emis, landCover)
    # Add up pft dependency.
    emis = np.sum(emis, axis=1)
    emis = np.multiply(emis,sec_per_month[:, np.newaxis, np.newaxis])
    if not monthly:
        emis = np.sum(emis, axis=0)
    emis[emis==np.inf]=0.
    
    fuel_consumption = np.multiply(emis, inv_BA)
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, landCover_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, landCover_data, time_data) 
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, landCover_data, time_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,landCover_data,time_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,landCover_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC  

