from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_CLM_pft = Dataset('../../model_data/CLM_S1_fFirepft.nc', 
                'r', format = 'NETCDF4')
# Total emissions, not per pft, because burnt area per pft is not available.                   
emis_CLM = Dataset('../../model_data/CLM_S1_CFFIRE.nc', 
                'r', format = 'NETCDF4')                
BA_CLM = Dataset('../../model_data/CLM_S1_BAF.nc', 
                'r', format = 'NETCDF4')                
grid_backup_CLM = Dataset('../../model_data/CLM-gridcell.nc', 
                    'r', format = 'NETCDF4')
grid_CLM = Dataset('../../model_data/CLM-gridarea-nomask.nc', 
                    'r', format = 'NETCDF4')            
# CLM has no separate time data, must use JSBACH data which matches CLM.                  
time_data = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')



#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, time_data, keep_time=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(time_data["time"]):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
        else:
            days_per_month.append(31)
    days_per_month = np.array(days_per_month)
    
    BA = BA_data["BAF"][year*12:year*12+month_period]
    BA = np.divide(BA,100.)
    BA = np.multiply(BA,days_per_month[:,np.newaxis,np.newaxis])
    
    BA = np.multiply(BA, grid_data["cell_area"])
    if keep_time:
        return BA
    BA = np.sum(BA, axis=0)
    return BA
    
def get_global_BA_yearly(year, BA_data, grid_data, time_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data, time_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#

                   
def get_grid_emissions(year, month_period, emis_data, grid_data, time_data,keep_time=False):
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
    
    emis = emis_data["CFFIRE"][time:time+month_period]
    emis = np.multiply(emis, grid_data["cell_area"])
    emis = np.multiply(emis, sec_per_month[:, np.newaxis, np.newaxis])
    if keep_time:
        return emis
    emissions = np.sum(emis, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data, time_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#


def get_grid_fuel_consumption(year, month_period, emis_data, BA_data, time_data, monthly=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(time_data["time"]):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
        else:
            days_per_month.append(31)
    sec_per_month = np.multiply(days_per_month, 86400)
    days_per_month = np.array(days_per_month)
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    
    
    
    BA = BA_data["BAF"][year*12:year*12+month_period]
    BA = np.divide(BA,100.)
    BA = np.multiply(BA,days_per_month[:,np.newaxis,np.newaxis])
    
    if not monthly:
        BA = np.sum(BA, axis=0)
    inv_BA = 1./BA
    # Remove infinities due to division by 0.
    inv_BA[inv_BA == np.inf] = 0
    
    # Burnt area data per pft is not available, so I use CFFIRE for this calculation.
    emis = emis_data["CFFIRE"][time:time+month_period]
    emis = np.multiply(emis, sec_per_month[:, np.newaxis, np.newaxis])
    if not monthly:
        emis = np.sum(emis, axis=0)
    
    fuel_consumption = np.multiply(emis, inv_BA)
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, time_data) 
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    

def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, time_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,time_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,time_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC    

