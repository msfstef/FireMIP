from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

# Total emissions, not per pft, because burnt area per pft is not available.                   
data_GFED = Dataset('../../model_data/GFED_DATA_1997-2013.nc', 
                'r', format = 'NETCDF4')            
grid_GFED = Dataset('../../model_data/GFED_grid.nc', 
                    'r', format = 'NETCDF4')                    
                

#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, data, grid_data, keep_time=False):
    time = int(year*12)
    BA = np.array(data["BA"][time:time+month_period])
    # Convert to fractional.
    BA = np.divide(BA,100.)
    
    BA = np.multiply(BA, grid_data["grid_cell_area"])
    if keep_time:
        return BA
    BA = np.sum(BA, axis=0)
    return BA
    
def get_global_BA_yearly(year, data, grid_data):
    BA_grid = get_grid_burnt_area(year, 12, data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#


def get_grid_emissions(year, month_period, data, grid_data, keep_time=False):
    time = int(year*12)
   
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    emis = np.array(data["C"][time:time+month_period])
    # Convert from g to kg.
    emis = np.divide(emis, 1000.)
    
    emis = np.multiply(emis, grid_data["grid_cell_area"])
    if keep_time:
        return emis
    emis = np.sum(emis, axis = 0)
    return emis

def get_global_emissions_yearly(year, data, grid_data):
    emissions_grid = get_grid_emissions(year, 12, data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year, month_period, data, monthly=False):
    time = int(year*12)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    BA = np.array(data["BA"][time:time+month_period])
    # Convert to fractional.
    BA = np.divide(BA,100.)
    if not monthly:
        BA = np.sum(BA, axis=0)
    inv_BA = 1./BA
    # Remove infinities due to division by 0.
    inv_BA[inv_BA == np.inf] = 0
    
    emis = data["C"][time:time+month_period]
    # Convert from g to kg.
    emis = np.divide(emis, 1000)
    if not monthly:
        emis = np.sum(emis, axis=0)
    
    fuel_consumption = np.multiply(emis, inv_BA)
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    return fuel_consumption
    
    
def get_global_mean_FC_yearly(year, data, grid_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["grid_cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, data, grid_data):
    total_emis = get_global_emissions_yearly(year,data,grid_data)
    total_BA = get_global_BA_yearly(year,data,grid_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC

