from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_MC2 = Dataset('../../model_data/MC2_GlobalFire_Cfire.nc', 
                'r', format = 'NETCDF4')
BA_MC2 = Dataset('../../model_data/MC2_GlobalFire_BA.nc', 
                'r', format = 'NETCDF4')                
grid_MC2 = Dataset('../../model_data/HalfDegree-gridarea-8975.nc', 
                    'r', format = 'NETCDF4')

#
# Burnt Area Analysis
#


def get_grid_burnt_area(year, year_period, BA_data, grid_data, keep_time=False):    
    BA = BA_data["BA"][year:year+year_period]
    BA = np.divide(BA,100.)
    
    BA = np.array(BA)
    BA[BA>1.]=0
    
    BA = np.multiply(BA, grid_data["cell_area"])    
    
    if keep_time:
        last_yr = year+year_period
        if last_yr > 107:
            nan_arr = np.empty((last_yr-108,360,720))*np.nan
            BA = np.concatenate((BA,nan_arr),axis=0)
        return BA
    BA= np.sum(BA, axis=0)
    return BA

def get_global_BA_yearly(year, BA_data, grid_data):
    BA_grid = get_grid_burnt_area(year, 1, BA_data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#

def get_grid_emissions(year, year_period, emis_data, grid_data, keep_time=False):
    # Using 1 year equal to 365.25 days.
    sec_per_year = 31557600.
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    emis = emis_data["Cfire"][year:year+year_period]
    emis = np.multiply(emis,grid_data["cell_area"])
    emis = np.multiply(emis, sec_per_year)
    emis = np.array(emis)
    emis[emis>1e30]=0.
    
    if keep_time:
        last_yr = year+year_period
        if last_yr > 107:
            nan_arr = np.empty((last_yr-108,360,720))*np.nan
            emis = np.concatenate((emis,nan_arr),axis=0)
        return emis
    emis = np.sum(emis, axis = 0)
    return emis

def get_global_emissions_yearly(year, emis_data, grid_data):
    emissions_grid = get_grid_emissions(year, 1, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year, year_period, emis_data, BA_data, monthly=False):
    # Using 1 year equal to 365.25 days.
    sec_per_year = 31557600.
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = BA_data["BA"][year:year+year_period]
    BA = np.divide(BA,100.)
    BA = np.array(BA)
    BA[BA>1.] = 0
    if not monthly:
        BA = np.sum(BA, axis=0)
    
    inv_BA = np.array(1./BA)
    # Remove infinities due to division by 0.
    inv_BA[inv_BA == np.inf] = 0
    
    emis = emis_data["Cfire"][year:year+year_period]
    emis = np.multiply(emis, sec_per_year)
    emis = np.array(emis)
    emis[emis>1e30]=0
    if not monthly:
        emis = np.sum(emis,axis=0)
    
    fuel_consumption = np.multiply(emis, inv_BA)
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data): 
    FC_data = get_grid_fuel_consumption(year,1,emis_data, BA_data, grid_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 


def get_global_mean_FC_yearly_rough(year, emis_data, BA_data,grid_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC

