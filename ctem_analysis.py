from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_CTEM = Dataset('../../model_data/CTEM_S1_fFirepft.nc', 
                'r', format = 'NETCDF4')
BA_CTEM = Dataset('../../model_data/CTEM_S1_burntArea.nc', 
                'r', format = 'NETCDF4')                
grid_CTEM = Dataset('../../model_data/CTEM-gridarea.nc', 
                    'r', format = 'NETCDF4')
landCover_CTEM = Dataset('../../model_data/CTEM_S1_landCoverFrac.nc',
                    'r', format = 'NETCDF4')
                               
#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, 
                        landCover_data, keep_time=False):
    BA = BA_data["burntArea"][year*12:year*12+month_period, :9]
    BA = np.array(BA)
    BA[BA>100]=0
    BA[BA<0.]=0
    BA = np.divide(BA,100.)
    landCover = landCover_data["landCoverFrac"][year*12:year*12+month_period]
    landCover = np.array(landCover)
    landCover[landCover>1.]=0
    BA = np.multiply(BA, landCover)
    
    BA = np.multiply(BA, grid_data["cell_area"])
    BA = np.sum(BA, axis=1)
    if keep_time:
        return BA
    BA = np.sum(BA, axis=0)
    return BA

    
def get_global_BA_yearly(year, BA_data, grid_data, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data,landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area


#
# Carbon Emissions Analysis
#

                
def get_grid_emissions(year, month_period, emis_data, grid_data, 
                    landCover_data, keep_time=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(emis_data["time"]):
            days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
        else:
            days_per_month.append(31)
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    landCover = landCover_data["landCoverFrac"][time:time+month_period]
    landCover = np.array(landCover)
    landCover[landCover>1.]=0
    complete_area_data = np.multiply(landCover, grid_data["cell_area"])
    emis = emis_data["fFirepft"][time:time+month_period, :9]
    emis = np.multiply(emis, complete_area_data)
    emis = np.sum(emis, axis = 1)
    emis = np.multiply(emis,sec_per_month[:, np.newaxis, np.newaxis])
    emis[emis<0]=0.
    if keep_time:
        return emis
    emis = np.sum(emis, axis = 0)
    return emis


def get_global_emissions_yearly(year, emis_data, grid_data, landCover_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year, month_period, emis_data, BA_data, 
                            landCover_data, monthly=False):
    time = int(year*12)
    days_per_month = []
    for i in range(month_period):
        if time+i+1 < len(emis_data["time"]):
            days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
        else:
            days_per_month.append(31)
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    landCover = landCover_data["landCoverFrac"][time:time+month_period]
    landCover = np.array(landCover)
    landCover[landCover>1.]=0
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = BA_data["burntArea"][year*12:year*12+month_period, :9]
    BA = np.array(BA)
    BA[BA>100.]=0
    BA[BA<0.]=0
    BA = np.divide(BA,100.)

    BA = np.multiply(BA, landCover)

    # Add up pft dependency.
    BA=np.sum(BA,axis=1)
    
    if not monthly:
        BA = np.sum(BA,axis=0)
    
    inv_BA = 1./BA
    # Remove infinities due to division by 0.
    inv_BA[inv_BA == np.inf] = 0.
    
    emis = emis_data["fFirepft"][time:time+month_period, :9]
    emis = np.multiply(landCover, emis)
    # Add up pft dependency.
    emis = np.sum(emis, axis=1)
    emis = np.multiply(emis, sec_per_month[:, np.newaxis, np.newaxis])
                
    if not monthly:
        emis=np.sum(emis,axis=0)
    emis[emis<0]=0.
    
    fuel_consumption = np.multiply(emis, inv_BA)  
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
    # Removing these values as they seem singularities.
    # Must ask modeller what is wrong.
    fuel_consumption[fuel_consumption>100]=0
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, landCover_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, landCover_data) 
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, landCover_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,landCover_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,landCover_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC

