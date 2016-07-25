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

def get_grid_burnt_area(year, month_period, data, grid_data):
    time = int(year*12)
    BA = np.array(data["BA"][time:time+month_period])
    # Convert to fractional.
    BA = np.divide(BA,100.)
    
    burnt_area_data = np.multiply(BA, grid_data["grid_cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, data, grid_data):
    BA_grid = get_grid_burnt_area(year, 12, data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, data, grid_data):
    years = int(len(data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1997 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,data,grid_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='GFED Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

    
#plot_global_BA_yearly(16,data_GFED, grid_GFED)
#print get_global_BA_yearly(10, data_GFED, grid_GFED)


#
# Carbon Emissions Analysis
#


def get_grid_emissions(year_start, month_period, data, grid_data):
    time = int(year_start*12)
   
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    emis_data = np.array(data["C"][time:time+month_period])
    # Convert from g to kg.
    emis_data = np.divide(emis_data, 1000.)
    
    emis_per_month = np.multiply(emis_data, grid_data["grid_cell_area"])
    emissions = np.sum(emis_per_month, axis = 0)
    return emissions

def get_global_emissions_yearly(year, data, grid_data):
    emissions_grid = get_grid_emissions(year, 12, data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions
   

def plot_global_emissions_yearly(no_years, data, grid_data):
    years = int(len(data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1997 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, data, grid_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='GFED Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
        

#plot_global_emissions_yearly(17,data_GFED,grid_GFED)
#print get_global_emissions_yearly(10, data_GFED,grid_GFED)


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, data, monthly=False):
    time = int(year_start*12)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    BA = np.array(data["BA"][time:time+month_period])
    # Convert to fractional.
    BA = np.divide(BA,100.)
    if not monthly:
        BA = np.sum(BA, axis=0)
    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    emis = data["C"][time:time+month_period]
    # Convert from g to kg.
    emis = np.divide(emis, 1000)
    if not monthly:
        emis = np.sum(emis, axis=0)
    
    fuel_consumption = np.multiply(emis, inverse_BA)
    if monthly:
        fuel_consumption = np.sum(FC_data, axis = 0)
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
    
    
def plot_global_mean_FC_yearly(no_years, data, grid_data):
    years = int(len(data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1997 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x, data, grid_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='GFED Results')
    plt.ylabel('Fuel Consumption ($kg\, C/m^2 \,burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    

#print get_global_mean_FC_yearly_rough(10,data_GFED,grid_GFED)
#plot_global_mean_FC_yearly(17,data_GFED,grid_GFED)
