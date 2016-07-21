from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm


emis_BLAZE = Dataset('../../model_data/LPJ-GUESS-BLAZE_SF1_Cfire.nc', 
                'r', format = 'NETCDF4')
BA_BLAZE = Dataset('../../model_data/LPJ-GUESS-BLAZE_SF1_BA.nc', 
                'r', format = 'NETCDF4')            
grid_BLAZE = Dataset('../../model_data/HalfDegree-gridarea-8950.nc', 
                    'r', format = 'NETCDF4')

# LPJ BLAZE has no useful separate time data, must use JSBACH data which matches it.             
time_data = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')


#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data):
    BA = BA_data["BA."][year*12:year*12+month_period]
    BA = np.divide(BA, 100.)
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data):
    years = int(len(BA_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ BLAZE Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

#plot_global_BA_yearly(16,BA_BLAZE, grid_BLAZE)
#print get_global_BA_yearly(300, BA_BLAZE, grid_BLAZE)

 
#
# Carbon Emissions Analysis
#


def get_grid_emissions(year_start, month_period, emis_data, grid_data, time_data):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(time_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(time_data["time"]):
                days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    emis_rate_data = np.multiply(emis_data["Cfire.monthly"][time:time+month_period], grid_data["cell_area"])
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis])
    emissions = np.sum(emis_per_month, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data, time_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions
   

def plot_global_emissions_yearly(no_years, emis_data, grid_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data, time_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ BLAZE Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()


#plot_global_emissions_yearly(20,emis_BLAZE,grid_BLAZE,time_data)
#print get_global_emissions_yearly(308, emis_BLAZE, grid_BLAZE, time_data) 


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, time_data):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(time_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(time_data["time"]):
                days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    fractional_BA = np.divide(BA_data["BA."][time:time+month_period], 100)
    inverse_BA = 1./np.array(fractional_BA)
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    FC_data = np.multiply(emis_data["Cfire.monthly"][time:time+month_period], 
                           inverse_BA)                     

    FC_per_month = np.multiply(FC_data, 
                sec_per_month[:, np.newaxis, np.newaxis])
    fuel_consumption = np.nansum(FC_per_month, axis = 0)
    return fuel_consumption
    
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, time_data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, time_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,time_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC

    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x, emis_data, BA_data, grid_data, time_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ BLAZE Results')
    plt.ylabel('Fuel Consumption ($kg\, C/m^2 \,burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()


#print get_global_mean_FC_yearly_rough(300,emis_BLAZE,BA_BLAZE,grid_BLAZE,time_data)
#print plot_global_mean_FC_yearly(20, emis_BLAZE, BA_BLAZE, grid_BLAZE, time_data)
#print get_global_mean_FC_yearly(300, emis_BLAZE, BA_BLAZE, grid_BLAZE, time_data)
       
