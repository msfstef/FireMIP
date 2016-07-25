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

def get_grid_burnt_area(year, month_period, BA_data, 
                    grid_data, landmask, landCover_data):
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
    
    actual_landCover_data = np.multiply(landCover_data["LandCoverFrac"][time:time+month_period, :9],landmask)
    
    # Possibly convert to decimals and multiply by landCover (?)
    BA = BA_data["burntArea"][time:time+month_period]
    BA[BA<0.]=0
    actual_BA = np.multiply(BA, actual_landCover_data)
    actual_BA = np.multiply(actual_BA, sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    burnt_area_data = np.multiply(actual_BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data, landmask, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data,
                                        landmask, landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, 
                            grid_data, landmask, landCover_data):
    years = int(len(BA_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data,landmask,landCover_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='INFERNO Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
 
    
#plot_global_BA_yearly(16,BA_INFERNO,grid_INFERNO,landmask_INFERNO,landCover_INFERNO)  
#print get_global_BA_yearly(300, BA_INFERNO,grid_INFERNO,landmask_INFERNO,landCover_INFERNO)

#
# Carbon Emissions Analysis
#
                    
def get_grid_emissions(year_start, month_period, 
                    emis_data, grid_data, landmask, landCover_data):
    time = int(year_start*12)
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
    landmask[landmask < 0.] = 0
    actual_grid_area = np.multiply(grid_data["cell_area"],landmask)
    
    landCover = landCover_data["LandCoverFrac"][time:time + month_period,:9]
    complete_area_data = np.multiply(landCover, actual_grid_area)
    
    emis_rate_data = np.multiply(emis_data["fFirepft"][time:time + month_period], complete_area_data)
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    emissions_pft = np.sum(emis_per_month, axis = 0)
    emissions = np.sum(emissions_pft, axis = 0)
    return emissions
    
def get_global_emissions_yearly(year, emis_data, grid_data, landmask, landCover_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landmask, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions


def plot_global_emissions_yearly(no_years, emis_data, grid_data, landmask, landCover_data):
    years = int(len(emis_data["time"])/12 +1)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        # Progress bar
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data, landmask, landCover_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='INFERNO Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    print "success"
    

#plot_global_emissions_yearly(20,emis_INFERNO, grid_INFERNO, landmask_INFERNO, landCover_INFERNO)
#print get_global_emissions_yearly(178, emis_INFERNO, grid_INFERNO, landmask_INFERNO, landCover_INFERNO)


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, landmask, landCover_data, monthly=False):
    time = int(year_start*12)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    # Remove fill values.
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.] = 0
    
    actual_landCover_data = np.multiply(landCover_data["LandCoverFrac"][time:time+month_period, :9],landmask)
    emis = np.multiply(actual_landCover_data, emis_data["fFirepft"][time:time+month_period])
    emis = np.sum(emis,axis=1)
    if not monthly:
        emis = np.sum(emis,axis=0)
    
    BA = BA_data["burntArea"][time:time+month_period]
    BA[BA<0.]=0
    BA = np.multiply(BA, actual_landCover_data)
    BA = np.sum(BA,axis=1)
    if not monthly:
        BA = np.sum(BA,axis=0)

    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    fuel_consumption = np.multiply(emis, inverse_BA)
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
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, 
                            grid_data, landmask, landCover_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data,landmask,landCover_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='INFERNO Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()


#print get_global_mean_FC_yearly_rough(300, emis_INFERNO,BA_INFERNO,grid_INFERNO,landmask_INFERNO, landCover_INFERNO)
#plot_global_mean_FC_yearly(20, emis_INFERNO, BA_INFERNO, grid_INFERNO,landmask_INFERNO, landCover_INFERNO)
#print get_global_mean_FC_yearly(310,emis_INFERNO,BA_INFERNO, grid_INFERNO, landmask_INFERNO, landCover_INFERNO)

   
