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

def get_grid_burnt_area(year, month_period, BA_data, grid_data, landCover_data):
    BA = BA_data["burntArea"][year*12:year*12+month_period, :9]
    BA[BA>100.]=0
    BA[BA<0.]=0
    BA = np.divide(BA,100)
    BA = np.multiply(BA, landCover_data["landCoverFrac"][year*12:year*12+month_period])
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data,landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data, landCover_data):
    years = int(len(BA_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1861 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data,landCover_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CTEM Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
 
    
#plot_global_BA_yearly(16,BA_CTEM, grid_CTEM, landCover_CTEM)
#print get_global_BA_yearly(136, BA_CTEM, grid_CTEM, landCover_CTEM)


#
# Carbon Emissions Analysis
#

                
def get_grid_emissions(year_start, month_period, emis_data, grid_data, landCover_data):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(emis_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(emis_data["time"]):
                days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
       
    complete_area_data = np.multiply(landCover_data["landCoverFrac"][time:time+month_period], grid_data["cell_area"])
    emis_rate_data = np.multiply(emis_data["fFirepft"][time:time+month_period, :9], complete_area_data)
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    emissions_pft = np.sum(emis_per_month, axis = 0)
    emissions = np.sum(emissions_pft, axis = 0)
    return emissions


def get_global_emissions_yearly(year, emis_data, grid_data, landCover_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions


def plot_global_emissions_yearly(no_years, emis_data, grid_data, landCover_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1861 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data, landCover_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CTEM Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    

#plot_global_emissions_yearly(16,emis_CTEM,grid_CTEM,landCover_CTEM)
#print get_global_emissions_yearly(151,emis_CTEM,grid_CTEM, landCover_CTEM)



#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, landCover_data, monthly=False):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(emis_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(emis_data["time"]):
                days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(emis_data["time"][time+i+1]-emis_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    landCover = landCover_data["landCoverFrac"][time:time+month_period]
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = np.divide(BA_data["burntArea"][time:time+month_period, :9], 100)
    BA=np.array(BA)
    BA[BA<0.] = 0.
    #BA[BA>1.] = 0.
    
    BA = np.multiply(landCover, BA)

    # Add up pft dependency.
    BA=np.sum(BA,axis=1)
    
    if not monthly:
        BA = np.sum(BA,axis=0)
    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0.
    
    emis = np.multiply(landCover, emis_data["fFirepft"][time:time+month_period, :9])
    # Add up pft dependency.
    emis = np.sum(emis, axis=1)
    emis = np.multiply(emis, 
                sec_per_month[:, np.newaxis, np.newaxis])
                
    if not monthly:
        emis=np.sum(emis,axis=0)
    
    fuel_consumption = np.multiply(emis, inverse_BA)  
    if monthly:
        fuel_consumption = np.sum(fuel_consumption, axis = 0)
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
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data, landCover_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1861 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data,landCover_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CTEM Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    

#print get_global_mean_FC_yearly_rough(136, emis_CTEM, BA_CTEM, grid_CTEM, landCover_CTEM)
#plot_global_mean_FC_yearly(20, emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)
#print get_global_mean_FC_yearly(150,emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)

