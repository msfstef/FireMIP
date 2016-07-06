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
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, year_period, emis_data, BA_data):
    # Using 1 year equal to 365.25 days.
    sec_per_year = 31557600.
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = np.divide(BA_data["BA"][year_start:year_start+year_period], 100)
    BA[BA>100.] = 0
    inverse_BA = np.array(1./BA)
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    actual_emission_data = emis_data["Cfire"][year_start:year_start+year_period]
    
    FC_data = np.multiply(actual_emission_data, inverse_BA)  
    FC_per_month= np.multiply(FC_data, sec_per_year)
    fuel_consumption = np.sum(FC_per_month, axis = 0)
    return fuel_consumption

    
def get_global_total_FC_yearly(year, emis_data, BA_data):
    FC_grid = get_grid_fuel_consumption(year, 1, emis_data, BA_data)
    fuel_consumption = np.sum(FC_grid)
    return fuel_consumption
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data): 
    FC_data = get_grid_fuel_consumption(year,1,emis_data, BA_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data):
    years = len(emis_data["year"])
    x_data = range(years)
    x_data_plot = [x+1901 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CTEM Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data):
    map_data = get_grid_fuel_consumption(year_start,year_period,emis_data, BA_data)
    lats = emis_data["lat"]
    lons = emis_data["lon"]
    emis_data = map_data
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, 
        urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, emis_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Carbon emitted per area burned ($kg\, C / m^2 \, burned$)")
    plt.title("Fuel Consumption, CTEM")
    plt.show()
    
#print get_global_FC_yearly(300, emis_MC2, BA_MC2, grid_MC2)
#plot_FC_map_period(100,1 , emis_MC2, BA_MC2)
#plot_global_mean_FC_yearly(20, emis_MC2, BA_MC2, grid_MC2)
#print get_global_mean_FC_yearly(100,emis_MC2, BA_MC2, grid_MC2)



#
# Carbon Emissions Analysis
#

                   
def get_grid_emissions(year_start, year_period, emis_data, grid_data):
    # Using 1 year equal to 365.25 days.
    sec_per_year = 31557600.
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    complete_area_data = grid_data["cell_area"]
    emis_rate_data = np.multiply(emis_data["Cfire"][year_start:year_start+year_period], complete_area_data)
    emis_per_year = np.multiply(emis_rate_data, sec_per_year)
    emissions = np.sum(emis_per_year, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data):
    emissions_grid = get_grid_emissions(year, 1, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions


def plot_global_emissions_yearly(no_years, emis_data, grid_data):
    years = len(emis_data["year"])
    x_data = range(years)
    x_data_plot = [x+1901 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='MC2 Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, emis_data, grid_data):
    map_data = get_grid_emissions(year_start,year_period, emis_data, grid_data)
    # Convert to billions of kg.
    map_data = np.multiply(map_data, (1./(10**9)))
    
    lats = grid_data["lat"]
    lons = grid_data["lon"]
    emis_data = map_data
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='aqua')
    cs = m.contourf(lons,lats, emis_data, 10, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Total Emissions, MC2")
    plt.show()
    
#plot_map_period(100,2,emis_MC2, grid_MC2)    
#plot_global_emissions_yearly(20, emis_MC2, grid_MC2)
#print get_global_emissions_yearly(108, emis_MC2, grid_MC2)
