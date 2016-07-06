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
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, landCover_data):
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
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = np.divide(BA_data["burntArea"][time:time+month_period, :9], 100)
    BA[BA<0.] = 0
    BA[BA>1.] = 0
    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    actual_emission_data= np.multiply(landCover_data["landCoverFrac"][time:time+month_period], 
                                        emis_data["fFirepft"][time:time+month_period, :9])
    
    FC_data = np.multiply(actual_emission_data, inverse_BA)  
    FC_per_month_pft = np.multiply(FC_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    FC_per_month = np.sum(FC_per_month_pft, axis = 1)
    fuel_consumption = np.sum(FC_per_month, axis = 0)
    return fuel_consumption

    
def get_global_total_FC_yearly(year, emis_data, BA_data, landCover_data):
    FC_grid = get_grid_fuel_consumption(year, 12, emis_data, BA_data, landCover_data)
    fuel_consumption = np.sum(FC_grid)
    return fuel_consumption
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, landCover_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, landCover_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
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
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data, landCover_data):
    month_period = int(year_period*12)
    map_data = get_grid_fuel_consumption(year_start,month_period,emis_data, BA_data, landCover_data)
    lats = landCover_data["latitude"]
    lons = landCover_data["longitude"]
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
    
#print get_global_FC_yearly(300, emis_CLM, BA_CLM, landCover_CLM)
#plot_FC_map_period(140,1 , emis_CTEM,BA_CTEM,landCover_CTEM)
#plot_global_mean_FC_yearly(20, emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)
#print get_global_mean_FC_yearly(150,emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)


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

def get_global_emissions_monthly(month, emis_data, grid_data, landCover_data):
    emissions_grid = get_grid_emissions(month/12., 1, emis_data, grid_data, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions
   
    
def plot_global_emissions_monthly(no_months, emis_data, grid_data, landCover_data):
    x_data = range(len(emis_data["time"])-1)
    y_data = []
    for x in x_data[-no_months:]:
        print("%.2f" % ((x-x_data[-no_months])/
                float(len(x_data[-no_months:]))))
        y_data.append(get_global_emissions_monthly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data[-no_months:], y_data, color='r', linewidth=2.0,
               label='CTEM Results')
    plt.ylabel('Carbon Emitted ($Pg/month$)')
    plt.xlabel('Month')
    plt.legend()
    plt.show()
   

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

def plot_map_period(year_start, year_period, emis_data, grid_data, landCover_data):
    month_period = int(year_period*12)
    map_data = get_grid_emissions(year_start,month_period, emis_data, grid_data, landCover_data)
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
    cs = m.contourf(lons,lats, emis_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Total Emissions, CTEM")
    plt.show()
    
#print get_global_emissions_yearly(151,emis_CTEM,grid_CTEM, landCover_CTEM)
