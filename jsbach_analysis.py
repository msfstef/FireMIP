from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm


emis_JSBACH = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')
BA_JSBACH = Dataset('../../model_data/JSBACH_SF1burntArea.nc', 
                'r', format = 'NETCDF4')                              
grid_JSBACH = Dataset('../../model_data/JSBACH_grid.nc', 
                    'r', format = 'NETCDF4')              
                    

#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data):
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
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    BA = (np.array(BA_JSBACH["burntArea"][time:time+month_period]))
    # Remove fill values.
    BA[BA>10**10] = 0
    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    

    FC_data = np.multiply(emis_data["fFirepft"][time:time+month_period], 
                           inverse_BA)                     

    FC_per_month = np.multiply(FC_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    FC_pft = np.nansum(FC_per_month, axis = 0)
    fuel_consumption = np.nansum(FC_pft, axis = 0)
    return fuel_consumption
    

def get_global_FC_monthly(month, emis_data, BA_data):
    FC_grid = get_grid_emissions(month/12., 1, emis_data, BA_data)
    fuel_consumption = np.nansum(FC_grid)
    return fuel_consumption

    
def get_global_total_FC_yearly(year, emis_data, BA_data):
    FC_grid = get_grid_fuel_consumption(year, 12, emis_data, BA_data)
    fuel_consumption = np.nansum(FC_grid)
    return fuel_consumption


def get_global_mean_FC_monthly(month, emis_data, BA_data, grid_data):
    FC_data = get_grid_fuel_consumption(month/12.,1,emis_data, BA_data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.ma.average(FC_data, weights=grid_data["area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["area"])
    return global_mean_FC 


def plot_global_FC_monthly(no_months, emis_data, BA_data, grid_data):
    x_data = range(len(emis_data["time"])-1)
    y_data = []
    for x in x_data[-no_months:]:
        print("%.2f" % ((x-x_data[-no_months])/
                float(len(x_data[-no_months:]))))
        y_data.append(get_global_mean_FC_monthly(x, emis_data, BA_data, grid_data))
    plt.plot(x_data[-no_months:], y_data, color='r', linewidth=2.0,
               label='JSBACH Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Month')
    plt.legend()
    plt.show()
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x, emis_data, BA_data, grid_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='JSBACH Results')
    plt.ylabel('Fuel Consumption ($kg\, C/m^2 \,burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data):
    month_period = int(year_period*12)
    map_data = get_grid_fuel_consumption(year_start,month_period,emis_data,BA_data)
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    
    lats = emis_data["latitude"]
    lons = emis_data["longitude"]
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, 
        urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Carbon emitted per area burned ($kg\, C / m^2 \, burned $)")
    plt.title("Mean Fuel Consumption, JSBACH")
    plt.show()

def plot_burnt_area_map_period(year_start, year_period, BA_data, grid_data):
    month_period = int(year_period*12)
    burnt_area = np.sum(BA_data["burntArea"][year_start*12:year_start*12+month_period],axis=1)
    map_data = np.multiply(np.sum(burnt_area, axis=0),grid_data["area"])
    # Change units.
    map_data = np.divide(map_data, 1e9)
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    
    lats = BA_data["latitude"]
    lons = BA_data["longitude"]
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, 
        urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Burnt Area (billions of $m^2$)")
    plt.title("Mean Burnt Area, JSBACH")
    plt.show()

#plot_burnt_area_map_period(297,16,BA_JSBACH,grid_JSBACH)
#print plot_global_mean_FC_yearly(20, emis_JSBACH, BA_JSBACH, grid_JSBACH)
#plot_FC_map_period(297, 16 , emis_JSBACH, BA_JSBACH)
#print get_global_mean_FC_yearly(300, emis_JSBACH, BA_JSBACH, grid_JSBACH)


#
# Carbon Emissions Analysis
#


def get_grid_emissions(year_start, month_period, emis_data, grid_data):
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
    emis_rate_data = np.multiply(emis_data["fFirepft"][time:time+month_period], grid_data["area"])
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    emissions_pft = np.sum(emis_per_month, axis = 0)
    emissions = np.sum(emissions_pft, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions

def get_global_emissions_monthly(month, emis_data, grid_data):
    emissions_grid = get_grid_emissions(month/12., 1, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions
   
    
def plot_global_emissions_monthly(no_months, emis_data, grid_data):
    x_data = range(len(emis_data["time"])-1)
    y_data = []
    for x in x_data[-no_months:]:
        print("%.2f" % ((x-x_data[-no_months])/
                float(len(x_data[-no_months:]))))
        y_data.append(get_global_emissions_monthly(x, emis_data, grid_data))
    plt.plot(x_data[-no_months:], y_data, color='r', linewidth=2.0,
               label='JSBACH Results')
    plt.ylabel('Carbon Emitted ($Pg/month$)')
    plt.xlabel('Month')
    plt.legend()
    plt.show()
   

def plot_global_emissions_yearly(no_years, emis_data, grid_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='JSBACH Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_emissions_map_period(year_start, year_period, emis_data, grid_data):
    month_period = int(year_period*12)
    map_data = get_grid_emissions(year_start,month_period,emis_data,grid_data)
    # Convert to billions of kg.
    map_data = np.multiply(map_data, (1./(10**9)))
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    lats = emis_data["latitude"]
    lons = emis_data["longitude"]
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
    cb.set_label("Billions of kg of Carbon")
    plt.title("Mean Emissions, JSBACH")
    plt.show()
        

#plot_emissions_map_period(297,16,emis_JSBACH,grid_JSBACH)
#print get_grid_emissions(300,12,data_JSBACH,grid_data_JSBACH)
       
#print get_global_emissions_monthly(3002, data_JSBACH, grid_data_JSBACH)
#plot_global_emissions_monthly(data_JSBACH, grid_data_JSBACH, 100)
#plot_global_emissions_yearly(data_JSBACH, grid_data_JSBACH,30)
#print get_global_emissions_yearly(312, emis_JSBACH, grid_JSBACH)

