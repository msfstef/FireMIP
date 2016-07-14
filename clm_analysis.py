from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_CLM_pft = Dataset('../../model_data/CLM_S1_fFirepft.nc', 
                'r', format = 'NETCDF4')
# Total emissions, not per pft, because burnt area per pft is not available.                   
emis_CLM = Dataset('../../model_data/CLM_S1_CFFIRE.nc', 
                'r', format = 'NETCDF4')                
BA_CLM = Dataset('../../model_data/CLM_S1_BAF.nc', 
                'r', format = 'NETCDF4')                
grid_backup_CLM = Dataset('../../model_data/CLM-gridcell.nc', 
                    'r', format = 'NETCDF4')
grid_CLM = Dataset('../../model_data/CLM-gridarea-nomask.nc', 
                    'r', format = 'NETCDF4')            
# CLM has no separate time data, must use JSBACH data which matches CLM.                  
time_data = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')



#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, time_data):
    time = int(year*12)
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
    days_per_month = np.array(days_per_month)
    
    BA = BA_data["BAF"][year*12:year*12+month_period]
    BA = np.multiply(BA,days_per_month[:,np.newaxis,np.newaxis]/100.)
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data, time_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data, time_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data,time_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CLM Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_burnt_area_map_period(year_start, year_period, BA_data, grid_data, time_data):
    map_data = get_grid_burnt_area(year_start,year_period*12,BA_data,grid_data,time_data)
    # Change units.
    map_data = np.divide(map_data, 1e9)
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    
    lats = grid_data["lat"]
    lons = grid_data["lon"]
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
    plt.title("Mean Burnt Area, CLM")
    plt.show()
 
    
#plot_global_BA_yearly(16,BA_CLM, grid_CLM, time_data)
#plot_burnt_area_map_period(297,16,BA_CLM, grid_CLM, time_data)
#print get_global_BA_yearly(300, BA_CLM, grid_CLM, time_data)


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
    
    emis_rate_data = np.multiply(emis_data["CFFIRE"][time:time+month_period], grid_data["cell_area"])
    emis_per_month = np.multiply(emis_rate_data, sec_per_month[:, np.newaxis, np.newaxis])
    emissions = np.sum(emis_per_month, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data, time_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions

def get_global_emissions_monthly(month, emis_data, grid_data, time_data):
    emissions_grid = get_grid_emissions(month/12., 1, emis_data, grid_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions
   
    
def plot_global_emissions_monthly(no_months, emis_data, grid_data, time_data):
    x_data = range(len(data["time"])-1)
    y_data = []
    for x in x_data[-no_months:]:
        print("%.2f" % ((x-x_data[-no_months])/
                float(len(x_data[-no_months:]))))
        y_data.append(get_global_emissions_monthly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data[-no_months:], y_data, color='r', linewidth=2.0,
               label='CLM Results')
    plt.ylabel('Carbon Emitted ($Pg/month$)')
    plt.xlabel('Month')
    plt.legend()
    plt.show()
   

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
               label='CLM Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, emis_data, grid_data, time_data):
    month_period = int(year_period*12)
    map_data = get_grid_emissions(year_start,month_period, emis_data, grid_data, time_data)
    # Convert to billions of kg.
    map_data = np.multiply(map_data, (1./(10**9)))
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    lats = grid_data["lat"]
    lons = grid_data["lon"]
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Mean Emissions, CLM")
    plt.show()
    

#plot_global_emissions_yearly(20,emis_CLM,grid_CLM,time_data)
#plot_map_period(297,16,emis_CLM,grid_CLM,time_data)    
#print get_global_emissions_yearly(310,emis_CLM,grid_CLM,time_data)



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
    days_per_month = np.array(days_per_month)
    
    # Ignore division by zero warning. Returns inf.
    np.seterr(divide='ignore')
    
    
    BA = np.multiply(BA_data["BAF"][time:time+month_period], days_per_month[:, np.newaxis,np.newaxis])
    BA = np.divide(BA,100)
    inverse_BA = 1./BA
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    # Burnt area data per pft is not available, so I use CFFIRE for this calculation.
    actual_emission_data = emis_data["CFFIRE"][time:time+month_period]
    FC_data = np.multiply(actual_emission_data, inverse_BA)

    FC_per_month = np.multiply(FC_data, sec_per_month[:, np.newaxis, np.newaxis])
    fuel_consumption = np.sum(FC_per_month, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, time_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    

def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, time_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,time_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,time_data)
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
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data,time_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CLM Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data, time_data):
    month_period = int(year_period*12)
    map_data = get_grid_fuel_consumption(year_start,month_period,emis_data, BA_data, time_data)
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    lats = BA_data["lat"]
    lons = BA_data["lon"]
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
    cb.set_label("Carbon emitted per area burned ($kg\, C / m^2 \, burned$)")
    plt.title("Mean Fuel Consumption, CLM")
    plt.show()
    

#print get_global_mean_FC_yearly_rough(300,emis_CLM,BA_CLM,grid_CLM,time_data)

#print get_global_FC_yearly(300, emis_CLM, BA_CLM, time_data)
#plot_FC_map_period(297,16 , emis_CLM, BA_CLM, time_data)
#plot_global_mean_FC_yearly(20, emis_CLM, BA_CLM, grid_CLM,time_data)
#print get_global_mean_FC_yearly(300,emis_CLM,BA_CLM,grid_CLM, time_data)
