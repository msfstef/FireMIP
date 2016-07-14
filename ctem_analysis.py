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

def plot_burnt_area_map_period(year_start, year_period, BA_data, grid_data, landCover_data):
    map_data = get_grid_burnt_area(year_start,year_period*12,BA_data,grid_data,landCover_data)
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
    plt.title("Mean Burnt Area, CTEM")
    plt.show()
 
    
#plot_global_BA_yearly(16,BA_CTEM, grid_CTEM, landCover_CTEM)
#plot_burnt_area_map_period(136,16,BA_CTEM, grid_CTEM, landCover_CTEM) 
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
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    lats = grid_data["lat"]
    lons = grid_data["lon"]
    emis_data = map_data
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, emis_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Mean Emissions, CTEM")
    plt.show()
    

#plot_global_emissions_yearly(16,emis_CTEM,grid_CTEM,landCover_CTEM)
#plot_map_period(136,15.999,emis_CTEM,grid_CTEM, landCover_CTEM)
#print get_global_emissions_yearly(151,emis_CTEM,grid_CTEM, landCover_CTEM)



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
    
    landCover = landCover_data["landCoverFrac"][time:time+month_period]
    
    # Convert from percentage to decimal and eliminate meaningless values.
    BA = np.divide(BA_data["burntArea"][time:time+month_period, :9], 100)
    BA[BA<0.] = 0
    BA[BA>1.] = 0
    #BA = np.multiply(landCover, BA)

    # Add up pft dependency.
    BA=np.sum(BA,axis=1)
    inverse_BA = 1./BA
    
    actual_emission_data= np.multiply(landCover, emis_data["fFirepft"][time:time+month_period, :9])
    # Add up pft dependency.
    actual_emission_data = np.sum(actual_emission_data, axis=1)
    
    FC_data = np.multiply(actual_emission_data, inverse_BA)  
    FC_per_month = np.multiply(FC_data, 
                sec_per_month[:, np.newaxis, np.newaxis])
    fuel_consumption = np.sum(FC_per_month, axis = 0)
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
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data, landCover_data):
    month_period = int(year_period*12)
    map_data = get_grid_fuel_consumption(year_start,month_period,emis_data, BA_data, landCover_data)
    
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
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
    plt.title("Mean Fuel Consumption, CTEM")
    plt.show()
    

#print get_global_mean_FC_yearly_rough(136, emis_CTEM, BA_CTEM, grid_CTEM, landCover_CTEM)

#plot_FC_map_period(136,15.999 , emis_CTEM,BA_CTEM,landCover_CTEM)
#plot_global_mean_FC_yearly(20, emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)
#print get_global_mean_FC_yearly(150,emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)

