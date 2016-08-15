import numpy as np
import scipy.stats as stats
from scipy import spatial
import matplotlib.pyplot as plt


import spatial_comparison as spt


#FIELD OBSERVATIONS DATA
observ_latlon = [[18.35,-95.05],[36,-79.1],[-33.93,115.46],
[-2.29,-60.09],[34.1,-117.47],[40,-2],[-7.98,-38.32],
[-15.84,-47.95],[-5.35,-49.15],[35.21,-83.48],[-4.3,-49.03],
[-2.34,-60.09],[-12.4,132.5],[-12.3,133],[60.49,-150.98],
[-5.3,-49.15],[-4.3,-49.03],[-2.61,-60.17],[-25.15,31.14],
[-9.2,-60.5],[-12.35,30.21],[-16.6,27.15],[8.56,-67.25],
[62.69,-141.77],[-15.85,-60.52],[36.6,-118.81],[-15.84,-47.95],
[60.45,89.25],[2.54,-61.28],[61.61,-149.04],[19.3,-105.3],
[-2.5,-48.12],[-9.17,-63.18],[19.3,-105.3],[-2.61,-60.17],
[-10.16,-60.81],[-9.11,-63.16],[-14.52,24.49],[63.38,-158.25],
[61.67,-117.65],[-9.52,-56.06],[-0.52,117.01],[17.59,81.55],
[64.4,-145.74],[55.85,-107.67],[-15,23],[17.65,81.75],[-43.22,146.54],
[54.93,-114.17],[42.4,-124.1],[-2.75,114.32],[38.9,-120.67],
[53.92,-105.7],[-12.38,133.55],[46.73,-117.18],[65,-146],
[-35.77,148.03],[65.03,-147.85],[-9.97,-56.34],[-12.43,131.49],
[12.22,-2.7],[19.5,-99.5],[-2.52,113.79],[-33.68,116.25],
[68.58,-149.72],[-34.2,116.34],[-12.53,-54.88],[-33.91,116.16],
[-12.38,133.55],[34.73,-120.57],[40,-98],[-7.9,-72.44],[34.01,-80.72],
[-2.52,113.79],[55.07,-114.03],[48.87,-85.28],[46.87,-83.33],
[61.69,-107.94],[59.31,-111.02],[61.37,-117.63],[48.92,-85.29],
[46.78,-83.33],[46.89,-83.43],[-10.53,31.14],[63,-142],[-4.5,-49.01],
[-9.2,-60.5],[64.87,-147.71],[-15.51,-47.53],[34.82,-94.13],
[3.75,-60.5],[33.56,-81.7],[-22.85,-47.6],[61.6,-117.2],[34.63,-77.4],
[24.73,-81.4],[-2.37,102.68],[38.9,-120.62],[33.94,118.33],
[58.58,98.92],[58.7,98.42],[55.74,-97.91],[55.74,-97.85],[59.4,-113.03],
[54.29,-107.78],[54.05,-105.81],[53.57,-88.62],[64.06,-139.43],
[-37.09,145.08],[34.8,-82.6],[37.5,-122],[34.29,-118.33],
[32.32,-117.15],[67.14,-150.18],[60.43,-149.17],[63.45,-145.12],
[33.33,-117.47],[63.12,-143.59],[63.08,-142.3],[34.8,-82.6],[-36,148],
[64.45,-148.05],[-3.37,-52.62],[63.5,-145.15]]

# Units are in tonnes of C per hectare burned.
observ_FC = [380,2.3,43,77.4,45,1.1,64.2,7.2,51.6,92.7,102.9,
111.4,4.5,5.1,32.9,44,55,106.6,3.5,4,4.2,4.5,5.5,48.8,165.5,
212,8.2,16.9,24,51.4,91,43.1,20.7,23.6,82.2,109.9,190.5,2.2,
36.7,42.7,191.5,120,4,35,42,2.9,7.7,299,42.5,38,109,111.4,26.6,
1.4,11,40.2,47,56.7,72.8,2.4,4.2,17.4,331.7,28,40,58,154.6,53,
4.8,116,2.1,226.2,6.3,500,1,6,9,9,13,17,24.3,23.2,33.7,42.7,89.5,
189,65,31.7,7.5,4.7,2.6,10.5,20,36.2,5.6,13.3,37,107.6,2.9,16.7,
18.1,19,20,22,31,35,38,39,27,5.2,14.5,15,16,25.2,30.4,40.6,50,50.6,
55.6,65,66.6,83,110.3,129.1]
# Convert to kg C per m^2 burned.
observ_FC = np.array(observ_FC)*0.1

observ_year = [1980,1983,1983,1984,1986,1989,1989,1990,1990,1990,
1990,1990,1991,1991,1991,1991,1991,1991,1992,1992,1992,1992,1992,
1992,1992,1992,1993,1993,1993,1993,1993,1994,1995,1995,1995,1995,
1995,1996,1997,1997,1997,1998,1999,1999,1999,2000,2000,2000,2001,
2002,2002,2002,2003,2004,2004,2004,2004,2004,2004,2005,2005,2006,
2006,2007,2007,2007,2007,2008,2009,2009,2010,2010,2011,'>1997',
'1970-2005','1970-2005','1970-2005','1970-2005','1970-2005',
'1970-2005','1973-1983','1975-1981','1978-1982','1980-1981',
'1987-2005','1990-1991','1990-1993','1990-1999','1992-1994',
'1994-1996','1995-1998','1996-1997','1997-1998','1997-2000',
'1997-2010','1998-2000','2001-2002','2001-2002','2001-2005',
'2002-2007','2002-2007','2003-2004','2003-2004','2003-2004',
'2003-2004','2003-2004','2003-2004','2003-2004','2008-2009',
'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A',
'N/A','N/A','N/A','N/A','N/A']


def find_indices(model):
    """
    For the given model, returns the indices which
    correspond to the fuel consumption field observations
    in the order given in the obeserv_FC list. Also returns
    the shape of the model, used to create a grid of the
    observations.
    """
    lons_list, lats_list = spt.get_lons_lats(model)
    lats, lons = np.meshgrid(lats_list,lons_list)
    latlon_data = zip(lats.ravel(), lons.ravel())
    model_grid = spatial.KDTree(latlon_data)
    query = model_grid.query(observ_latlon)
    
    indices = query[1]
    model_shape = (len(lats_list), len(lons_list))
    lat_ind=np.mod(indices,len(lats_list))
    lon_ind= np.divide(indices,len(lats_list))
    model_indices = [(lat,lon) for lat,lon 
                        in zip(lat_ind,lon_ind)]
    return model_indices, model_shape

def get_observ_grid(model):
    """
    For the given model, creates a grid of the
    same shape as the model containing all the
    FC field observations. 
    
    If two or more field observations correspond
    to the same grid, they are averaged. Although
    this is entirely correct for the sake of comparison
    it is a reasonable way to go about it.
    """
    indices, shape = find_indices(model)
    grid = np.zeros(shape,dtype=object)
    for i in range(len(indices)):
        index = indices[i]
        value = observ_FC[i]
        if grid[index] == 0:
            grid[index] = [value]
        else:
            grid[index].append(value)
    for index,FC_list in np.ndenumerate(grid):
        grid[index] = np.mean(FC_list)
    grid = grid.astype('float64')
    return grid
    

def compare_points(model, region=0):
    """
    Creates a grid of the relative differences between
    the given model's output and field observations. Can
    give integer number to region from 1 to 12 for regional
    comparison. For more details on the regions, check
    spatial_comparison module. 
    """
    observ_grid = get_observ_grid(model)
    observ_grid = spt.get_regional_var_grid(1970,43,
                                    region,model,'FC',
                                    grid = observ_grid)
    model_grid = spt.get_regional_var_grid(1970,43,
                                    region,model,'FC')
                                    
    diff_grid = model_grid-observ_grid
    # Avoid division by 0.
    observ_grid[observ_grid==0]=np.inf
    diff_grid = np.divide(diff_grid,observ_grid)
    diff_grid = np.multiply(diff_grid, 100)
    return diff_grid


def calc_mean_dev_points(model, region=0, no_vals=False):
    """
    For given model and region, calculates the
    mean relative deviation from a grid of relative
    deviations, along with its associated statistical
    error.
    
    Argument no_vals can be set to True to return the
    number of observations in the given region. Used in
    plots.
    """
    diff_grid = compare_points(model, region)
    diff_vals = diff_grid[diff_grid!=0]
    mean_dev = np.mean(diff_vals)
    no_obs = len(diff_vals)
    stderr_dev = np.std(diff_vals)/np.sqrt(no_obs)
    if no_vals:
        return (mean_dev, stderr_dev), no_obs
    else:
        return (mean_dev, stderr_dev)


def calc_mean_dev_total(model, region=0):
    """
    For a given model and region, calculates the mean
    relative deviation between the model and field
    observations by comparing the mean of all the
    observations in the region and all grid cells
    of the model's output in the region.
    
    This method is not preferred, as the field observations
    are extremely sparse compared to the model outputs. 
    """
    observ_grid = get_observ_grid(model)
    observ_grid = spt.get_regional_var_grid(1970,43,
                                    region,model,'FC',
                                    grid = observ_grid)
    model_grid = spt.get_regional_var_grid(1970,43,
                                    region,model,'FC')
    
    observ_vals = observ_grid[observ_grid>0]
    model_vals = model_grid[model_grid>0]
    no_obs = len(observ_vals)
    no_model = len(model_vals)
    observ_mean = np.mean(observ_vals)
    observ_stderr = np.std(observ_vals)/np.sqrt(no_obs)
    model_mean = np.mean(model_vals)
    model_stderr = np.std(model_vals)/np.sqrt(no_model)
    
    mean_dev = ((model_mean-observ_mean)/observ_mean)*100
    stderr_dev = ((((model_stderr+observ_stderr)/
              (model_mean+observ_mean)) +
              (observ_stderr/observ_mean))*mean_dev)
    return (mean_dev, stderr_dev)


def plot_bar_chart(model='all', save=False):
    """
    Plots a bar chart of relative deviations for each
    region for the given model/s. To plot all models
    together, leave the argument model as 'all', which
    is the default setting.
    
    Argument save is set to False by default. Used by the
    generate_figures module.
    """
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno','spitfire',
                 'mc2','globfirm']
    colour_map=iter(plt.cm.Dark2(np.linspace(0,1,len(model_list))))
    region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                   'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    no_obs_list = []
    
    fig, ax = plt.subplots(figsize=(14,10))
    
    error_config = {'ecolor': '0.3'}
    no_of_bars = len(region_names)
    ind = np.arange(no_of_bars)*2
    bar_width = 0.15
    center = 0.7
    if model != 'all':
        model_list=[model]
        bar_width = 1.
        center = 0.5
     
    highest=[]
    for i in range(len(model_list)):
        # Progress Bar
        print("%.2f" % (float(i)/(len(model_list)-1)))
        model_name = model_list[i]
        if i==0:
            colour = 'k'          
        else:
            colour = next(colour_map)
        means = []
        stderrs = []
        for j in range(13):
            if i==0:
                dev, no_obs = calc_mean_dev_points(model_name,
                                       region=j, no_vals=True)
                no_obs_list.append(no_obs)
            else:
                dev = calc_mean_dev_points(model_name,
                                            region=j)
            
            means.append(dev[0])
            stderrs.append(dev[1])
        
        ax.bar(ind+i*bar_width,means,bar_width,
             color=colour,yerr=stderrs, 
             error_kw=error_config,label=model_name.upper())
        
        highest.append(np.max(means))
    
    
    xticks = [region + '\n' + str(no_obs) for region,no_obs in
               zip(region_names,no_obs_list)]
    
    ax.set_title('Mean Regional Model Deviations'+
                ' from Field Observations of FC')
    ax.set_ylabel('Mean Relative Deviation (%)')
    ax.set_xticks(ind+center)
    ax.set_xticklabels(xticks)
    ax.set_xlim([np.min(ind)-.5,np.max(ind)+bar_width*len(model_list)+.5])
    ax.set_ylim([-100.,np.max(highest)*1.1])
    ax.legend(ncol=2)
    plt.tight_layout()
    
    if save:
        return fig
    else:
        plt.show()


#plot_bar_chart()
#print calc_mean_dev_points('orchidee',11,True)
#print calc_mean_dev_total('clm',12)

