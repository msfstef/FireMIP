import spatial_comparison as spatial
import temporal_comparison as temporal
import field_observations as field_obs
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def get_var_name(var):
    """
    Returns string of full name of variable for the
    given abbreviation.
    
    e.g. for input var='FC', returns 'fuel_consumption'.
    """
    if var == 'FC':
        title = 'fuel_consumption'
    elif var == 'emis':
        title = 'emissions'
    elif var == 'BA':
        title = 'burnt_area'
    return title


def generate_maps():
    """
    Generates the maps for fuel consumption, burnt area,
    and emissions for all models in their respective folders.
    """
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno', 'spitfire']
    var_list = ['FC','emis','BA']
    
    for model in model_list:
        for var in var_list:
            fig=spatial.plot_map(1997,16,model,var,
                                binned=True, save=True)
            fig.savefig('./figures/'+model+'/'+
                        get_var_name(var)+'_map.png')
            plt.close(fig)
        print(model+' finished!')
    print "~Variable Maps Generated~"


def generate_diff_maps():
    """
    Generates the difference maps for fuel consumption, burnt area,
    and emissions for all models vs GFED in their respective folders.
    """
    model_list = ['jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno', 'spitfire']
    var_list = ['FC','emis','BA']
    
    for model in model_list:
        for var in var_list:
            fig=spatial.plot_diff_map(1997,16,model,var,
                                binned=True, save=True)
            fig.savefig('./figures/'+model+'/'+
                        get_var_name(var)+'_diff_map.png')
            plt.close(fig)
        print(model+' finished!')
    print "~Difference Maps Generated~"


def generate_regional_box_plots():
    """
    Generates regional boxplots for fuel consumption, burnt area,
    and emissions for all models in their respective folders.
    """
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno', 'spitfire']
    var_list = ['FC','emis','BA']
    
    for model in model_list:
        for var in var_list:
            fig=spatial.plot_multimodel_box(1997,16,var,
                                model=model,save=True)
            fig.savefig('./figures/'+model+'/'+
                        get_var_name(var)+'_regional_boxplot.png')
            plt.close(fig)
        print(model+' finished!')
    print "~Box Plots Generated~"
    
    
def generate_model_specific_plots():
    """
    Generates all model specific plots.
    """
    generate_maps()
    generate_diff_maps()
    generate_regional_box_plots()
    print '~Model Specific Plots Generated~'


def generate_multimodel_box_plots():
    """
    Generates global boxplots for fuel consumption, burnt area,
    and emissions for comparison of all models.
    """
    var_list = ['FC','emis','BA']
    
    for var in var_list:
        fig=spatial.plot_multimodel_box(1997,16,var,
                            save=True)
        fig.savefig('./figures/spatial_comparison/'+
              'multimodel_'+get_var_name(var)+'_boxplot.png')
        plt.close(fig)
        print('Boxplot of '+get_var_name(var)+' finished!')
    print "~Box Plots Generated~"


def generate_global_temporal_plots():
    """
    Generates temporal plots for all variables, along with
    tables of means and standard deviations, correlations with
    GFED, and correlations with the multimodel mean.
    """
    var_list = ['FC','emis','BA']
    for var in var_list:
        fig=temporal.plot_time_series(1997,16,var,0,
                means=True,corr_gfed=True,corr_multimodel=True,
                save=True)
        fig.savefig('./figures/temporal_comparison/present_'+
                   get_var_name(var)+'_global.png')
        plt.close(fig)
        print 'Plot of ', get_var_name(var), ' generated!'
    print '~Global Temporal Plots Generated~'


def generate_regional_temporal_plots():
    """
    Generates temporal plots for all variables globally and
    for each individual region, all together in a single figure.
    """
    var_list = ['FC','emis','BA']
    for var in var_list:
        fig=temporal.plot_time_series(1997,16,var,
                        all_regions=True,save=True)
        fig.savefig('./figures/temporal_comparison/present_'+
                   get_var_name(var)+'_global.png')
        plt.close(fig)
        print 'Plot of ', get_var_name(var), ' generated!'
    print '~Regional Temporal Plots Generated~'


def generate_standard_deviation_map():
    """
    Generates two maps of the standard deviations in each
    grid cell between all the models to illustrate regions
    of high intermodel variability for each variable. 
    
    One map is in the lowest available resolution, 
    i.e. that of CTEM, and the other is in the highest, 
    i.e. that of GFED (0.5x0.5 degrees).
    """
    var_list = ['FC','emis','BA']
    # GFED resolution.
    for var in var_list:
        fig = spatial.plot_std_map(1997,16,ref_grid='gfed',
                                    save=True)
        fig.savefig('./figures/spatial_comparison/'+
              get_var_name(var)+'_standard_dev_map_HIRES.png')
        plt.close(fig)
        print 'HiRes stdev plot of ', get_var_name(var), ' generated!'
    
    # CTEM resolution. 
    for var in var_list:
        fig = spatial.plot_std_map(1997,16,ref_grid='ctem',
                                    save=True)
        fig.savefig('./figures/spatial_comparison/'+
              get_var_name(var)+'_standard_dev_map_LORES.png')
        plt.close(fig)
        print 'LoRes stdev plot of ', get_var_name(var), ' generated!'
    print '~Standard Dev. Maps Generated~'

def generate_spatial_correlations_table():
    """
    Generates .csv file with table of spatial
    correlations between all the models and
    GFED.
    """
    var_list = ['FC','emis','BA']
    for var in var_list:
        spatial.get_spatial_correlations(1997,16,var)
    print '~Spatial Correlations Table Generated~'

def generate_field_observations_histogram():
    """
    Generates bar chart of mean deviations of model
    outputs of fuel consumption from field observations
    taken from the van Leeuwen et al (2014) paper.
    """
    fig = field_obs.plot_bar_chart(save=True)
    fig.savefig('./figures/spatial_comparison/'+
          'field_observations_deviations.png')
    plt.close(fig)
    print '~Field Observations Bar Chart Generated~'


def generate_multimodel_plots():
    """
    Generates all multimodel plots.
    """
    generate_multimodel_box_plots()
    generate_global_temporal_plots()
    generate_regional_temporal_plots()
    generate_standard_deviation_map()
    generate_spatial_correlations_table()
    generate_field_observations_histogram()
    print '~Multimodel Plots Generated~'

#generate_multimodel_plots()
