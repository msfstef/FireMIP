import spatial_comparison as spatial
import temporal_comparison as temporal
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
    print "~Maps Generated~"


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



        

