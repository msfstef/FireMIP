import spatial_comparison as spatial
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def get_var_name(var):
    if var == 'FC':
        title = 'fuel_consumption'
    elif var == 'emis':
        title = 'emissions'
    elif var == 'BA':
        title = 'burnt_area'
    return title

def generate_maps():
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno']
    var_list = ['FC','emis','BA']
    
    for model in model_list:
        for var in var_list:
            fig=spatial.plot_map(1997,15,model,var,
                                binned=True, save=True)
            fig.savefig('./figures/'+model+'/'+
                        get_var_name(var)+'_map.png')
            plt.close(fig)
        print(model+' finished!')
    print "~Maps Generated~"

