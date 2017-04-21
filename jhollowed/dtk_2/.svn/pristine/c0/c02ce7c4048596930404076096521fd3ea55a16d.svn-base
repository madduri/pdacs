
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

def cat_strings(strings):
    string_result = "";
    for st in strings:
        string_result +=st;
    return string_result;


def replace_strings(strings,target,new_val):
    st_result = []
    for st in strings:
        if(st==target):
            st=str(new_val)
        st_result.append(st)
    return st_result
        
def cat_replace_strings(strings,target,new_val):
    return cat_strings(replace_strings(strings,target,new_val))


def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

dtk_figcount =1
dtk_figs_path = ""

def set_fig_path(path):
    global dtk_figs_path
    dtk_figs_path = path
    ensure_dir(path)

def save_figs(path=None,extension='.pdf'):
    global dtk_figcount
    if(path is None):
        path = dtk_figs_path
    ensure_dir(path)
    for i in plt.get_fignums():
        plt.figure(i)
        plt.savefig(""+path+'fig%d' % (dtk_figcount)+extension)
        dtk_figcount+=1

def get_colors(data,cmap,log=False,vmax=None, vmin=None):
    if(vmax == None):
        vmax = np.max(data)
    if(vmin == None):
        vmin = np.min(data)
    if(log):
        cnorm = colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        cnorm = colors.Normalize(vmin=vmin,vmax=vmax)
    smp = cm.ScalarMappable(norm=cnorm,cmap=cmap)
    clrs = smp.to_rgba(data)
    return clrs
