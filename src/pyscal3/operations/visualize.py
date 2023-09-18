import numpy as np
import itertools
import warnings
import matplotlib.pyplot as plt
import matplotlib

def plot_by_selection(sys, height=500, width=500, size=5, ):

    try:
        import ipyvolume as ipv
    except ImportError:
        print("Install ipyvolume for visualisation")

    cFalse = '#b0bec5'
    cTrue = '#ef5350'

    scatter = ipv.scatter(sys.atoms.positions[:,0][sys.atoms.selection], 
                          sys.atoms.positions[:,1][sys.atoms.selection], 
                          sys.atoms.positions[:,2][sys.atoms.selection], 
                          marker='sphere',
                          size=size,
                          lighting=True,
                          color=cTrue)

    scatter = ipv.scatter(sys.atoms.positions[:,0][np.logical_not(sys.atoms.selection)], 
                          sys.atoms.positions[:,1][np.logical_not(sys.atoms.selection)], 
                          sys.atoms.positions[:,2][np.logical_not(sys.atoms.selection)], 
                          marker='sphere',
                          size=size,
                          lighting=True,
                          color=cFalse)

    #plotting box
    box = sys.box.copy()
    origin = np.array([0,0,0])
    combos = list(itertools.combinations(range(3), 2))
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]        
        unn = np.vstack((f1, f2))
        ipv.plot(unn[:,0], unn[:,1], unn[:,2], color='#546e7a')
    
    ipv.style.box_off()
    ipv.style.axes_off()
    ipv.squarelim()
    
    #done
    return ipv.gcf()

def plot_by_property(sys, colorby, ids=None, 
    indices=None, 
    condition=None, 
    cmap = 'viridis', 
    height=500, width=500, size=5, ):


    try:
        import ipyvolume as ipv
    except ImportError:
        print("Install ipyvolume for visualisation")
    
    colorby = colorby.copy().astype(int)

    sys.apply_selection(ids=ids, indices=indices, condition=condition)
    colorby = [x for count, x in enumerate(colorby) if sys.atoms.selection[count]]

    fig = ipv.figure(debug=False, width=width, height=height)
 
    prop = (colorby.copy() - min(colorby))/(max(colorby) - min(colorby))
    cmap = matplotlib.cm.get_cmap(cmap)
    colors = [cmap(x) for x in prop]

    scatter = ipv.scatter(sys.atoms.positions[:,0][sys.atoms.selection], 
                          sys.atoms.positions[:,1][sys.atoms.selection], 
                          sys.atoms.positions[:,2][sys.atoms.selection], 
                          marker='sphere',
                          size=size,
                          lighting=True,
                          color=colors)

    sys.remove_selection()

    #plotting box
    box = sys.box.copy()
    origin = np.array([0,0,0])
    combos = list(itertools.combinations(range(3), 2))
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]        
        unn = np.vstack((f1, f2))
        ipv.plot(unn[:,0], unn[:,1], unn[:,2], color='#546e7a')
    
    ipv.style.box_off()
    ipv.style.axes_off()
    ipv.squarelim()
    
    #done
    return ipv.gcf()


def plot_simple(sys, colors=None, height=500, width=500, size=5, ):

    try:
        import ipyvolume as ipv
    except ImportError:
        print("Install ipyvolume for visualisation")
    if colors is None:
        colors =[ '#33a02c', '#fb9a99', '#e31a1c',
                 '#a6cee3', '#1f78b4', '#b2df8a', 
                 '#fdbf6f', '#ff7f00', '#cab2d6', 
                 '#6a3d9a', '#ffff99', '#b15928']
    comp_ints = sys.atoms.composition_ints
    if len(comp_ints) > len(colors):
        warnings.warn('less colors than number of species, expect repetitions, or provide more colors')
        diff = np.ceil((len(comp_ints)-len(colors))/len(colors))
        colors = colors*diff
        
    fig = ipv.figure(debug=False, width=width, height=height)
    for key in comp_ints.keys():
        scatter = ipv.scatter(sys.atoms.positions[:,0][sys.atoms.types == key], 
                              sys.atoms.positions[:,1][sys.atoms.types == key], 
                              sys.atoms.positions[:,2][sys.atoms.types == key], 
                              marker='sphere',
                              size=size,
                              lighting=True,
                              color=colors[key-1])

    
    #plotting box
    box = sys.box.copy()
    origin = np.array([0,0,0])
    combos = list(itertools.combinations(range(3), 2))
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]        
        unn = np.vstack((f1, f2))
        ipv.plot(unn[:,0], unn[:,1], unn[:,2], color='#546e7a')
    
    ipv.style.box_off()
    ipv.style.axes_off()
    ipv.squarelim()
    
    #done
    return ipv.gcf()