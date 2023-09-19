import numpy as np
import itertools
import warnings
import matplotlib.pyplot as plt
import matplotlib

def create_box_plot(box, origin=[0,0,0]):
    """
    Create traces which correspond to the simulation cell

    Parameters
    ----------
    box : list
        dimensions of the simulation box

    origin : list, optional
        Origin of the simulation box. Default [0, 0, 0]

    Returns
    -------
    traces : list of Scatter3d objects

    """

    try:
        from plotly import graph_objs as go
        import ipywidgets as widgets
    except ImportError:
        print("Install ipyvolume for visualisation")

    box = np.array(box)
    origin = np.array(origin)
    combos = list(itertools.combinations(range(3), 2))
    faces = []
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]
        faces.append(np.array(f1))
        faces.append(np.array(f2))
    traces = []
    for face in faces:
        trace = go.Scatter3d(
            x=face[:,0],
            y=face[:,1],
            z=face[:,2],
            mode='lines',
            name='lines',
            line=dict(width=2.0, color='#263238'),
            showlegend=False
        )
        traces.append(trace)
    return traces

def plot_by_selection(sys, radius=10, 
    opacity=1.0):

    try:
        from plotly import graph_objs as go
    except ImportError:
        print("Install plotly for visualisation")

    cFalse = '#b0bec5'
    cTrue = '#ef5350'

    #plotting box
    box = sys.box.copy()
    origin = np.array([0,0,0])
    traces = create_box_plot(box, origin)

    data=go.Scatter3d(
        x=sys.atoms.positions[:,0][sys.atoms.selection],
        y=sys.atoms.positions[:,1][sys.atoms.selection],
        z=sys.atoms.positions[:,2][sys.atoms.selection],
        mode='markers',
        opacity=1.0,
        marker=dict(
            sizemode='diameter',
            sizeref=750,
            size=radius,
            color = cTrue,
            opacity = opacity,
            #colorscale = colorscale,
            #colorbar=dict(thickness=20, title=cmap_title),
            line=dict(width=0.5, color='#455A64'),            
        )
    )

    traces.append(data)

    data=go.Scatter3d(
        x=sys.atoms.positions[:,0][np.logical_not(sys.atoms.selection)],
        y=sys.atoms.positions[:,1][np.logical_not(sys.atoms.selection)],
        z=sys.atoms.positions[:,2][np.logical_not(sys.atoms.selection)],
        mode='markers',
        opacity=1.0,
        marker=dict(
            sizemode='diameter',
            sizeref=750,
            size=radius,
            color = cFalse,
            opacity = opacity,
            #colorscale = colorscale,
            #colorbar=dict(thickness=20),
            line=dict(width=0.5, color='#455A64')
        )
    )

    traces.append(data)

    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig['layout'].update(scene=dict(aspectmode="data"))
    return fig.show()


def plot_by_property(sys, colorby, 
    cmap = 'viridis', 
    radius=10, 
    opacity=1.0,
    hide_zero=False):


    try:
        from plotly import graph_objs as go
    except ImportError:
        print("Install plotly for visualisation")
    
    #colorby = colorby.copy().astype(int)

    #sys.apply_selection(ids=ids, indices=indices, condition=condition)
    #colorby = [x for count, x in enumerate(colorby) if sys.atoms.selection[count]]
    
    #colorby = colorby.copy()
    #a = min(colorby)
    #b = max(colorby)
    #prop = (colorby - a)/(b - a)

    box = sys.box.copy()
    origin = np.array([0,0,0])
    traces = create_box_plot(box, origin)

    if hide_zero:
        data=go.Scatter3d(
            x=sys.atoms.positions[:,0][np.where(colorby > 0)],
            y=sys.atoms.positions[:,1][np.where(colorby > 0)],
            z=sys.atoms.positions[:,2][np.where(colorby > 0)],
            mode='markers',
            opacity=1.0,
            marker=dict(
                sizemode='diameter',
                sizeref=750,
                size=radius,
                color = colorby,
                opacity = opacity,
                colorscale = cmap,
                colorbar=dict(thickness=20),
                line=dict(width=0.5, color='#455A64'),            
            )
        )
    else:
        data=go.Scatter3d(
            x=sys.atoms.positions[:,0],
            y=sys.atoms.positions[:,1],
            z=sys.atoms.positions[:,2],
            mode='markers',
            opacity=1.0,
            marker=dict(
                sizemode='diameter',
                sizeref=750,
                size=radius,
                color = colorby,
                opacity = opacity,
                colorscale = cmap,
                colorbar=dict(thickness=20),
                line=dict(width=0.5, color='#455A64'),            
            )
        )

    traces.append(data)

    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig['layout'].update(scene=dict(aspectmode="data"))
    #add plot
    #sys.remove_selection()
    return fig.show()


def plot_by_boolean(sys, colorby, 
    color = '#ff7f00', 
    radius=10, 
    opacity=1.0,
    hide_zero=False):


    try:
        from plotly import graph_objs as go
    except ImportError:
        print("Install plotly for visualisation")
    
    #colorby = colorby.copy().astype(int)

    #sys.apply_selection(ids=ids, indices=indices, condition=condition)
    #colorby = [x for count, x in enumerate(colorby) if sys.atoms.selection[count]]
    
    #colorby = colorby.copy()
    #a = min(colorby)
    #b = max(colorby)
    #prop = (colorby - a)/(b - a)

    box = sys.box.copy()
    origin = np.array([0,0,0])
    traces = create_box_plot(box, origin)

    data=go.Scatter3d(
        x=sys.atoms.positions[:,0][np.where(colorby == True)],
        y=sys.atoms.positions[:,1][np.where(colorby == True)],
        z=sys.atoms.positions[:,2][np.where(colorby == True)],
        mode='markers',
        opacity=1.0,
        marker=dict(
            sizemode='diameter',
            sizeref=750,
            size=radius,
            color = color,
            opacity = opacity,
            #colorscale = cmap,
            #colorbar=dict(thickness=20),
            line=dict(width=0.5, color='#455A64'),            
        )
    )

    traces.append(data)

    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig['layout'].update(scene=dict(aspectmode="data"))
    #add plot
    #sys.remove_selection()
    return fig.show()



def plot_simple(sys, colors=None,  
    radius=10, 
    opacity=1.0 ):

    try:
        from plotly import graph_objs as go
    except ImportError:
        print("Install plotly for visualisation")
    
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
    

    box = sys.box.copy()
    origin = np.array([0,0,0])
    traces = create_box_plot(box, origin)

    for key in comp_ints.keys():
        data=go.Scatter3d(
            x=sys.atoms.positions[:,0][sys.atoms.types == key],
            y=sys.atoms.positions[:,1][sys.atoms.types == key],
            z=sys.atoms.positions[:,2][sys.atoms.types == key],
            mode='markers',
            opacity=1.0,
            marker=dict(
                sizemode='diameter',
                sizeref=750,
                size=radius,
                color = colors[key-1],
                opacity = opacity,
                #colorbar=dict(thickness=20, title=cmap_title),
                line=dict(width=0.5, color='#455A64'),            
            )
        )

        traces.append(data)

    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig['layout'].update(scene=dict(aspectmode="data"))
    return fig.show()

    
