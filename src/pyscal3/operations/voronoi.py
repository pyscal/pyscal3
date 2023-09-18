import numpy as np

def calculate_vorovector(system, edge_cutoff=0.05, area_cutoff=0.01, edge_length=False):
    """
    get the voronoi structure identification vector.

    Parameters
    ----------

    edge_cutoff : float, optional
        cutoff for edge length. Default 0.05.


    area_cutoff : float, optional
        cutoff for face area. Default 0.01.

    edge_length : bool, optional
        if True, a list of unrefined edge lengths are returned. Default false.

    Returns
    -------
    vorovector : array like, int
        array of the form (n3, n4, n5, n6)

    Notes
    -----
    Returns a vector of the form `(n3, n4, n5, n6)`, where `n3` is the number
    of faces with 3 vertices, `n4` is the number of faces with 4
    vertices and so on. This can be used to identify structures [1] [2].

    The keywords `edge_cutoff` and `area_cutoff` can be used to tune the values to minimise
    the effect of thermal distortions. Edges are only considered in the analysis if the
    `edge_length/sum(edge_lengths)` is at least `edge_cutoff`. Similarly, faces are only
    considered in the analysis if the  `face_area/sum(face_areas)` is at least `face_cutoff`.

    References
    ----------
    .. [1] Finney, JL, Proc. Royal Soc. Lond. A 319, 1970
    .. [2] Tanemura, M, Hiwatari, Y, Matsuda, H,Ogawa, T, Ogita, N, Ueda, A. Prog. Theor. Phys. 58, 1977

    """
    complete_edge_lengths = []
    vorovectors = []

    for x in range(len(system.atoms['positions'])):
        st = 1
        refined_edges = []
        edge_lengths = []
        for vno in system.atoms['face_vertices'][x]:
            vphase = system.atoms['vertex_numbers'][x][st:st+vno]
            edgecount = 0
            dummy_edge_lengths = []
            #now calculate the length f each edge
            for i in range(-1, len(vphase)-1):
                #get pairs of indices
                #verts are i, i+1
                ipos = system.atoms['vertex_vectors'][x][vphase[i]*3:vphase[i]*3+3]
                jpos = system.atoms['vertex_vectors'][x][vphase[i+1]*3:vphase[i+1]*3+3]

                #now calculate edge length
                edgeln = np.sqrt((ipos[0]-jpos[0])**2 + (ipos[1]-jpos[1])**2 + (ipos[2]-jpos[2])**2)
                dummy_edge_lengths.append(edgeln)

            edge_lengths.append(dummy_edge_lengths)
            st += (vno+1)

        #now all the edge lengths are saved
        for c, ed in enumerate(edge_lengths):
            #normalise the edge lengths
            norm = (ed/np.sum(ed))
            #apply face area cutoff
            if (system.atoms['neighborweight'][x][c] > area_cutoff):
                #check for edge length cutoff
                edgecount = len([cc for cc,x in enumerate(norm) if x > edge_cutoff])
                refined_edges.append(edgecount)

        #now loop over refined edges and collect n3, n4, n5, n6
        vorovector = [0, 0, 0, 0]

        for ed in refined_edges:
            if ed == 3:
                vorovector[0] += 1
            elif ed == 4:
                vorovector[1] += 1
            elif ed == 5:
                vorovector[2] += 1
            elif ed == 6:
                vorovector[3] += 1

        complete_edge_lengths.append(edge_lengths)
        vorovectors.append(vorovector)
    
    system.atoms["edge_lengths"] = np.array(complete_edge_lengths)
    system.atoms["vorovector"] = np.array(vorovectors)

    mapdict = {}
    mapdict["voronoi"] ={}
    mapdict["voronoi"]["vector"] = "vorovector"
    mapdict["voronoi"]["face"] = {}
    mapdict["voronoi"]["face"]["edge_lengths"] = "edge_lengths"
    system.atoms._add_attribute(mapdict)

