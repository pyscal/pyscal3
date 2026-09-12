#include "system.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <map>
#include <string>
#include <any>
#include "voro++.hh"

using namespace voro;

void get_all_neighbors_voronoi(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const double face_area_exponent)
    {

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj, pos;
    int tnx,tny,tnz, ti, nverts;
    double vol, weightsum;
    // offset added by voro++ to every vertex coordinate; we want vertices
    // relative to the particle, so it must be exactly zero
    const double x = 0.0, y = 0.0, z = 0.0;

    vector<int> neigh,f_vert, vert_nos;
    vector<double> facearea, v, faceperimeters;
    voronoicell_neighbor c;


    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    //vector<bool> mask_1 = atoms[py::str("mask_1")].cast<vector<bool>>();
    //vector<bool> mask_2 = atoms[py::str("mask_2")].cast<vector<bool>>();
    vector<bool> ghost = atoms[py::str("ghost")].cast<vector<bool>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);

    //specific properties related to  Voronoi
    vector<double> volume(nop);
    vector<vector<int>> face_vertices(nop);
    vector<vector<double>> face_perimeters(nop);
    vector<vector<double>> vertex_vectors(nop);
    vector<vector<int>> vertex_numbers(nop);
    vector<vector<vector<double>>> vertex_positions(nop);
    vector<vector<bool>> vertex_unique(nop);

    // ------------------------------------------------------------------
    // Cell geometry for voro++.
    //
    // voro++'s periodic container expects the cell in lower-triangular
    // form, a = (bx,0,0), b = (bxy,by,0), c = (bxz,byz,bz).  For a
    // triclinic cell we build the orthonormal frame Q (rows e1,e2,e3)
    // that brings the cell into that form, rotate the positions into it
    // and rotate the vertex vectors back afterwards.  Distances, face
    // areas and volumes are frame independent.
    // ------------------------------------------------------------------
    double Q[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double bx = box[0], bxy = 0.0, by = box[1], bxz = 0.0, byz = 0.0, bz = box[2];
    if (triclinic == 1){
        // cell vectors are the columns of rot (rot = cell.T)
        double av[3] = {rot[0][0], rot[1][0], rot[2][0]};
        double bv[3] = {rot[0][1], rot[1][1], rot[2][1]};
        double cv[3] = {rot[0][2], rot[1][2], rot[2][2]};
        bx = sqrt(av[0]*av[0] + av[1]*av[1] + av[2]*av[2]);
        for (int k=0; k<3; k++) Q[0][k] = av[k]/bx;
        bxy = bv[0]*Q[0][0] + bv[1]*Q[0][1] + bv[2]*Q[0][2];
        double bp[3];
        for (int k=0; k<3; k++) bp[k] = bv[k] - bxy*Q[0][k];
        by = sqrt(bp[0]*bp[0] + bp[1]*bp[1] + bp[2]*bp[2]);
        for (int k=0; k<3; k++) Q[1][k] = bp[k]/by;
        bxz = cv[0]*Q[0][0] + cv[1]*Q[0][1] + cv[2]*Q[0][2];
        byz = cv[0]*Q[1][0] + cv[1]*Q[1][1] + cv[2]*Q[1][2];
        double cp[3];
        for (int k=0; k<3; k++) cp[k] = cv[k] - bxz*Q[0][k] - byz*Q[1][k];
        bz = sqrt(cp[0]*cp[0] + cp[1]*cp[1] + cp[2]*cp[2]);
        for (int k=0; k<3; k++) Q[2][k] = cp[k]/bz;
    }

    // Process one Voronoi cell computed by voro++ for atom ti
    auto process_cell = [&](voronoicell_neighbor& c, int ti){
        c.face_areas(facearea);
        c.neighbors(neigh);
        c.face_orders(f_vert);
        c.face_vertices(vert_nos);
        c.vertices(x,y,z,v);
        c.face_perimeters(faceperimeters);

        vol = c.volume();

        weightsum = 0.0;
        for (size_t i=0; i<facearea.size(); i++){
            weightsum += pow(facearea[i], face_area_exponent);
        }

        nverts = int(v.size())/3;
        if (triclinic == 1){
            // rotate vertex vectors back into the original Cartesian frame
            for(int si=0; si<nverts; si++){
                double vx = v[3*si], vy = v[3*si+1], vz = v[3*si+2];
                for (int k=0; k<3; k++){
                    v[3*si+k] = Q[0][k]*vx + Q[1][k]*vy + Q[2][k]*vz;
                }
            }
        }

        volume[ti] = vol;
        vertex_vectors[ti] = v;
        vertex_numbers[ti] = vert_nos;

        //clean up and add vertex positions
        pos = positions[ti];
        for(int si=0; si<nverts; si++){
            vector<double> temp;
            for(int k=0; k<3; k++){
                temp.emplace_back(v[3*si+k]+pos[k]);
            }
            vertex_positions[ti].emplace_back(temp);
            vertex_unique[ti].emplace_back(!ghost[ti]);
        }

        for (size_t tj=0; tj<neigh.size(); tj++){
            d = get_abs_distance(positions[ti], positions[neigh[tj]],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            neighbors[ti].emplace_back(neigh[tj]);
            neighbordist[ti].emplace_back(d);
            neighborweight[ti].emplace_back(pow(facearea[tj], face_area_exponent)/weightsum);

            face_vertices[ti].emplace_back(f_vert[tj]);
            face_perimeters[ti].emplace_back(faceperimeters[tj]);

            diffi.clear();
            diffi.emplace_back(diffx);
            diffi.emplace_back(diffy);
            diffi.emplace_back(diffz);

            diff[ti].emplace_back(diffi);

            convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);

            r[ti].emplace_back(tempr);
            phi[ti].emplace_back(tempphi);
            theta[ti].emplace_back(temptheta);
        }
        // per-atom cutoff = distance to the farthest Voronoi neighbour, so
        // that every Voronoi neighbour counts for clustering, the ACE
        // radial cutoff and the local density
        double dmax = 0.0;
        for (size_t tj=0; tj<neighbordist[ti].size(); tj++){
            if (neighbordist[ti][tj] > dmax) dmax = neighbordist[ti][tj];
        }
        cutoff[ti] = dmax;
    };

    if (triclinic == 1){
        // block counts as in pre_container::guess_optimal
        double ilscale = pow(double(nop)/(optimal_particles*bx*by*bz), 1.0/3.0);
        tnx = int(bx*ilscale + 1);
        tny = int(by*ilscale + 1);
        tnz = int(bz*ilscale + 1);
        container_periodic con(bx, bxy, by, bxz, byz, bz, tnx, tny, tnz, 8);
        for(int i=0; i<nop; i++){
            const vector<double>& p = positions[i];
            double px = Q[0][0]*p[0] + Q[0][1]*p[1] + Q[0][2]*p[2];
            double py = Q[1][0]*p[0] + Q[1][1]*p[1] + Q[1][2]*p[2];
            double pz = Q[2][0]*p[0] + Q[2][1]*p[1] + Q[2][2]*p[2];
            // put() remaps the particle into the primary cell
            con.put(i, px, py, pz);
        }
        c_loop_all_periodic cl(con);
        if (cl.start()) do if(con.compute_cell(c,cl)) {
            process_cell(c, cl.pid());
        } while (cl.inc());
    }
    else{
        pre_container pcon(0.00, box[0], 0.00, box[1], 0.0, box[2], true, true, true);
        for(int i=0; i<nop; i++){
            pos = positions[i];
            pos = remap_atom_into_box(pos, triclinic, rot, rotinv, box);
            pcon.put(i, pos[0], pos[1], pos[2]);
        }
        pcon.guess_optimal(tnx, tny, tnz);
        container con(0.00, box[0], 0.00, box[1], 0.0, box[2], tnx, tny, tnz, true, true, true, nop);
        pcon.setup(con);

        c_loop_all cl(con);
        if (cl.start()) do if(con.compute_cell(c,cl)) {
            process_cell(c, cl.pid());
        } while (cl.inc());
    }


    //calculation over lets assign
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;
    atoms[py::str("voronoi_volume")] = volume;
    atoms[py::str("face_vertices")] = face_vertices;
    atoms[py::str("face_perimeters")] = face_perimeters;
    atoms[py::str("vertex_vectors")] = vertex_vectors;
    atoms[py::str("vertex_numbers")] = vertex_numbers;
    atoms[py::str("vertex_is_unique")] = vertex_unique;
    atoms[py::str("vertex_positions")] = vertex_positions;
} 


bool check_if_in_box(const vector<double>& pos,
    const vector<double>& box){
    if ((pos[0] < -0.01) || (pos[0] > box[0]+0.01)) return false;
    else if ((pos[1] < -0.0001) || (pos[1] > box[1])) return false;
    else if ((pos[2] < -0.0001) || (pos[2] > box[2])) return false;
    else return true;
}

vector<vector<double>> clean_voronoi_vertices(py::dict& atoms,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const double distance_cutoff){

    vector<vector<vector<double>>> positions = atoms[py::str("vertex_positions")].cast<vector<vector<vector<double>>>>();
    vector<vector<bool>> vertex_unique = atoms[py::str("vertex_is_unique")].cast<vector<vector<bool>>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    //vector<bool> ghost = atoms[py::str("ghost")].cast<vector<bool>>();
    
    int nop = positions.size();

    double d, diffx, diffy, diffz;
    int nn;

    for(int ti=0; ti<nop; ti++){
        //if (ghost[ti]) continue;
        for(int vi=0; vi<positions[ti].size(); vi++){
            if (!vertex_unique[ti][vi]) continue;
            if (!check_if_in_box(positions[ti][vi], box)){
                vertex_unique[ti][vi] = false;
                continue;
            }
            for(int vj=vi+1; vj<positions[ti].size(); vj++){
                if (!vertex_unique[ti][vj]) continue;
                d = get_abs_distance(positions[ti][vi], positions[ti][vj],
                    triclinic, rot, rotinv, box, 
                    diffx, diffy, diffz);
                if (d < distance_cutoff){
                    vertex_unique[ti][vj] = false;
                }                
            }
            for(int tj=0; tj<neighbors[ti].size(); tj++){
                nn = neighbors[ti][tj];
                if (ti==nn) continue;
                //if (ghost[nn]) continue;
                for(int vj=0; vj<positions[nn].size(); vj++){
                    if (!vertex_unique[nn][vj]) continue;
                    if (!check_if_in_box(positions[nn][vj], box)){
                        vertex_unique[nn][vj] = false;
                        continue;
                    }
                    d = get_abs_distance(positions[ti][vi], positions[nn][vj],
                        triclinic, rot, rotinv, box, 
                        diffx, diffy, diffz);
                    if (d < distance_cutoff){
                        vertex_unique[nn][vj] = false;
                    }                                        
                }
            }
        }
    }

    vector<vector<double>> unique_positions;
    for(int ti=0; ti<nop; ti++){
        for(int tj=0; tj<vertex_unique[ti].size(); tj++){
            if(vertex_unique[ti][tj]){
                unique_positions.emplace_back(positions[ti][tj]);
            }
        }
    }
    
    return unique_positions;    
}	