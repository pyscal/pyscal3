#include "system.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
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

double get_abs_distance(vector<double> pos1, vector<double> pos2, 
	const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
	const vector<double>& box,
    double& diffx,
    double& diffy,
    double& diffz){
    /*
    Get absolute distance between two atoms
    */

    double abs, ax, ay, az;
    diffx = pos1[0] - pos2[0];
    diffy = pos1[1] - pos2[1];
    diffz = pos1[2] - pos2[2];


    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

        //scale to match the triclinic box size
        diffx = ax*box[0];
        diffy = ay*box[1];
        diffz = az*box[2];

        //now check pbc
        //nearest image
        diffx -= box[0]*round(diffx/box[0]);   // nearest image, any distance
        diffy -= box[1]*round(diffy/box[1]);   // nearest image, any distance
        diffz -= box[2]*round(diffz/box[2]);   // nearest image, any distance

        //now divide by box vals - scale down the size
        diffx = diffx/box[0];
        diffy = diffy/box[1];
        diffz = diffz/box[2];

        //now transform back to normal system
        ax = rot[0][0]*diffx + rot[0][1]*diffy + rot[0][2]*diffz;
        ay = rot[1][0]*diffx + rot[1][1]*diffy + rot[1][2]*diffz;
        az = rot[2][0]*diffx + rot[2][1]*diffy + rot[2][2]*diffz;

        //now assign to diffs and calculate distnace
        diffx = ax;
        diffy = ay;
        diffz = az;

        //finally distance
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

    }
    else{
        //nearest image
        diffx -= box[0]*round(diffx/box[0]);   // nearest image, any distance
        diffy -= box[1]*round(diffy/box[1]);   // nearest image, any distance
        diffz -= box[2]*round(diffz/box[2]);   // nearest image, any distance
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }
    return abs;
}

vector<double> get_distance_vector(vector<double> pos1, 
    vector<double> pos2, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    double diffx, diffy, diffz, dist;

    dist = get_abs_distance(pos1, pos2, triclinic, rot, rotinv, box, diffx, diffy, diffz);

    vector<double> dvec;
    dvec.emplace_back(diffx);
    dvec.emplace_back(diffy);
    dvec.emplace_back(diffz);
    return dvec;
} 

vector<double> remap_atom_into_box(vector<double> pos, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    double dx, dy, dz;
    double abs, ax, ay, az;

    dx = pos[0];
    dy = pos[1];
    dz = pos[2];

    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*dx + rotinv[0][1]*dy + rotinv[0][2]*dz;
        ay = rotinv[1][0]*dx + rotinv[1][1]*dy + rotinv[1][2]*dz;
        az = rotinv[2][0]*dx + rotinv[2][1]*dy + rotinv[2][2]*dz;

        //scale to match the triclinic box size
        dx = ax*box[0];
        dy = ay*box[1];
        dz = az*box[2];

        //now check pbc
        //nearest image
        dx -= box[0]*floor(dx/box[0]);   // wrap into [0, L)
        dy -= box[1]*floor(dy/box[1]);   // wrap into [0, L)
        dz -= box[2]*floor(dz/box[2]);   // wrap into [0, L)

        //now divide by box vals - scale down the size
        dx = dx/box[0];
        dy = dy/box[1];
        dz = dz/box[2];

        //now transform back to normal system
        ax = rot[0][0]*dx + rot[0][1]*dy + rot[0][2]*dz;
        ay = rot[1][0]*dx + rot[1][1]*dy + rot[1][2]*dz;
        az = rot[2][0]*dx + rot[2][1]*dy + rot[2][2]*dz;

        //now assign to diffs and calculate distnace
        dx = ax;
        dy = ay;
        dz = az;
    }
    else{
        //nearest image
        dx -= box[0]*floor(dx/box[0]);   // wrap into [0, L)
        dy -= box[1]*floor(dy/box[1]);   // wrap into [0, L)
        dz -= box[2]*floor(dz/box[2]);   // wrap into [0, L)
    }
    
    vector<double> rpos;
    rpos.emplace_back(dx);
    rpos.emplace_back(dy);
    rpos.emplace_back(dz);
    return rpos;
} 


vector<double> remap_and_displace_atom(vector<double> pos, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    const vector<double>& perturbation){

    double dx, dy, dz;
    double abs, ax, ay, az;

    dx = pos[0];
    dy = pos[1];
    dz = pos[2];

    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*dx + rotinv[0][1]*dy + rotinv[0][2]*dz;
        ay = rotinv[1][0]*dx + rotinv[1][1]*dy + rotinv[1][2]*dz;
        az = rotinv[2][0]*dx + rotinv[2][1]*dy + rotinv[2][2]*dz;

        //scale to match the triclinic box size
        dx = ax*box[0];
        dy = ay*box[1];
        dz = az*box[2];

        //now check pbc
        //nearest image
        dx -= box[0]*floor(dx/box[0]);   // wrap into [0, L)
        dy -= box[1]*floor(dy/box[1]);   // wrap into [0, L)
        dz -= box[2]*floor(dz/box[2]);   // wrap into [0, L)

        //now divide by box vals - scale down the size
        dx = dx/box[0];
        dy = dy/box[1];
        dz = dz/box[2];

        dx += perturbation[0];
        dy += perturbation[1];
        dz += perturbation[2];

        //now transform back to normal system
        ax = rot[0][0]*dx + rot[0][1]*dy + rot[0][2]*dz;
        ay = rot[1][0]*dx + rot[1][1]*dy + rot[1][2]*dz;
        az = rot[2][0]*dx + rot[2][1]*dy + rot[2][2]*dz;

        //now assign to diffs and calculate distnace
        dx = ax;
        dy = ay;
        dz = az;
    }
    else{
        //nearest image
        dx -= box[0]*floor(dx/box[0]);   // wrap into [0, L)
        dy -= box[1]*floor(dy/box[1]);   // wrap into [0, L)
        dz -= box[2]*floor(dz/box[2]);   // wrap into [0, L)

        dx = dx/box[0];
        dy = dy/box[1];
        dz = dz/box[2];

        dx += perturbation[0];
        dy += perturbation[1];
        dz += perturbation[2];

        dx = dx*box[0];
        dy = dy*box[1];
        dz = dz*box[2];

    }
    
    vector<double> rpos;
    rpos.emplace_back(dx);
    rpos.emplace_back(dy);
    rpos.emplace_back(dz);
    return rpos;
} 

void convert_to_spherical_coordinates(double x, 
    double y, 
    double z, 
    double &r, 
    double &phi, 
    double &theta){

    r = sqrt(x*x+y*y+z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
}
