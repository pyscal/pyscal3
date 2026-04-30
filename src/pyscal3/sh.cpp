#include <iostream>
#include <cmath>
#include <vector>
#include "system.h"

void calculate_factors(const int lm, 
	vector<vector<double>>& alm,
	vector<vector<double>>& blm, 
	vector<vector<double>>& clm,
	vector<double>& dl,
	vector<double>& el)
	{
	
	double t1, t2, t3;
	double temp1, temp2, temp3, temp4;

	//resize and set everything to zero
	alm.resize(lm+1);
	blm.resize(lm+1);
	clm.resize(lm+1);
	dl.resize(lm+1);
	el.resize(lm+1);

	for(int l=0; l<=lm; l++){
		for(int m=0; m<=l; m++){
			alm[l].emplace_back(0);
			blm[l].emplace_back(0);
			clm[l].emplace_back(0);
		}
		dl[l] = 0;
		el[l] = 0;
	}


	for(int l=0; l<=lm; l++){
		
		temp1 = 4*l*l-1;
		temp2 = l*l-2*l+1;
		temp3 = 2*l+1;
		temp4 = 2*l-1;

		for(int m=0; m<l-1; m++){
			
			t1 = sqrt((temp1)/(l*l-m*m));
			t2 = -sqrt((temp2-m*m)/(4*temp2-1));
			t3 = sqrt(((l-m)*(l+m)*temp3)/(temp4));

			alm[l][m] = t1;
			blm[l][m] = t2;
			clm[l][m] = t3;

		}
	}

	for(int l=2; l<=lm; l++){
		dl[l] = sqrt(2*(l-1)+3);
		el[l] = sqrt(1+0.5/double(l));
	}
}


vector<vector<double>> calculate_plm(const int lm,
	const double costheta, 
	const double sintheta)
	{

	vector<vector<double>> alm;
	vector<vector<double>> blm; 
	vector<vector<double>> clm;
	vector<double> dl;
	vector<double> el;

	calculate_factors(lm, alm, blm, clm, dl, el);

	double a1 = sqrt(0.5/3.141592653589793);
	double a2 = sqrt(3.0);
	double a3 = sqrt(1.5);

	vector<vector<double>> plm;
	plm.resize(lm+1);
	for(int l=0; l<=lm; l++){
		for(int m=0; m<=l; m++){
			plm[l].emplace_back(0);
		}
	}


	plm[0][0] = a1;
	plm[1][0] = a1*costheta*a2;
	a1 = a1*sintheta*-a3;
	plm[1][1] = a1;


	for(int l=2; l<=lm; l++){

		int m = 0;
		for(int m=0; m<l-1; m++){
			plm[l][m] = alm[l][m]*(costheta*plm[l-1][m]+blm[l][m]*plm[l-2][m]);
		}

		plm[l][l-1] = dl[l]*costheta*a1;
		a1 *= -el[l]*sintheta;
		plm[l][l] = a1;
	}

	return plm;
}


double dfactorial(int l,
	int m){

    double fac = 1.00;
    for(int i=0;i<2*m;i++){
        fac*=double(l+m-i);
    }
    return (1.00/fac);
}

vector<vector<vector<double>>> calculate_ylm(const int lm,
	const double costheta,
	const double sintheta,
	const double cosphi,
	const double sinphi)
	{
	vector<vector<vector<double>>> ylm;
	vector<vector<double>> plm;

	ylm.resize(lm+1);
	
	for(int l=0; l<=lm; l++){
		ylm[l].resize(2*l+1);
		for(int m=0; m<(2*l+1); m++){
			ylm[l][m].emplace_back(0.0);
			ylm[l][m].emplace_back(0.0);
		}
	}
	plm = calculate_plm(lm, costheta, sintheta);
	double zv = 1.0/sqrt(2.0);
	double cosi = 1.0;
	double cosf = cosphi;
	double sini = 0.0;
	double sinf = -sinphi; 
	double ci, si;
	double fa, fb, factor;
	for(int l=0; l<=lm; l++){

		ylm[l][l][0] = plm[l][0]*zv;

	}
	for(int m=1; m<=lm; m++){

		ci = 2*cosphi*cosi-cosf;
		si = 2*cosphi*sini-sinf;
		sinf = sini;
		sini = si;
		cosf = cosi;
		cosi = ci;

		for(int l=m; l<=lm; l++){
			
			fa = plm[l][m]*ci*zv;
			fb = plm[l][m]*si*zv;
			//factor = sqrt(((2.0*double(l) + 1.0)/ (2.0*PI))*dfactorial(l,m));
			factor = 1.00;

			ylm[l][l+m][0] = factor*fa;
			ylm[l][l+m][1] = factor*fb;
			ylm[l][l-m][0] = factor*fa*pow(-1.0,-m);
			ylm[l][l-m][1] = factor*fb*pow(-1.0,-m);
		}
	}

	return ylm;

}

vector<vector<vector<vector<double>>>> calculate_q_atom(const int lm,
	const vector<double>& theta,
	const vector<double>& phi)
{
	vector<vector<vector<vector<double>>>> ylm_atom;
	//this is like a loop over neighbors
	for(int i=0; i<theta.size(); i++){
		ylm_atom.emplace_back(calculate_ylm(lm, cos(theta[i]), sin(theta[i]), cos(phi[i]), sin(phi[i])));
	}
	return ylm_atom;
}

void calculate_q(py::dict& atoms,
	const int lm)
{
	//we need theta and pi
    vector<vector<double>> theta = atoms[py::str("theta")].cast<vector<vector<double>>>();
    vector<vector<double>> phi = atoms[py::str("phi")].cast<vector<vector<double>>>();
    vector<vector<double>> weights = atoms[py::str("neighborweight")].cast<vector<vector<double>>>();
    int nop = theta.size();
    vector<vector<vector<double>>> qlm_real(nop);
    vector<vector<vector<double>>> qlm_img(nop);
    vector<vector<double>> q(nop);

    int nn;
    double summ, weightsum;
    double realti, imgti;
	for (int ti=0; ti<nop; ti++){
		//calculate ylm first
		auto ylm_atom = calculate_q_atom(lm, theta[ti], phi[ti]);	
		qlm_real[ti].resize(lm+1);
		qlm_img[ti].resize(lm+1);
		
		for(int l=0; l<=lm; l++){
			summ = 0;
			for(int m=0; m<(2*l+1); m++){
				realti = 0;
				imgti = 0;
				weightsum = 0;
				for(int ci=0; ci<theta[ti].size(); ci++){
					//TODO: add condition
					realti += weights[ti][ci]*ylm_atom[ci][l][m][0];
					imgti += weights[ti][ci]*ylm_atom[ci][l][m][1];					
					weightsum += weights[ti][ci];
				}
				//TODO: turn off for Voronoi
				realti = realti/float(weightsum);
				imgti = imgti/float(weightsum);

				qlm_real[ti][l].emplace_back(realti);
				qlm_img[ti][l].emplace_back(imgti);

				summ += realti*realti + imgti*imgti;
			}
			summ = pow(((4.0*PI/(2*l+1)) * summ),0.5);
			q[ti].emplace_back(summ);
		}
	}

	string key1, key2, key3;
	vector<double> qtemp;
	vector<vector<double>> qtemp1, qtemp2;

	for(int l=0; l<=lm; l++){
		key1 = "q"+to_string(l);
		key2 = "q"+to_string(l)+"_real";
		key3 = "q"+to_string(l)+"_imag";

		qtemp.clear();
		qtemp1.clear();
		qtemp2.clear();

		for (int ti=0; ti<nop; ti++){
			qtemp.emplace_back(q[ti][l]);
			qtemp1.emplace_back(qlm_real[ti][l]);
			qtemp2.emplace_back(qlm_img[ti][l]);
		}
	    atoms[py::str(key1)] = qtemp;
	    atoms[py::str(key2)] = qtemp1;
	    atoms[py::str(key3)] = qtemp2;
	}
}

/**********************************************************************
New set of functions that use the old algorithm
**********************************************************************/
double plm(const int l,
	const int m,
	const double theta){

	double x = cos(theta);
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    pll = 0.0;
    if (m < 0 || m > l || fabs(x) > 1.0)
        cerr << "impossible combination of l and m" << "\n";
    pmm=1.0;
    if (m > 0){
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
        for (i=1;i<=m;i++){
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }

    if (l == m)
        return pmm;
    else{
        pmmp1=x*(2*m+1)*pmm;
        if (l == (m+1))
            return pmmp1;
        else{
            for (ll=m+2;ll<=l;ll++){
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
            pmm=pmmp1;
            pmmp1=pll;
            }
        return pll;
        }
    }	
}

double sph_legendre(const int l,
	const int m,
	const double theta){

	double factor = ((2.0*double(l) + 1.0)/ (4.0*PI))*dfactorial(l,m);
	double m_plm = plm(l, m, theta);
	return sqrt(factor)*m_plm;

}

void calculate_qlm(const int l, 
	const int m, 
	const double theta, 
	const double phi, 
	double &ylm_real, 
	double &ylm_imag){

    double factor;
    double m_plm;

    m_plm = sph_legendre(l, abs(m), theta);
    //factor = ((2.0*double(l) + 1.0)/(4.0*PI))*dfactorial(l,m);
    //factor = (1.0/dfactorial(l,m));
    ylm_real = m_plm*cos(double(m)*phi);
    ylm_imag  = m_plm*sin(double(m)*phi);
}


void calculate_q_single(py::dict& atoms,
	const int lm){

	//we need theta and pi
    vector<vector<double>> theta = atoms[py::str("theta")].cast<vector<vector<double>>>();
    vector<vector<double>> phi = atoms[py::str("phi")].cast<vector<vector<double>>>();
    vector<vector<double>> weights = atoms[py::str("neighborweight")].cast<vector<vector<double>>>();
    
    int nop = theta.size();
    vector<vector<double>> qlm_real(nop);
    vector<vector<double>> qlm_img(nop);
    vector<double> q;

    int nn;
    double summ, weightsum;
    double realti, imgti;
    double realylm, imgylm;
	
    for (int ti=0; ti<nop; ti++){
		summ = 0;
		for (int mi=-lm; mi<lm+1; mi++){
            realti = 0.0;
            imgti = 0.0;
            weightsum = 0;
			for(int ci=0; ci<theta[ti].size(); ci++){
				calculate_qlm(lm, mi, theta[ti][ci], phi[ti][ci], realylm, imgylm);
				realti += weights[ti][ci]*realylm;
				imgti += weights[ti][ci]*imgylm;
				weightsum += weights[ti][ci];
			}
			realti = realti/float(weightsum);
			imgti = imgti/float(weightsum);

			qlm_real[ti].emplace_back(realti);
			qlm_img[ti].emplace_back(imgti);

			summ += realti*realti + imgti*imgti;
		}
		//cout<<summ<<endl;
		summ = pow(((4.0*PI/(2*lm+1))*summ),0.5);

		//cout<<summ<<endl;

		q.emplace_back(summ);
    }

	string key1, key2, key3;

	key1 = "q"+to_string(lm);
	key2 = "q"+to_string(lm)+"_real";
	key3 = "q"+to_string(lm)+"_imag";

    atoms[py::str(key1)] = q;
    atoms[py::str(key2)] = qlm_real;
    atoms[py::str(key3)] = qlm_img;
	
}

void calculate_aq_single(py::dict& atoms,
	const int lm){

    double realti, imgti;
    double summ;
    int nns;

    string key1, key2;
	key1 = "q"+to_string(lm)+"_real";
	key2 = "q"+to_string(lm)+"_imag";

    vector<vector<double>> q_real = atoms[py::str(key1)].cast<vector<vector<double>>>();
    vector<vector<double>> q_imag = atoms[py::str(key2)].cast<vector<vector<double>>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<double> q;
    int nop = neighbors.size();
 
    for (int ti= 0;ti<nop;ti++){
        summ = 0;
        for (int mi = 0;mi<2*lm+1;mi++){
            realti = q_real[ti][mi];
            imgti = q_imag[ti][mi];
            
            nns = 0;
            for (int ci = 0;ci<neighbors[ti].size();ci++){
            	realti += q_real[neighbors[ti][ci]][mi];
            	imgti += q_imag[neighbors[ti][ci]][mi];
                nns += 1;
        }
        
        realti = realti/(double(nns+1));
        imgti = imgti/(double(nns+1));

        summ+= realti*realti + imgti*imgti;
        
        }
        
        //normalise summ
        summ = pow(((4.0*PI/(2*lm+1)) * summ),0.5);
        q.emplace_back(summ);
    }
    key1 = "avg_q"+to_string(lm);
    atoms[py::str(key1)] = q;
}


/*-----------------------------------------------------
    Wigner 3j symbol and W_l parameter
    Ref: Steinhardt, Nelson & Ronchetti, Phys. Rev. B 28, 784 (1983)
    Averaged: Lechner & Dellago, J. Chem. Phys. 129, 114707 (2008)
-----------------------------------------------------*/

// Factorial for integers up to ~25 (sufficient for l <= 12)
static double factorial_int(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; i++) result *= i;
    return result;
}

// Wigner 3j symbol via the Racah formula (j1 = j2 = j3 = l)
//
//   ( l   l   l  )
//   ( m1  m2  m3 )
//
// = (-1)^{l - m3} * sqrt( Delta(l,l,l) * (l+m1)!(l-m1)!...(l-m3)! )
//   * sum_t (-1)^t / [t! (l-t)! (l-m1-t)! (l+m2-t)! (m1+t)! (-m2+t)!]
//
// where Delta(l,l,l) = (l!)^3 / (3l+1)!
static double wigner3j(int l, int m1, int m2, int m3) {
    // Selection rules
    if (m1 + m2 + m3 != 0) return 0.0;
    if (abs(m1) > l || abs(m2) > l || abs(m3) > l) return 0.0;

    int J = 3 * l;
    if (J % 2 != 0) return 0.0;  // 3j vanishes when 3l is odd (j1=j2=j3=l)

    // Phase: (-1)^{j1 - j2 - m3} = (-1)^{-m3} = (-1)^{|m3|}
    double phase = (abs(m3) % 2 == 0) ? 1.0 : -1.0;

    // Triangle coefficient  Delta(l,l,l) = (l!)^3 / (3l+1)!
    double delta = factorial_int(l) * factorial_int(l) * factorial_int(l)
                   / factorial_int(J + 1);

    // m-dependent factor: product of (l +/- mi)! for i = 1,2,3
    double mfact = factorial_int(l + m1) * factorial_int(l - m1)
                 * factorial_int(l + m2) * factorial_int(l - m2)
                 * factorial_int(l + m3) * factorial_int(l - m3);

    // Racah sum bounds (all factorial arguments must be >= 0)
    // Denominator factorials: t!, (l-t)!, (l-m1-t)!, (l+m2-t)!, (m1+t)!, (-m2+t)!
    int t_min = 0;
    t_min = max(t_min, -m1);    // m1 + t >= 0
    t_min = max(t_min, m2);     // -m2 + t >= 0

    int t_max = l;               // l - t >= 0
    t_max = min(t_max, l - m1);  // l - m1 - t >= 0
    t_max = min(t_max, l + m2);  // l + m2 - t >= 0

    if (t_min > t_max) return 0.0;

    double sum = 0.0;
    for (int t = t_min; t <= t_max; t++) {
        double sign = (t % 2 == 0) ? 1.0 : -1.0;
        double denom = factorial_int(t) * factorial_int(l - t)
                     * factorial_int(l - m1 - t) * factorial_int(l + m2 - t)
                     * factorial_int(m1 + t) * factorial_int(-m2 + t);
        sum += sign / denom;
    }

    return phase * sqrt(delta * mfact) * sum;
}


void calculate_w_single(py::dict& atoms,
    const int lm) {

    string key_real = "q" + to_string(lm) + "_real";
    string key_imag = "q" + to_string(lm) + "_imag";

    vector<vector<double>> q_real = atoms[py::str(key_real)].cast<vector<vector<double>>>();
    vector<vector<double>> q_imag = atoms[py::str(key_imag)].cast<vector<vector<double>>>();

    int nop = q_real.size();
    vector<double> w_values(nop, 0.0);
    vector<double> wbar_values(nop, 0.0);

    // For odd l, W_l = 0 (3j symbol vanishes when 3l is odd)
    if ((3 * lm) % 2 != 0) {
        atoms[py::str("w" + to_string(lm))] = w_values;
        atoms[py::str("what" + to_string(lm))] = wbar_values;
        return;
    }

    // Precompute all needed Wigner 3j symbols
    // m1 runs from -l to l, m2 from -l to l, m3 = -(m1+m2)
    vector<vector<double>> w3j_table((2 * lm + 1), vector<double>(2 * lm + 1, 0.0));
    for (int m1 = -lm; m1 <= lm; m1++) {
        for (int m2 = -lm; m2 <= lm; m2++) {
            int m3 = -(m1 + m2);
            if (abs(m3) <= lm) {
                w3j_table[m1 + lm][m2 + lm] = wigner3j(lm, m1, m2, m3);
            }
        }
    }

    for (int ti = 0; ti < nop; ti++) {
        double w_val = 0.0;
        double norm_sq = 0.0;

        // Compute |q_lm|^2 sum for normalization
        for (int mi = 0; mi < 2 * lm + 1; mi++) {
            norm_sq += q_real[ti][mi] * q_real[ti][mi] + q_imag[ti][mi] * q_imag[ti][mi];
        }

        // W_l = sum_{m1+m2+m3=0} (l l l / m1 m2 m3) * q_lm1 * q_lm2 * q_lm3
        // where q_lm is complex: q_real + i*q_imag
        for (int m1 = -lm; m1 <= lm; m1++) {
            int idx1 = m1 + lm;  // index into q arrays
            for (int m2 = -lm; m2 <= lm; m2++) {
                int m3 = -(m1 + m2);
                if (abs(m3) > lm) continue;
                int idx2 = m2 + lm;
                int idx3 = m3 + lm;

                double w3j = w3j_table[m1 + lm][m2 + lm];
                if (w3j == 0.0) continue;

                // Complex product: q1 * q2 * q3
                // (a1+ib1)(a2+ib2) = (a1a2-b1b2) + i(a1b2+a2b1)
                double r1 = q_real[ti][idx1], i1 = q_imag[ti][idx1];
                double r2 = q_real[ti][idx2], i2 = q_imag[ti][idx2];
                double r3 = q_real[ti][idx3], i3 = q_imag[ti][idx3];

                double r12 = r1 * r2 - i1 * i2;
                double i12 = r1 * i2 + i1 * r2;
                double r123 = r12 * r3 - i12 * i3;
                // imaginary part of triple product should be ~0 for real W_l
                // double i123 = r12 * i3 + i12 * r3;

                w_val += w3j * r123;
            }
        }

        w_values[ti] = w_val;
        // Normalized: W-hat_l = W_l / (sum |q_lm|^2)^(3/2)
        double norm_cubed = pow(norm_sq, 1.5);
        wbar_values[ti] = (norm_cubed > 1e-30) ? w_val / norm_cubed : 0.0;
    }

    atoms[py::str("w" + to_string(lm))] = w_values;
    atoms[py::str("what" + to_string(lm))] = wbar_values;
}


void calculate_aw_single(py::dict& atoms,
    const int lm) {

    // Compute averaged W_l using neighbor-averaged q_lm values
    // First, compute neighbor-averaged q_lm (same as calculate_aq_single but
    // we keep the full complex components, then compute W_l from those)

    string key_real = "q" + to_string(lm) + "_real";
    string key_imag = "q" + to_string(lm) + "_imag";

    vector<vector<double>> q_real = atoms[py::str(key_real)].cast<vector<vector<double>>>();
    vector<vector<double>> q_imag = atoms[py::str(key_imag)].cast<vector<vector<double>>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();

    int nop = q_real.size();

    // For odd l, W_l = 0
    if ((3 * lm) % 2 != 0) {
        vector<double> zeros(nop, 0.0);
        atoms[py::str("avg_w" + to_string(lm))] = zeros;
        atoms[py::str("avg_what" + to_string(lm))] = zeros;
        return;
    }

    // Precompute 3j table
    vector<vector<double>> w3j_table((2 * lm + 1), vector<double>(2 * lm + 1, 0.0));
    for (int m1 = -lm; m1 <= lm; m1++) {
        for (int m2 = -lm; m2 <= lm; m2++) {
            int m3 = -(m1 + m2);
            if (abs(m3) <= lm) {
                w3j_table[m1 + lm][m2 + lm] = wigner3j(lm, m1, m2, m3);
            }
        }
    }

    vector<double> w_values(nop, 0.0);
    vector<double> wbar_values(nop, 0.0);

    for (int ti = 0; ti < nop; ti++) {
        // First compute averaged q_lm for this atom
        int nns = neighbors[ti].size();
        vector<double> avg_real(2 * lm + 1, 0.0);
        vector<double> avg_imag(2 * lm + 1, 0.0);

        for (int mi = 0; mi < 2 * lm + 1; mi++) {
            avg_real[mi] = q_real[ti][mi];
            avg_imag[mi] = q_imag[ti][mi];
            for (int ci = 0; ci < nns; ci++) {
                avg_real[mi] += q_real[neighbors[ti][ci]][mi];
                avg_imag[mi] += q_imag[neighbors[ti][ci]][mi];
            }
            avg_real[mi] /= double(nns + 1);
            avg_imag[mi] /= double(nns + 1);
        }

        // Compute W_l from averaged q_lm
        double w_val = 0.0;
        double norm_sq = 0.0;
        for (int mi = 0; mi < 2 * lm + 1; mi++) {
            norm_sq += avg_real[mi] * avg_real[mi] + avg_imag[mi] * avg_imag[mi];
        }

        for (int m1 = -lm; m1 <= lm; m1++) {
            for (int m2 = -lm; m2 <= lm; m2++) {
                int m3 = -(m1 + m2);
                if (abs(m3) > lm) continue;
                int idx1 = m1 + lm, idx2 = m2 + lm, idx3 = m3 + lm;
                double w3j = w3j_table[m1 + lm][m2 + lm];
                if (w3j == 0.0) continue;

                double r12 = avg_real[idx1] * avg_real[idx2] - avg_imag[idx1] * avg_imag[idx2];
                double i12 = avg_real[idx1] * avg_imag[idx2] + avg_imag[idx1] * avg_real[idx2];
                double r123 = r12 * avg_real[idx3] - i12 * avg_imag[idx3];
                w_val += w3j * r123;
            }
        }

        w_values[ti] = w_val;
        double norm_cubed = pow(norm_sq, 1.5);
        wbar_values[ti] = (norm_cubed > 1e-30) ? w_val / norm_cubed : 0.0;
    }

    atoms[py::str("avg_w" + to_string(lm))] = w_values;
    atoms[py::str("avg_what" + to_string(lm))] = wbar_values;
}


void calculate_disorder(py::dict& atoms,
	const int lm){
    
    double sum2ti, sum2tj;
    double realdotproduct, imgdotproduct;
    double connection;
    double dis;

    string key1, key2;
	key1 = "q"+to_string(lm)+"_real";
	key2 = "q"+to_string(lm)+"_imag";

    vector<vector<double>> q_real = atoms[py::str(key1)].cast<vector<vector<double>>>();
    vector<vector<double>> q_imag = atoms[py::str(key2)].cast<vector<vector<double>>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<double> sii;
    vector<double> disorder;
    int nop = neighbors.size();

    for(int ti=0; ti<nop; ti++){

        sum2ti = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;

        for (int mi = 0;mi < 2*lm+1 ; mi++){
            sum2ti += q_real[ti][mi]*q_real[ti][mi] + q_imag[ti][mi]*q_imag[ti][mi];
            realdotproduct += q_real[ti][mi]*q_real[ti][mi];
            imgdotproduct  += q_imag[ti][mi]*q_imag[ti][mi];
        }
        connection = (realdotproduct+imgdotproduct)/(sqrt(sum2ti)*sqrt(sum2ti));
        sii.emplace_back(connection);
    }

    //first round is over
    //now find cross terms
    for(int ti=0; ti<nop; ti++){

        sum2ti = 0.0;
        sum2tj = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;
        dis = 0;

        for(int tj=0; tj<neighbors[ti].size(); tj++){
            for (int mi = 0; mi<2*lm+1 ; mi++){
                sum2ti += q_real[ti][mi]*q_real[ti][mi] + q_imag[ti][mi]*q_imag[ti][mi];
                sum2tj += q_real[tj][mi]*q_real[tj][mi] + q_imag[tj][mi]*q_imag[tj][mi];
                realdotproduct += q_real[ti][mi]*q_real[tj][mi];
                imgdotproduct  += q_imag[ti][mi]*q_imag[tj][mi];
            }
            connection = (realdotproduct+imgdotproduct)/(sqrt(sum2tj)*sqrt(sum2ti));
            dis += (sii[ti] + sii[tj] - 2*connection);
        }
        disorder.emplace_back(dis/float(neighbors[ti].size()));
    }

    atoms[py::str("disorder")] = disorder;
}