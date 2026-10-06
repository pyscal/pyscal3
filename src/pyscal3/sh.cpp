#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "system.h"
#include "parallel.h"

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
				realti = realti/weightsum;
				imgti = imgti/weightsum;

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
Spherical harmonics of one bond, all m at once
**********************************************************************/
vector<double> ylm_norms(const int l){
    // sqrt((2l + 1) / (4 pi) * (l - m)! / (l + m)!) for m = 0 .. l
    vector<double> norm(l + 1);
    for (int m = 0; m <= l; m++)
        norm[m] = sqrt(((2.0*double(l) + 1.0)/ (4.0*PI))*dfactorial(l, m));
    return norm;
}

void ylm_all_m(const int l,
    const double theta,
    const double phi,
    const vector<double>& norm,
    double* ylm_real,
    double* ylm_imag){

    // Y_lm for m = -l .. l, stored at index m + l. For each m the associated
    // Legendre function P_l^m(cos theta) uses the usual recurrence in l,
    // cos(m phi) and sin(m phi) follow from the angle-addition formulas.
    const double x = cos(theta);
    const double somx2 = sqrt((1.0 - x)*(1.0 + x));
    const double c1 = cos(phi);
    const double s1 = sin(phi);
    double cm = 1.0, sm = 0.0;
    double pmm = 1.0, fact = 1.0;

    for (int m = 0; m <= l; m++){
        if (m > 0){
            pmm *= -fact*somx2;
            fact += 2.0;
            const double c = cm*c1 - sm*s1;
            sm = sm*c1 + cm*s1;
            cm = c;
        }
        double p;
        if (l == m){
            p = pmm;
        }
        else{
            double pa = pmm, pb = x*(2*m + 1)*pmm;
            for (int ll = m + 2; ll <= l; ll++){
                const double pc = (x*(2*ll - 1)*pb - (ll + m - 1)*pa)/(ll - m);
                pa = pb;
                pb = pc;
            }
            p = pb;
        }
        p *= norm[m];
        ylm_real[l + m] = p*cm;
        ylm_imag[l + m] = p*sm;
        if (m > 0){
            // Condon-Shortley phase for negative m:  Y_{l,-m} = (-1)^m conj(Y_{lm}).
            // It cancels in |q_lm|^2 (q_l) but not in the Wigner-3j contraction
            // used for W_l, which is not rotationally invariant without it.
            const double sign = (m % 2 == 1) ? -1.0 : 1.0;
            ylm_real[l - m] = sign*p*cm;
            ylm_imag[l - m] = -sign*p*sm;
        }
    }
}


py::tuple calculate_q_single(const nl_index& offsets,
    const nl_values& theta,
    const nl_values& phi,
    const nl_values& weights,
    const int lm){

    // offsets[i] .. offsets[i + 1] are the bonds of atom i in theta, phi, weights
    const std::int64_t* off = offsets.data();
    const double* th = theta.data();
    const double* ph = phi.data();
    const double* w = weights.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    const int nm = 2*lm + 1;

    py::array_t<double> q(nop);
    py::array_t<double> qlm_real(vector<py::ssize_t>{nop, nm});
    py::array_t<double> qlm_img(vector<py::ssize_t>{nop, nm});
    double* qo = q.mutable_data();
    double* qr = qlm_real.mutable_data();
    double* qi = qlm_img.mutable_data();

    const vector<double> norm = ylm_norms(lm);

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    vector<double> ylm_real(nm), ylm_imag(nm), sum_real(nm), sum_imag(nm);
    double summ, weightsum;
    double realti, imgti;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        fill(sum_real.begin(), sum_real.end(), 0.0);
        fill(sum_imag.begin(), sum_imag.end(), 0.0);
        weightsum = 0;
        for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
            ylm_all_m(lm, th[ci], ph[ci], norm, ylm_real.data(), ylm_imag.data());
            for (int k=0; k<nm; k++){
                sum_real[k] += w[ci]*ylm_real[k];
                sum_imag[k] += w[ci]*ylm_imag[k];
            }
            weightsum += w[ci];
        }
        summ = 0;
        for (int k=0; k<nm; k++){
            realti = sum_real[k]/weightsum;
            imgti = sum_imag[k]/weightsum;
            qr[ti*nm + k] = realti;
            qi[ti*nm + k] = imgti;
            summ += realti*realti + imgti*imgti;
        }
        qo[ti] = pow(((4.0*PI/(2*lm+1))*summ),0.5);
    }
    }, 256);
    }
    return py::make_tuple(q, qlm_real, qlm_img);
}

py::array_t<double> calculate_aq_single(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm){

    // average q_lm over each atom and its neighbours
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* qr = q_real.data();
    const double* qi = q_imag.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    const int nm = 2*lm + 1;

    py::array_t<double> q(nop);
    double* qo = q.mutable_data();

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    double realti, imgti, summ;
    int nns;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        summ = 0;
        for (int mi=0; mi<nm; mi++){
            realti = qr[ti*nm + mi];
            imgti = qi[ti*nm + mi];
            nns = 0;
            for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
                realti += qr[nb[ci]*nm + mi];
                imgti += qi[nb[ci]*nm + mi];
                nns += 1;
            }
            realti = realti/(double(nns+1));
            imgti = imgti/(double(nns+1));
            summ += realti*realti + imgti*imgti;
        }
        qo[ti] = pow(((4.0*PI/(2*lm+1)) * summ),0.5);
    }
    }, 256);
    }
    return q;
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


// 3j symbols (l l l; m1 m2 -(m1+m2)), indexed [m1 + l][m2 + l]
static vector<vector<double>> wigner3j_table(const int lm) {
    vector<vector<double>> w3j_table((2 * lm + 1), vector<double>(2 * lm + 1, 0.0));
    for (int m1 = -lm; m1 <= lm; m1++) {
        for (int m2 = -lm; m2 <= lm; m2++) {
            int m3 = -(m1 + m2);
            if (abs(m3) <= lm) {
                w3j_table[m1 + lm][m2 + lm] = wigner3j(lm, m1, m2, m3);
            }
        }
    }
    return w3j_table;
}

// W_l = sum_{m1+m2+m3=0} (l l l / m1 m2 m3) * q_lm1 * q_lm2 * q_lm3,
// with q_lm = re[m + l] + i im[m + l]
static double wigner_w(const int lm, const double* re, const double* im,
    const vector<vector<double>>& w3j_table) {
    double w_val = 0.0;
    for (int m1 = -lm; m1 <= lm; m1++) {
        int idx1 = m1 + lm;
        for (int m2 = -lm; m2 <= lm; m2++) {
            int m3 = -(m1 + m2);
            if (abs(m3) > lm) continue;
            int idx2 = m2 + lm;
            int idx3 = m3 + lm;
            double w3j = w3j_table[m1 + lm][m2 + lm];
            if (w3j == 0.0) continue;
            double r12 = re[idx1] * re[idx2] - im[idx1] * im[idx2];
            double i12 = re[idx1] * im[idx2] + im[idx1] * re[idx2];
            double r123 = r12 * re[idx3] - i12 * im[idx3];
            w_val += w3j * r123;
        }
    }
    return w_val;
}

py::tuple calculate_w_single(const nl_values& q_real,
    const nl_values& q_imag,
    const int lm) {

    // W_l and its normalised form from the q_lm of each atom, (nop, 2l+1)
    const py::ssize_t nop = q_real.shape(0);
    const int nm = 2*lm + 1;
    const double* qr = q_real.data();
    const double* qi = q_imag.data();
    py::array_t<double> w_values(nop), wbar_values(nop);
    double* wv = w_values.mutable_data();
    double* wb = wbar_values.mutable_data();
    fill(wv, wv + nop, 0.0);
    fill(wb, wb + nop, 0.0);

    // For odd l, W_l = 0 (3j symbol vanishes when 3l is odd)
    if ((3 * lm) % 2 != 0) {
        return py::make_tuple(w_values, wbar_values);
    }

    const vector<vector<double>> w3j_table = wigner3j_table(lm);

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        const double* re = qr + ti*nm;
        const double* im = qi + ti*nm;
        wv[ti] = wigner_w(lm, re, im, w3j_table);
        double norm_sq = 0.0;
        for (int mi = 0; mi < nm; mi++) {
            norm_sq += re[mi] * re[mi] + im[mi] * im[mi];
        }
        // Normalized: W-hat_l = W_l / (sum |q_lm|^2)^(3/2)
        double norm_cubed = pow(norm_sq, 1.5);
        wb[ti] = (norm_cubed > 1e-30) ? wv[ti] / norm_cubed : 0.0;
    }
    }, 256);
    }
    return py::make_tuple(w_values, wbar_values);
}


py::tuple calculate_aw_single(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm) {

    // W_l of the q_lm averaged over each atom and its neighbours
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* qr = q_real.data();
    const double* qi = q_imag.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    const int nm = 2*lm + 1;
    py::array_t<double> w_values(nop), wbar_values(nop);
    double* wv = w_values.mutable_data();
    double* wb = wbar_values.mutable_data();
    fill(wv, wv + nop, 0.0);
    fill(wb, wb + nop, 0.0);

    if ((3 * lm) % 2 != 0) {
        return py::make_tuple(w_values, wbar_values);
    }

    const vector<vector<double>> w3j_table = wigner3j_table(lm);

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    vector<double> avg_real(nm), avg_imag(nm);
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        const std::int64_t nns = off[ti+1] - off[ti];
        for (int mi = 0; mi < nm; mi++) {
            avg_real[mi] = qr[ti*nm + mi];
            avg_imag[mi] = qi[ti*nm + mi];
            for (std::int64_t ci = off[ti]; ci < off[ti+1]; ci++) {
                avg_real[mi] += qr[nb[ci]*nm + mi];
                avg_imag[mi] += qi[nb[ci]*nm + mi];
            }
            avg_real[mi] /= double(nns + 1);
            avg_imag[mi] /= double(nns + 1);
        }

        double norm_sq = 0.0;
        for (int mi = 0; mi < nm; mi++) {
            norm_sq += avg_real[mi] * avg_real[mi] + avg_imag[mi] * avg_imag[mi];
        }
        wv[ti] = wigner_w(lm, avg_real.data(), avg_imag.data(), w3j_table);
        double norm_cubed = pow(norm_sq, 1.5);
        wb[ti] = (norm_cubed > 1e-30) ? wv[ti] / norm_cubed : 0.0;
    }
    }, 256);
    }
    return py::make_tuple(w_values, wbar_values);
}


py::array_t<double> calculate_disorder(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm){

    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* qr = q_real.data();
    const double* qi = q_imag.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    const int nm = 2*lm + 1;

    vector<double> sii(nop);
    py::array_t<double> disorder(nop);
    double* dout = disorder.mutable_data();

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    double sum2ti;
    double realdotproduct, imgdotproduct;
    double connection;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        sum2ti = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;
        for (int mi = 0; mi < nm; mi++){
            const double a = qr[ti*nm + mi], b = qi[ti*nm + mi];
            sum2ti += a*a + b*b;
            realdotproduct += a*a;
            imgdotproduct  += b*b;
        }
        connection = (realdotproduct+imgdotproduct)/(sqrt(sum2ti)*sqrt(sum2ti));
        sii[ti] = connection;
    }
    }, 256);
    }

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    double sum2ti, sum2tj;
    double realdotproduct, imgdotproduct;
    double connection;
    double dis;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        dis = 0;
        for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
            const std::int64_t tj = nb[ci];
            sum2ti = 0.0;
            sum2tj = 0.0;
            realdotproduct = 0.0;
            imgdotproduct = 0.0;
            for (int mi = 0; mi < nm; mi++){
                sum2ti += qr[ti*nm + mi]*qr[ti*nm + mi] + qi[ti*nm + mi]*qi[ti*nm + mi];
                sum2tj += qr[tj*nm + mi]*qr[tj*nm + mi] + qi[tj*nm + mi]*qi[tj*nm + mi];
                realdotproduct += qr[ti*nm + mi]*qr[tj*nm + mi];
                imgdotproduct  += qi[ti*nm + mi]*qi[tj*nm + mi];
            }
            connection = (realdotproduct+imgdotproduct)/(sqrt(sum2tj)*sqrt(sum2ti));
            dis += (sii[ti] + sii[tj] - 2*connection);
        }
        const std::int64_t nn = off[ti+1] - off[ti];
        dout[ti] = (nn > 0) ? dis/double(nn) : 0.0;
    }
    }, 256);
    }
    return disorder;
}