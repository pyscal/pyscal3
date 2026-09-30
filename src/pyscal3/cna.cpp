/*
Common neighbour analysis (adaptive and with a lattice constant) and the
diamond structure identification.

The neighbours of each atom are its nearest candidates from the neighbour
search (neighbor_backend.cpp), with their bond vectors. Two neighbours m, n
of atom i are bonded if |v_m - v_n| <= cutoff_i, so no minimum image and no
padded supercell is needed. The neighbour graph of an atom is kept as
bitmasks (at most 16 neighbours).

Signature of neighbour k of atom i:
    c0  number of neighbours of i bonded to k (common neighbours)
    c1  number of bonds among these common neighbours
    c2  largest number of such bonds of one common neighbour
    c3  smallest number of such bonds of one common neighbour (8 if none)
fcc: 12 x (4,2,1,1); hcp: 6 x (4,2,1,1) + 6 x (4,2,2,0); ico: 12 x (5,5,2,2);
bcc (14 neighbours): 6 x (4,4,2,2) + 8 x (6,6,2,2).
Labels: 0 others, 1 fcc, 2 hcp, 3 bcc, 4 ico.
*/
#include "system.h"
#include "neighbor_backend.h"

#include <cmath>
#include <cstdint>

namespace {

using nlb::idx;
using nlb::Rows;

constexpr int MAXNB = 16;

inline int popcount(std::uint32_t x) {
    int c = 0;
    while (x) {
        x &= x - 1;
        c++;
    }
    return c;
}

// v: nn bond vectors (3 values each), nn <= MAXNB
void signatures(const double *v, int nn, double cutoff, int sig[][4]) {
    std::uint32_t adj[MAXNB] = {0};
    for (int a = 0; a + 1 < nn; a++) {
        for (int b = a + 1; b < nn; b++) {
            const double dx = v[3 * a] - v[3 * b];
            const double dy = v[3 * a + 1] - v[3 * b + 1];
            const double dz = v[3 * a + 2] - v[3 * b + 2];
            if (sqrt(dx * dx + dy * dy + dz * dz) <= cutoff) {
                adj[a] |= 1u << b;
                adj[b] |= 1u << a;
            }
        }
    }
    for (int k = 0; k < nn; k++) {
        const std::uint32_t common = adj[k];
        int bonds = 0, maxbonds = 0, minbonds = 8;
        for (int l = 0; l < nn; l++) {
            if (!(common >> l & 1u)) continue;
            const int b = popcount(adj[l] & common);
            bonds += b;
            maxbonds = std::max(maxbonds, b);
            minbonds = std::min(minbonds, b);
        }
        sig[k][0] = popcount(common);
        sig[k][1] = bonds / 2;
        sig[k][2] = maxbonds;
        sig[k][3] = minbonds;
    }
}

bool is(const int s[4], int a, int b, int c, int d) {
    return s[0] == a && s[1] == b && s[2] == c && s[3] == d;
}

// 1 fcc, 2 hcp, 4 ico or 0 from the signatures of 12 (or more) neighbours
int classify_cn12(const int sig[][4], int nn) {
    int nfcc = 0, nhcp = 0, nico = 0;
    for (int k = 0; k < nn; k++) {
        if (is(sig[k], 4, 2, 1, 1)) nfcc++;
        else if (is(sig[k], 4, 2, 2, 0)) nhcp++;
        else if (is(sig[k], 5, 5, 2, 2)) nico++;
    }
    if (nfcc == 12) return 1;
    if (nfcc == 6 && nhcp == 6) return 2;
    if (nico == 12) return 4;
    return 0;
}

// 3 bcc or 0 from the signatures of 14 neighbours
int classify_cn14(const int sig[][4], int nn) {
    int nbcc1 = 0, nbcc2 = 0;
    for (int k = 0; k < nn; k++) {
        if (is(sig[k], 4, 4, 2, 2)) nbcc1++;
        else if (is(sig[k], 6, 6, 2, 2)) nbcc2++;
    }
    return (nbcc1 == 6 && nbcc2 == 8) ? 3 : 0;
}

py::array_t<std::int64_t> to_array(const std::vector<std::int64_t> &s) {
    py::array_t<std::int64_t> out(static_cast<py::ssize_t>(s.size()));
    std::copy(s.begin(), s.end(), out.mutable_data());
    return out;
}

bool enough_candidates(const Rows &c, idx n, int nmin) {
    for (idx a = 0; a < n; a++)
        if (c.offsets[a + 1] - c.offsets[a] < nmin) return false;
    return true;
}

}  // namespace

py::tuple cna_structure(const nl_positions &positions, const nl_positions &cell,
                        const vector<bool> &pbc, double prefactor,
                        double lattice_constant, int nmin) {
    // adaptive cutoffs if lattice_constant is 0; returns (labels, whether
    // every atom had at least nmin candidates)
    const nlb::Geometry g = nlb::make_geometry(positions, cell, pbc);
    // only the 14 nearest candidates are used, so most atoms are searched
    // with a smaller radius (the result is the same)
    const double guess = nlb::guess_radius(g, prefactor);
    const Rows c = nlb::nearest_candidates(g, 0.85 * guess, guess, 14);
    std::vector<std::int64_t> structure(g.n, 0);
    int sig[MAXNB][4];

    for (idx ti = 0; ti < g.n; ti++) {
        const idx lo = c.offsets[ti];
        const idx cnt = c.offsets[ti + 1] - lo;
        const double *v = &c.v[3 * lo];
        const double *d = &c.d[lo];

        if (cnt >= 12) {
            double cutoff;
            if (lattice_constant == 0.0) {
                double ssum = 0;
                for (int i = 0; i < 12; i++) ssum += d[i];
                cutoff = 1.207 * ssum / 12.00;
            } else {
                cutoff = 0.854 * lattice_constant;
            }
            signatures(v, 12, cutoff, sig);
            structure[ti] = classify_cn12(sig, 12);
        }
        if (structure[ti] == 0 && cnt >= 14) {
            double cutoff;
            if (lattice_constant == 0.0) {
                double ssum = 0;
                for (int i = 0; i < 8; i++) ssum += 1.1547 * d[i];
                for (int i = 8; i < 14; i++) ssum += d[i];
                cutoff = 1.207 * ssum / 14.00;
            } else {
                cutoff = 1.207 * lattice_constant;
            }
            signatures(v, 14, cutoff, sig);
            structure[ti] = classify_cn14(sig, 14);
        }
    }
    return py::make_tuple(to_array(structure), enough_candidates(c, g.n, nmin));
}

py::tuple diamond_structure_cna(const nl_positions &positions, const nl_positions &cell,
                                const vector<bool> &pbc, double prefactor) {
    // Diamond identification: CNA on the second shell (the four nearest
    // candidates of each of the four nearest candidates, without the atom
    // itself). Labels: 0 others, 1 cubic diamond, 2 cubic diamond first
    // neighbour, 3 cubic diamond second neighbour, 4 hexagonal diamond,
    // 5 and 6 its first and second neighbours.
    const nlb::Geometry g = nlb::make_geometry(positions, cell, pbc);
    // only the 4 nearest candidates are used (see cna_structure)
    const double guess = nlb::guess_radius(g, prefactor);
    const Rows c = nlb::nearest_candidates(g, 0.65 * guess, guess, 4);
    const idx n = g.n;
    std::vector<std::int64_t> structure(n, 0);
    // second-shell neighbour indices of each atom, at most MAXNB
    std::vector<idx> second(n * MAXNB);
    std::vector<int> nsecond(n, 0);
    int sig[MAXNB][4];

    for (idx ti = 0; ti < n; ti++) {
        const idx lo = c.offsets[ti];
        if (c.offsets[ti + 1] - lo < 4) continue;
        double v2[3 * MAXNB], d2[MAXNB];
        int nn = 0;
        for (int i = 0; i < 4; i++) {
            const idx tj = c.j[lo + i];
            const idx lj = c.offsets[tj];
            if (c.offsets[tj + 1] - lj < 4) continue;
            for (int k = 0; k < 4; k++) {
                const idx tk = c.j[lj + k];
                // r_i - r_k through the bond vectors, which picks the right image
                const double x = c.v[3 * (lo + i)] + c.v[3 * (lj + k)];
                const double y = c.v[3 * (lo + i) + 1] + c.v[3 * (lj + k) + 1];
                const double z = c.v[3 * (lo + i) + 2] + c.v[3 * (lj + k) + 2];
                const double dist = sqrt(x * x + y * y + z * z);
                if (tk == ti && dist < 1e-8) continue;   // back to the atom itself
                v2[3 * nn] = x;
                v2[3 * nn + 1] = y;
                v2[3 * nn + 2] = z;
                d2[nn] = dist;
                second[ti * MAXNB + nn] = tk;
                nn++;
            }
        }
        nsecond[ti] = nn;
        if (nn < 12) continue;
        double ssum = 0;
        for (int i = 0; i < 12; i++) ssum += d2[i];
        signatures(v2, nn, 1.207 * ssum / 12.00, sig);
        const int s = classify_cn12(sig, nn);
        structure[ti] = (s == 1) ? 5 : (s == 2) ? 8 : s;
    }

    // first neighbours of a diamond site
    for (idx ti = 0; ti < n; ti++) {
        if (structure[ti] >= 5) continue;
        const idx lo = c.offsets[ti];
        // an atom with fewer than four candidates has no first shell
        const idx nf = (c.offsets[ti + 1] - lo >= 4) ? 4 : 0;
        for (idx i = 0; i < nf; i++) {
            const std::int64_t sj = structure[c.j[lo + i]];
            if (sj == 5) { structure[ti] = 6; break; }
            if (sj == 8) { structure[ti] = 9; break; }
        }
    }
    // second neighbours of a diamond site
    for (idx ti = 0; ti < n; ti++) {
        if (structure[ti] >= 5) continue;
        for (int k = 0; k < nsecond[ti]; k++) {
            const std::int64_t sk = structure[second[ti * MAXNB + k]];
            if (sk == 5) { structure[ti] = 7; break; }
            if (sk == 8) { structure[ti] = 10; break; }
        }
    }
    for (idx ti = 0; ti < n; ti++) {
        if (structure[ti] >= 5) structure[ti] -= 4;
    }
    return py::make_tuple(to_array(structure), enough_candidates(c, n, 4));
}
