#include "Positioner.hh"
#include <cstdio>

//         (-2)(-1)( 0)  x
//-2       1   2   3 ( 1)
//-1     4   5   6   7 ( 2)
// 0   8   9  10  11  12
// 1    13  14  15  16
// 2     17  18  19
// y

Positioner::Positioner() {
    v[0][0] = 1.;
    v[0][1] = 0.;
    v[1][0] = 0.5;
    v[1][1] = sqrt(3.)/2;
    
    const int nrow = 5;
    const int y0 = -2;
    const int x0[nrow] = {0,-1,-2,-2,-2};
    const int xn[nrow] = {3,4,5,4,3};
    unsigned int n = 0;
    for(int y = y0; y < y0+nrow; y++) {
        for(int x = x0[y-y0]; x < x0[y-y0]+xn[y-y0]; x++) {
            vpos[0][n] = x;
            vpos[1][n] = y;
            ++n;
        }
    }
    
    for(n=0; n<N; n++) {
        pos[0][n] = vpos[0][n]*v[0][0] + vpos[1][n]*v[1][0];
        pos[1][n] = vpos[0][n]*v[0][1] + vpos[1][n]*v[1][1];
    }
}

void Positioner::calcPos(const double* L) {
    for(int i=0; i<2; i++) px[i] = sx[i] = 0;
    double w = 0;
    for(unsigned int n=0; n<N; n++) {
        w += L[n];
        for(int i=0; i<2; i++) {
            px[i] += L[n]*pos[i][n];
            sx[i] += L[n]*pos[i][n]*pos[i][n];
        }
    }
    if(w > 0) {
        for(int i=0; i<2; i++) {
            px[i] /= w;
            sx[i] = sx[i]/w - px[i]*px[i];
            sx[i] = sx[i]>0? sqrt(sx[i]) : 0;
        }
    } else {
        for(int i=0; i<2; i++) px[i] = sx[i] = -1000;
    }
}
