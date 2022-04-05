//
// Created by charlonia on 2022/4/4.
//

#include <iostream>
#include <complex.h>
#include <cstdio>
#include <math.h>
#include <array>
#include <fstream>
#include <omp.h>
using namespace std;
int main(){
    clock_t c1 = clock();
    int l =200;
    complex<long double>A[l][l];
    for(int i=0;i<l;++i){

        for(int j =0;j<l;++j){
            A[i][j]=complex<long double>{static_cast<long double>(i),static_cast<long double>(j)};
        }
    }
    omp_set_num_threads(1);
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < l; ++j) {
                A[i][j] = pow(A[i][j], 2);
            }
        }
    }
    clock_t c2 = clock();
    cout<<(double)(c2 - c1) / CLOCKS_PER_SEC <<endl;

};