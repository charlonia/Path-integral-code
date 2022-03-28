//
// Created by charlonia on 2022/3/26.
//
//
// Created by charlonia on 2022/3/20.
//
#include <iostream>
#include <complex.h>
#include <math.h>
#include <array>
#include <fstream>
using namespace std;
int main(){
    clock_t start_time=clock();
    char  watoin;
    cin>>watoin;
    double grid =1;
    long double grid_div=0.01;
    double time_div = 0.002; //时间间隔
    double time_all = 1;     //总的演化时间
    long double d = 0.4/grid_div;          //狭缝宽度
    double v=10; // 设置电子入射速度
    double hb=1;
    double m=1;
    if (watoin == 'y'){
        cout<<"网格长度"<<endl;
        cin>>grid;
        cout<<"网格间隔"<<endl;
        cin>>grid_div;
        cout<<"时间长度"<<endl;
        cin>>time_all;
        cout<<"时间间隔"<<endl;
        cin>>time_div;
        cout<<"狭缝宽度"<<endl;
        cin>>d;
        d =d/grid_div;
        cout<<"hbar"<<endl;
        cin>>hb;
        cout<<"粒子速度"<<endl;
        cin>>v;
        cout<<"粒子质量"<<endl;
        cin>>m;
    }
    const int totallen = round(grid/grid_div);
    const int totaltime = round(time_all/time_div);
    int low;
    int high;
    int midx = round(totallen/2)-1;
    complex<long double> A[totallen][totallen];
    complex<long double> B[totallen][totallen];
    const long double z=-pow(grid_div,2)/(2*time_div)*m/hb ;
    for (int j=0;j<totaltime;++j){
        for (int kx=0;kx < round(totallen/2)-1;++kx) {
            for (int ky = 0; ky < totallen; ++ky) {
                for (int lx = 0; lx < round(totallen / 2) - 1; ++lx) {
                    for (int ly = 0; ly < totallen; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<long double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<long double>){pow(grid_div,2),0};
                    }
                }
                for (int lx = round(totallen / 2)-1; lx < totallen - 1; ++lx) {
                    low = round(((totallen - d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    low = max(low, 0);
                    high = round(((totallen + d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    high = min(high, totallen);
                    for (int ly = low-1; ly < high; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<long double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<long double>){pow(grid_div,2),0};
                    }
                }
            }
        }
        for (int ky=round((totallen - d) / 2);ky<round((totallen + d) / 2);++ky){
            for(int lx = 0;lx < totallen - 1;++lx){
                for(int ly = 0;ly < totallen ;++ly){
                    B[ky][midx] +=A[ly][lx] * exp((complex<long double>){0, z * (pow(lx - midx, 2) + pow(ky - ly, 2))})*(complex<long double>){pow(grid_div,2),0};

                }
            }
        }
        for (int kx=round(totallen/2);kx < totallen;++kx){
            for(int ky=0;ky < totallen;++ky){
                for (int lx = round(totallen/2); lx < totallen - 1; ++lx) {
                    for (int ly = 0; ly < totallen; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<long double>){0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<long double>){pow(grid_div,2),0};
                    }
                }
                for (int lx = 0; lx <round(totallen / 2); ++lx) {
                    low = round(((totallen - d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    low = max(low, 0);
                    high = round(((totallen + d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    high = min(high, totallen);
                    for (int ly = low; ly < high; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<long double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<long double>){pow(grid_div,2),0};
                    }
                }
            }
        }
        for (int y =0;totallen;++y){
            for(int x = 0; totallen;++x){
                A[y][x]=B[y][x];
            }
        }
    }
    clock_t end_time=clock();
    cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    ofstream outfile("output.xls");
    for (int i = 1; i <= 10; i++)
        outfile << i << '\t' << "Hello" << endl;
    return 0;
}

