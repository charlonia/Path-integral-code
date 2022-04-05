//
// Created by charlonia on 2022/4/4.
//
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
#include <pthread.h>
#include<sstream>
using namespace std;
int totallen;
complex<long double>**A;
complex<long double>**B;
complex<long double>**C;
long double d;
struct mape{
    int ky;
    int kx;
};
void * leftcalcu(void* mapee
)
{
    mape *mapo;
    mapo = (mape*) mapee;
    int & ky = mapo->ky;
    int & kx = mapo->kx;
    long double lo;
    long double hi;
    int high;
    int low;
    for (int lx = 0; lx < round(totallen / 2) - 1; ++lx) {
        for (int ly = 0; ly < totallen; ++ly) {
            B[ky][kx] +=A[ly][lx] * C[abs(ly-ky)][abs(lx-kx)];
        }
    }
    for (int lx = round(totallen / 2)-1; lx < totallen - 1; ++lx) {
        lo = ((totallen - d) / 2 - ky) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
        low = floor(lo);
        if (abs( round(lo)- (lo))<1/(long double)totallen){
            low = round(lo);
        }
        low = max(low, 1);
        hi = ((totallen + d) / 2 - ky-1) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
        high = ceil(hi);
        if (abs(round(hi)-(hi))<1/(long double)totallen){
            high = round(hi);
        }
        high = min(high, totallen);
        for (int ly = low-1; ly < high; ++ly) {
            B[ky][kx] +=A[ly][lx] * C[abs(ly-ky)][abs(lx-kx)];
        }
    }
    return 0;
}
void *rightcalcu(void* mapee
){
    mape *mapo;
    mapo = (mape*) mapee;
    int & ky = mapo->ky;
    int & kx = mapo->kx;
    long double lo;
    long double hi;
    int high;
    int low;
    for (int lx = round(totallen/2); lx < totallen - 1; ++lx) {
        for (int ly = 0; ly < totallen; ++ly) {
            B[ky][kx] +=A[ly][lx] * C[abs(ly-ky)][abs(lx-kx)];
        }
    }
    for (int lx = 0; lx <round(totallen / 2); ++lx) {
        lo = ((totallen - d) / 2 - ky) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
        low = floor(lo);
        if (abs( round(lo)- (lo))<1/(long double)totallen){
            low = round(lo);
        }
        low = max(low, 1);
        hi = ((totallen + d) / 2 - ky-1) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
        high = ceil(hi);
        if (abs(round(hi)-(hi))<1/(long double)totallen){
            high = round(hi);
        }
        high = min(high, totallen);
        for (int ly = low-1; ly < high; ++ly) {
            B[ky][kx] +=A[ly][lx] * C[abs(ly-ky)][abs(lx-kx)];
        }
    }
    return 0;
}
int main(){
    clock_t start_time=clock();
    cout << "waiting for input"<<endl;
    char  watoin;
    cin>>watoin;
    double grid =1;
    long double grid_div=0.01;
    double time_div = 0.002; //时间间隔
    double time_all = 0.003;     //总的演化时间
    d = 0.4/grid_div;          //狭缝宽度
    double v=10; // 设置电子入射速度a
    double hb=1;
    double m=1;
    if (watoin == 'y'){
        cout<<"grid total length"<<endl;
        cin>>grid;
        cout<<"grid division"<<endl;
        cin>>grid_div;
        cout<<"total time length"<<endl;
        cin>>time_all;
        cout<<"time division"<<endl;
        cin>>time_div;
        cout<<"slit width"<<endl;
        cin>>d;
        d =d/grid_div;
        cout<<"hbar"<<endl;
        cin>>hb;
        cout<<"speed of the particle"<<endl;
        cin>>v;
        cout<<"mass of the particle"<<endl;
        cin>>m;
    }if (watoin == '1'){
        cout<<"grid total length"<<endl;
        cin>>grid;
        cout<<"grid division"<<endl;
        cin>>grid_div;
        cout<<"total time length"<<endl;
        cin>>time_all;
        cout<<"time division"<<endl;
        cin>>time_div;
    }
    clock_t start_time1=clock();
    totallen = round(grid/grid_div);
    const int totaltime = round(time_all/time_div);
    int midx = round(totallen/2)-1;
    A = new complex<long double>*[totallen];
    B = new complex<long double>*[totallen];
    C = new complex<long double>*[totallen];
    complex<long double> *buffera = new complex<long double>[totallen*totallen];
    complex<long double> *bufferb = new complex<long double>[totallen*totallen];
    complex<long double> *bufferc = new complex<long double>[totallen*totallen];
    for(int i=0;i<totallen;i++) {
        A[i] = buffera + i * totallen;
        B[i] = bufferb + i * totallen;
        C[i] = bufferc + i * totallen;
    }
    const long double z=-pow(grid_div,2)/(2*time_div)*m/hb ;
    cout << "start calculation"<<endl;
    mape ma[totallen][totallen];
    for(int j=0;j<totallen;j++) {
        for (int i = 0; i <totallen; i++) {
            A[j][i]=(complex<long double>){0,0};
            B[j][i]=(complex<long double>){0,0};
            C[j][i]=exp((complex<long double>) {0, z * (pow(j, 2) + pow(i, 2))})*(complex<long double>){pow(grid_div,2),0};
            ma[j][i].kx=i;
            ma[j][i].ky=j;
        }
    }
    for(int j=0;j<totallen;j++) {
        for (int i = round(totallen/5)-1; i <round(2*totallen/5)-1; i++) {
            A[j][i]=exp((complex<long double>){0,(i+1)*m*v*grid_div/hb});
        }
    }
    int NUM_THREADS =100;
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    void *status;
//    pthread_attr_init(&attr);
//    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
//    这里+1是为了能和matlab程序相互对应
    clock_t end_time1=clock();
    cout << "initialization"<<"The run time is: " <<(double)(end_time1 - start_time1) / CLOCKS_PER_SEC <<endl;
    for (int j=0;j<totaltime;++j){
        clock_t start_time=clock();
        for (int kx=0;kx < round(totallen/2)-1;++kx) {
            for (int ky = 0; ky < totallen; ++ky) {
                auto rc = pthread_create(&threads[ky], NULL, leftcalcu, (void *) &ma[ky][kx]);
                if (rc) {
                    cout << "Error:unable to create thread," << rc << endl;
                    exit(-1);
                }
            }
            for (int ky = 0; ky < totallen; ++ky) {
                auto rc = pthread_join(threads[ky], &status);
                if (rc) {
                    cout << "Error:unable to join," << rc << endl;
                    exit(-1);
                }
            }
        }
        for (int ky=round((totallen - d) / 2);ky<round((totallen + d) / 2);++ky){
            for(int lx = 0;lx < totallen - 1;++lx){
                for(int ly = 0;ly < totallen ;++ly){
                    B[ky][midx] +=A[ly][lx] * C[abs(ly-ky)][abs(lx-midx)];
                }
            }
        }

        for (int kx=round(totallen/2);kx < totallen;++kx){
            for (int ky = 0; ky < totallen; ++ky) {
                auto rc = pthread_create(&threads[ky], NULL,rightcalcu, (void *) &ma[ky][kx]);
                if (rc){
                    cout << "Error:unable to create thread," << rc << endl;
                    exit(-1);
                }
            }
            for (int ky = 0; ky < totallen; ++ky) {
                auto rc = pthread_join(threads[ky], &status);
                if (rc){
                    cout << "Error:unable to join," << rc << endl;
                    exit(-1);
                }
            }
        }

        for (int y =0;y<totallen;++y){
            for(int x = 0; x<totallen;++x){
                A[y][x]=B[y][x];
            }
        }

        for (int y =0;y<totallen;++y){
            for(int x = 0; x<totallen-1;++x){
                B[y][x]=(complex<long double>){0,0};
            }
        }

        clock_t end_time=clock();
        cout << j<<" loop "<<"The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    }
    clock_t end_time=clock();
    cout <<"Total "<< "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    ofstream opt;
    opt.precision(20);
    stringstream ss1;
    ss1 << grid<< "_"<<grid_div<<"g_"<<time<< "_" <<time_div<<"t_"<<d<<"d" <<".csv";
    string file = ss1.str();
    cout<<file<<endl;
    opt.open(file);
    for(int j=0;j<totallen;j++) {
        for (int i = 0; i <totallen; i++) {
            opt << real(A[j][i]);
            if (imag(A[j][i])<0){
                opt<<imag(A[j][i])<<'i';}
            else{
                opt<<'+'<<imag(A[j][i])<<'i';}
            opt<<',';

        }
        opt<<endl;
    }
    delete[] buffera,delete[] bufferb,delete[] bufferc,delete[]A,delete[]B,delete[]C;
    buffera = nullptr,bufferb= nullptr,bufferc= nullptr,A= nullptr,B= nullptr,C= nullptr;
//    pthread_exit(NULL);
}
