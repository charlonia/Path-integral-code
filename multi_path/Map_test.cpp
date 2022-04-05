//
// Created by charlonia on 2022/4/4.
//
//备注：pthread库不是Linux系统默认的库，连接时需要使用静态库libpthread.a，所以在线程函数在编译时，需要连接库函数，如上图    gcc pthread_create.c -o pthread_create -lpthread
#include <iostream>
#include <complex.h>
#include <math.h>
#include <array>
#include <fstream>
#include <pthread.h>
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
int main(){
    double grid =1;
    long double grid_div=0.01;
    double time_div = 0.002; //时间间隔
    double time_all = 0.003;     //总的演化时间
    d = 0.4/grid_div;          //狭缝宽度
    double v=10; // 设置电子入射速度a
    double hb=1;
    double m=1;
    totallen = round(grid/grid_div);
    const int totaltime = round(time_all/time_div);
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
    for(int j=0;j<totallen;j++) {
        for (int i = 0; i <totallen; i++) {
            A[j][i]=(complex<long double>){0,0};
            B[j][i]=(complex<long double>){0,0};
            C[j][i]=exp((complex<long double>) {0, z * (pow(j, 2) + pow(i, 2))})*(complex<long double>){pow(grid_div,2),0};
        }
    }
    int kx =1;
    int ky =1;
    int NUM_THREADS=100;
    pthread_t threads[NUM_THREADS];
    mape ma[NUM_THREADS][NUM_THREADS];
    pthread_attr_t attr;
    void *status;
//    pthread_attr_init(&attr);
//    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for(int i=0; i < NUM_THREADS; i++ ){
        cout <<"main() : creating thread, " << i << endl;
        ma[i][1].kx = kx;
        ma[i][1].ky = ky;
        auto rc = pthread_create(&threads[i], NULL,leftcalcu, (void *) &ma[i][1]);
        if (rc){
            cout << "Error:unable to create thread," << rc << endl;
            exit(-1);
        }
    }
//    pthread_attr_destroy(&attr);
    for( int i=0; i < NUM_THREADS; i++ ){
        auto rc = pthread_join(threads[i], &status);
        if (rc){
            cout << "Error:unable to join," << rc << endl;
            exit(-1);
        }
        cout << "Main: completed thread id :" << i ;
        cout << "  exiting with status :" << status << endl;
    }
    delete[] buffera,delete[] bufferb,delete[] bufferc,delete[]A,delete[]B,delete[]C;
    buffera = nullptr,bufferb= nullptr,bufferc= nullptr,A= nullptr,B= nullptr,C= nullptr;
    pthread_exit(NULL);

    return 0;
}
