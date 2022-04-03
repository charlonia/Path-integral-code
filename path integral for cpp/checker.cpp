//
// Created by charlonia on 2022/4/2.
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
using namespace std;
int main(){
    clock_t start_time=clock();
    cout << "�ȴ�����"<<endl;
    char  watoin;
    cin>>watoin;
    double grid =1;
    long double grid_div=0.01;
    double time_div = 0.002; //ʱ����
    double time_all = 0.001;     //�ܵ��ݻ�ʱ��
    long double d = 0.4/grid_div;          //������
    double v=10; // ���õ��������ٶ�a
    double hb=1;
    double m=1;
    if (watoin == 'y'){
        cout<<"���񳤶�"<<endl;
        cin>>grid;
        cout<<"������"<<endl;
        cin>>grid_div;
        cout<<"ʱ�䳤��"<<endl;
        cin>>time_all;
        cout<<"ʱ����"<<endl;
        cin>>time_div;
        cout<<"������"<<endl;
        cin>>d;
        d =d/grid_div;
        cout<<"hbar"<<endl;
        cin>>hb;
        cout<<"�����ٶ�"<<endl;
        cin>>v;
        cout<<"��������"<<endl;
        cin>>m;
    }
    const int totallen = round(grid/grid_div);
    const int totaltime = round(time_all/time_div);
    int low;
    int high;
    int midx = round(totallen/2)-1;
    long double hi;
    long double lo;
    complex<long double> A[totallen][totallen];
    complex<long double> B[totallen][totallen];
    complex<long double> C[totallen][totallen];
    const long double z=-pow(grid_div,2)/(2*time_div)*m/hb ;
    cout << "���㿪ʼ"<<endl;
    clock_t start_time1=clock();
    for(int j=0;j<totallen;j++) {
        for (int i = 0; i <totallen; i++) {
            A[j][i]=(complex<long double>){0,0};
            B[j][i]=(complex<long double>){0,0};
            C[j][i]=exp((complex<long double>) {0, z * (pow(j, 2) + pow(i, 2))})*(complex<long double>){pow(grid_div,2),0};
        }
    }
    for(int j=0;j<totallen;j++) {
        for (int i = round(totallen/5)-1; i <round(2*totallen/5)-1; i++) {
            A[j][i]=exp((complex<long double>){0,(i+1)*m*v*grid_div/hb});
        }
    }
//    ����+1��Ϊ���ܺ�matlab�����໥��Ӧ
    clock_t end_time1=clock();
    cout << "��ʼ��"<<" ѭ�� "<<"The run time is: " <<(double)(end_time1 - start_time1) / CLOCKS_PER_SEC <<endl;
    for (int j=0;j<totaltime;++j){
        cout <<j<<endl;
        clock_t start_time=clock();
        int ky =23;
        int kx =70;

        for (int lx = round(totallen / 2)-1; lx < totallen - 1; ++lx) {
            lo = ((totallen - d) / 2 - ky) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
            low = floor(lo);
            if (abs( round(lo)- (lo))<1/(long double)totallen){
                low = round(lo);
                cout<<"����"<<endl;
            }
            low = max(low, 1);
            hi = ((totallen + d) / 2 - ky-1) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
            high = ceil(hi);
            if (abs(round(hi)-(hi))<1/(long double)totallen){
                high = round(hi);
            }
            high = min(high, totallen);
            for (int ly = low-1; ly < high; ++ly) {
                        B[ly][lx] +=(complex<long double>){1.1111111,1.111111111111111};
                    }
                }
        int lx = 67;
        lo = ((totallen - d) / 2 - ky) / (totallen / 2 - kx-1) * (lx - kx) + ky+1;
        cout<<"���"<<((round(lo)-lo)<1/(long double)totallen)<<endl;
        cout<<"���"<<(1/(long double)totallen)<<endl;
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
        cout << j<<" ѭ�� "<<"The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    }
    clock_t end_time=clock();
    cout <<"��"<< "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    ofstream opt;
    opt.precision(20);
    opt.open("checker.csv");
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
}