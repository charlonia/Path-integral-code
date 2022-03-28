//
// Created by charlonia on 2022/3/20.
//
#include <iostream>
#include <complex.h>
#include <math.h>
#include <array>
using namespace std;
int main(){
    clock_t start_time=clock();

    char  watoin;
    cin>>watoin;
    double grid =1;
    double grid_div=0.01;
    double time_div = 0.002; //ʱ����
    double time_all = 0.004;     //�ܵ��ݻ�ʱ��
    double d = 0.4/grid_div;          //������
    double v=10; // ���õ��������ٶ�
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
    complex<double> A[totallen][totallen];
    complex<double> B[totallen][totallen];
    const double z=-pow(grid_div,2)/(2*time_div)*m/hb ;
    for (int j=0;j<totaltime;++j){
        for (int kx=0;kx < round(totallen/2)-1;++kx) {
            for (int ky = 0; ky < totallen; ++ky) {
                for (int lx = 0; lx < round(totallen / 2) - 1; ++lx) {
                    for (int ly = 0; ly < totallen; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<double>){pow(grid_div,2),0};
                    }
                }
                for (int lx = round(totallen / 2)-1; lx < totallen - 1; ++lx) {
                    low = round(((totallen - d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    low = max(low, 0);
                    high = round(((totallen + d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    high = min(high, totallen);
                    for (int ly = low-1; ly < high; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<double>){pow(grid_div,2),0};
                    }
                }
            }
        }
        for (int ky=round((totallen - d) / 2);ky<round((totallen + d) / 2);++ky){
            for(int lx = 0;lx < totallen - 1;++lx){
                for(int ly = 0;ly < totallen ;++ly){
                    B[ky][midx] +=A[ly][lx] * exp((complex<double>){0, z * (pow(lx - midx, 2) + pow(ky - ly, 2))})*(complex<double>){pow(grid_div,2),0};

                }
            }
        }
        for (int kx=round(totallen/2);kx < totallen;++kx){
            for(int ky=0;ky < totallen;++ky){
                for (int lx = round(totallen/2); lx < totallen - 1; ++lx) {
                    for (int ly = 0; ly < totallen; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<double>){0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<double>){pow(grid_div,2),0};
                    }
                }
                for (int lx = 0; lx <round(totallen / 2); ++lx) {
                    low = round(((totallen - d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    low = max(low, 0);
                    high = round(((totallen + d) / 2 - ky) / (totallen / 2 - kx) * (lx - kx) + ky);
                    high = min(high, totallen);
                    for (int ly = low; ly < high; ++ly) {
                        B[ky][kx] +=A[ly][lx] * exp((complex<double>) {0, z * (pow(lx - kx, 2) + pow(ky - ly, 2))})*(complex<double>){pow(grid_div,2),0};
                    }
                }
            }
        }
        for (int y =0;totallen;++y){
            for(int x = 0; totallen;++x){
                A[y][x]=B[y][x];
            }
        }
        cout<<j+1<<"finished"<<endl;
    }
    clock_t end_time=clock();
    cout << "The run time is: " <<(double)(end_time - start_time) / CLOCKS_PER_SEC <<endl;
    return 0;
}
