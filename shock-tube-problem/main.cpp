#include <bits/stdc++.h>
#include<iostream>
#include<vector>
#include "func.h"
using namespace std;

int main() {
    // 定数系定義
    double x_min = -5.0;
    double x_max = 5.0;
    double dx = 0.1;
    double dt = 0.1;
    // 格子数
    int nx = (x_max-x_min)/dx;
    // ステップ数
    int nt = 100;
    // 初期条件
    double pl = 1.0;
    double pr = 0.1;
    double rhol = 1.0;
    double rhor = 0.125;
    double Ul = 0.0;
    double Ur = 0.0;
    // 配列系定義
    double x[nx];
    double Q[nx][3];
    double Q_kome[nx][3];
    double E[nx][3];
    double E_kome[nx][3];
    for(int i=0; i<nx; i++){
        x[i] = x_min;
        x_min += dx;
        q=Q[i];
        if(x[i]<=0){
        *q=rhol;
        *(q+1)=Ul;
        *(q+2) = rhol;
        }
        cout << Q[i] << endl;
    }

	return 0;
}

