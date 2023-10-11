#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



using namespace std;
ifstream fin;
ofstream fout;

double u(double x, double t) {
	return 2 * pow(x, 4) - 3 * pow(t, 3) + 3 * pow(t, 2) * x - 2 * exp(x);
}
double u0(double t) {
	return u(0.0, t);
}
double u1(double t) {
	return u(1.0, t);
}
double f(double x, double t) {
	return -pow(t, 2)*9 + 6 * t * x - 0.027 * pow(x, 2) *24 + 2*0.027 * exp(x);
}


int main() {
	fin.open("input.txt");
	fout.open("output.txt");
	double t, h; // <1
	int N, M;
	cin >> t;
	cout << "enter a number greater than" << sqrt(2 * 0.027 * t);
	cin >> h;
	N = 1+(1 / t);
	M = 1+(1 / h);
    double coef = 0.027*t/pow(h, 2); 
	double *U = (double*)calloc(M, sizeof(double));
	double *U_next = (double*)calloc(M, sizeof(double)); // Следующий слой
    
    for (int i = 0; i < M; i++) {  // U(x, 0)
		U[i] = u(i * h, 0.0);
		fout << U[i] <<" ";
	}
	fout << "\n";
	


    double S = 0;
	for (int j = 0; j < N-1; j++) {
		U_next[0] = u0(t*(j+1.0));
		U_next[M-1] = u1(t*(j+1.0));
		for (int k = 1; k < M-1; k++) {
			U_next[k] = coef * (U[k + 1] + U[k - 1] - 2 * U[k]) + U[k] + t*f(k*h, j*t);
			double s = fabs(u(k*h, (j+1.0)*t) - U_next[k]);
			S = (s > S ? s :S);
		}
		for (int k = 0; k < M; k++) {
			fout << U_next[k] << " ";
			U[k] = U_next[k];
		}
		fout << "\n";
	}
	fout << "Ïîãðåøíîñòü S = " << S;
	free(U);
	free(U_next);
	fin.close();
	fout.close();
	return 0;
}
