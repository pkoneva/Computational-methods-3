#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



using namespace std;
ifstream fin;
ofstream fout;
const float a = 0.3;
double u(double x, double t) {
	return a * t + log(x - a * t + 2.2) - x + 0.4;
}
double u0(double t) {
	return u(0.0, t);
}



int main() {
	fout.open("output.txt");
	double t, h; // <1
	int N, M;
	cin >> h;
	cout << "enter a number less than " << h/a<<"\n";
	cin >> t;
	N = 1 + (1 / t);
	M = 1 + (1 / h);
	double coef = 0.5 * a * t / h;
	double* U = (double*)calloc(M, sizeof(double));
	double* U_next = (double*)calloc(M, sizeof(double)); 

	for (int i = 0; i < M; i++) {  // U(x, 0)
		U[i] = u(i * h, 0.0);
		fout << U[i] << " ";
	}
	fout << "\n";



	double S = 0;
	for (int j = 0; j < N - 1; j++) {
		double s;
		U_next[0] = u0(t * (j + 1.0));
		for (int k = 1; k < M - 1; k++) {
			U_next[k] = -coef * (U[k + 1] - U[k - 1]) + 0.5*(U[k + 1] + U[k - 1]);
			s = fabs(u(k * h, (j + 1.0) * t) - U_next[k]);
			S = (s > S ? s : S);
		}
		
		U_next[M - 1] = 3 * (U_next[M - 2] - U_next[M - 3]) + U_next[M - 4];
		s = fabs(u((M-1) * h, (j + 1.0) * t) - U_next[M - 1]);
		S = (s > S ? s : S);
		for (int k = 0; k < M; k++) {
			fout << U_next[k] << " ";
			U[k] = U_next[k];
		}
		fout << "\n";
	}
	fout << "Error of order "<<t+h<<"\nError S = " << S;
	free(U);
	free(U_next);
	fout.close();
	return 0;
}
