#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
using namespace std;
ifstream fin;
ofstream fout;
const double a = -1.1;
const double b = -0.8;
void norm(double** Y, int M) {
	double S = 0;
	for (int i = 0; i <= M; i++) {
		for (int j = 0; j <= M; j++) S += pow(Y[i][j], 2);
	}
	S = sqrt(S) / (M + 2);
	for (int i = 0; i <= M; i++) {
		for (int j = 0; j <= M; j++) Y[i][j] /= S;
	}

}
double scalar_mult(double** y1, double** y2, int M) {
	double S = 0;
	
	for (int j = 0; j <= M; j++) {
		for (int i = 0; i <= M; i++)
			S += y1[i][j] * y2[i][j];
	}
	
	return S;
}
double f(double x, double y) {
	return 1.1 * sin(x) + (3.2 * pow(x, 2) + 4.4 * pow(y, 2)) * cos(2 * x * y);
}
double fi(double x, double y) {
	return sin(x) + cos(2 * x * y);
}
double A(double y01, double y11, double y21, double y10, double y12, double h) {
	return (a * (y01 - 2 * y11 + y21) + b * (y10 - 2 * y11 + y12)) ;
}


void AY(double** y, double** yn, int M) {
	double h = pow(M, -1);
	
	yn[0][0] = A(fi(0, h), y[0][0], y[1][0], fi(h, 0), y[0][1], h);
	for (int j = 1; j < M; j++) yn[j][0] = A(y[j - 1][0], y[j][0], y[j][1], fi(h * (j + 1), 0), y[j + 1][0], h);
	yn[M][0] = A(y[M - 1][0], y[M][0], fi(1, h), y[M][1], fi(h * (M + 1), 0), h);

	
	for (int i = 1; i <M; i++) {
		yn[0][i] = A(fi(0, h * (i + 1)), y[0][i], y[1][i], y[0][i + 1], y[0][i - 1], h);
		for (int j = 1; j < M; j++) {
			yn[j][i] = A(y[j - 1][i], y[j][i], y[j + 1][i], y[j][i + 1], y[j][i - 1], h);
		}
		yn[M][i] = A(y[M - 1][i], y[M][i], fi(1, h * (i + 1)), y[M][i + 1], y[M][i - 1], h);
	}


	yn[0][M] = A(fi(0, h * (M + 1)), y[1 + M / 2][M], y[2 + M / 2][M], y[1 + M / 2][M - 1], fi(h, 1), h);
	for (int j = 1; j < M; j++) yn[j][M] = A(y[j - 1][M], y[j][M], y[j + 1][M], y[j][M - 1], fi(h * (j + 1), 1), h);
	yn[M][M] = A(y[M - 1][M], y[M][M], fi(1, h * (M + 1)), y[M][M - 1], fi(h * (M + 1), 1), h);

}
double coef(double** err, double** err_next, int M) {
	return scalar_mult(err, err, M) / scalar_mult(err, err_next, M);
}
int main() {

	fin.open("input.txt");
	fout.open("output.txt");
	double h, e;
	fin >> h >> e;
	int M = int(1 / h) - 2;
	int n = 0;
	//задание всех матриц
	double** yn = new double* [M + 1];
	double** y = new double* [M + 1];
	double** Ay = new double* [M + 1];
	double** err = new double* [M + 1];
	double** Aerr = new double* [M + 1];
	double** F = new double* [M + 1];
	
	
	for (int i = 0; i <= M; i++) {
		y[i] = new double[M + 1]; yn[i] = new double[M + 1];
		Ay[i] = new double[M + 1];
		err[i] = new double[M + 1];
		Aerr[i] = new double[M + 1];
		for (int j = 0; j <= M; j++) {
			err[i][j] = 0;
			Aerr[i][j] = 0;
		}
		F[i] = new double[M];
	}
	for (int j = 0; j <= M; j++) {
		for (int i = 0; i <= M; i++) {
			y[i][j] = 1;
			F[i][j] = pow(h,2)*f((i + 1) * h, (j + 1) * h);
		}

	}
	


	double t;
	double dif = e + 1, norma = 0, nor = 0;
	while (dif > e) {
		n++;
		norm(y, M);
		AY(y, Ay, M);
		for (int j = 0; j <= M; j++) {
			for (int i = 0; i <= M; i++)
				err[i][j] = Ay[i][j] - F[i][j];
		}
		
		AY(err, Aerr, M);
		t = scalar_mult(err, err, M) / scalar_mult(Aerr, err, M);

		for (int j = 0; j <=M; j++) {
			for (int i = 0; i <= M; i++) {
				yn[i][j] = y[i][j] - t * err[i][j];
				
			}
		}
		for (int j = 0; j <= M; j++) {
			for (int i = 0; i <= M; i++) {
				err[i][j] = yn[i][j] - y[i][j];
				nor += pow(/*t * */err[i][j], 2);
				y[i][j] = yn[i][j];
				
			}
		}
		dif = h * fabs(t) * sqrt(nor);
	};

	for (int j = 0; j <= M; j++) {
		for (int i = 0; i <= M; i++)
			F[i][j] = fi((i + 1) * h, (j + 1) * h);;
	}
	for (int i = 0; i <= M; i++)
		for (int j = 0; j <= M; j++) {
			Ay[i][j] = y[i][j] - F[i][j];
		}
	
	fout << n << "; " << h * sqrt(scalar_mult(Ay, Ay, M));
	fin.close();
	fout.close();
	delete[] y;
	delete[] F;
	delete[] Ay;
	delete[] err;
	delete[] Aerr;

	return 0;
}
