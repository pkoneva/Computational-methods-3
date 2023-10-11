#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
const double a = -1.1;
const double b = -0.8;

ifstream fin;
ofstream fout;
double next_max(double y01, double y11, double y21, double y10, double y12) {
	return a * (y01 - 2 *y11 + y21) + b * (y10 - 2 * y11 + y12);
}

double next_min( double y01, double y11, double y21, double y10, double y12, double lmax, double h) {
	return lmax* y11 - (a * (y01 - 2 * y11 + y21) + b * (y10 - 2 * y11 + y12))*pow(h, -2);
}
void normY(double** Y, int M) {
	double S=0;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) S += pow(Y[i][j], 2);
	}
	S = sqrt(S);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) Y[i][j] /= S;
	}

}

double eigenvalue(double** Y0, double** Y1, int M) {
	double S = 0;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) 
			S += Y0[i][j]* Y1[i][j];
	}
	return S;
}

int main() {
	fout.open("output.txt");
	fin.open("input.txt");
	double h, e;
	fin >> h >> e;
	int M = int(1 / h)/*-1*/;
	int n = 0;
	double** Y = new double* [M+1];
	double** Y_next = new double* [M+1];
	for (int i = 0; i < M+1; i++) {
		Y[i] = new double[M+1];
		Y_next[i] = new double[M+1];
		for (int j = 0; j < M+1; j++) Y[i][j] = 0;
	}
	for (int j = 1; j < M / 2; j++) {
		for (int i = 1; i < M; i++)
			Y[i][j] = 1;
	}
	for (int j = M / 2; j < M; j++) {
		for (int i = 1 + M / 2; i < M; i++) {
			Y[i][j] = 1;
		}
	}
	double L = 0, Lnext = 0;
	do {
		n++;
		normY(Y, M+1);

		for (int i = 0; i < M+1; i++) {
			for (int j = 0; j < M+1; j++)
				Y_next[i][j] = Y[i][j];
			
		}
		for (int j = 1; j < M / 2; j++) {
			for (int i = 1; i < M; i++)
				Y_next[i][j] = next_max(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1]);
		}
		for (int j = M / 2; j < M; j++) {
			for (int i = 1 + M / 2; i < M; i++) {
				Y_next[i][j] = next_max(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1]);
			}
		}
		
		L = Lnext;
		Lnext = eigenvalue(Y, Y_next, M+1) / pow(h, 2);
		for (int i = 0; i < M+1; i++)
			for (int j = 0; j < M+1; j++)
				Y[i][j] = Y_next[i][j];
	} while ((fabs(L - Lnext) / L) > e);
	const double lmax = Lnext;
	
        //поиск наименьшего собственного числа
	fout <<n << ";" << lmax << "\n\n";
	for (int j = 1; j < M / 2; j++) {
		for (int i = 1; i < M; i++)
			Y[i][j] = 1;
	}
	for (int j = M / 2; j < M; j++) {
		for (int i = 1 + M / 2; i < M; i++) {
			Y[i][j] = 1;
		}
	}
	L = 0; Lnext = 0;
	n = 0;
	do{
		n++;
		normY(Y, M+1);

		for (int i = 0; i < M+1; i++) {
			for (int j = 0; j < M+1; j++)
				Y_next[i][j] = Y[i][j];
			
		}
		for (int j = 1; j < M / 2; j++) {
			for (int i = 1; i < M; i++)
				Y_next[i][j] = next_min(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], lmax, h);
		}
		for (int j = M / 2; j < M; j++) {
			for (int i = 1 + M / 2; i < M; i++) {
				Y_next[i][j] = next_min(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], lmax, h);
			}
		}

		L = Lnext;
		Lnext = eigenvalue(Y, Y_next, M+1);
		for (int i = 0; i < M+1; i++)
			for (int j = 0; j < M+1; j++)
				Y[i][j] = Y_next[i][j];
		
	} while (fabs((L - Lnext) / (lmax - L)) > e);


	fout << n << ";" << lmax - Lnext;
	fout.close();
	fin.close();
	delete[] Y;
	delete[] Y_next;
	return 0;
}
