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

//double** nextY(double** Y, int M, double A, double B) {
//	double** Y_next = new double* [M];
//	for (int i = 0; i < M; i++) 
//		Y_next[i] = new double[M];
//	Y_next[0][0] = next(0, Y[0][0], Y[1][0], 0, Y[0][1], A, B);
//	for (int j = 1; j < M - 1; j++) Y_next[0][j] = next(0, Y[0][j], Y[1][j], Y[0][j - 1], Y[0][j + 1], A, B);
//	Y_next[0][M - 1] = next(0, Y[0][M - 1], Y[1][M - 1], Y[0][M - 2], 0, A, B);
//	for (int i = 1; i < M - 1; i++) {
//		if (i == 0) {
//			for (int j = 0; j < M; j++) {
//				Y_next[i][j] = A * (Y[i - 1][j] - 2 * Y[i][j] + Y[i + 1][j]) + B * (Y[i][j - 1] - 2 * Y[i][j] + Y[i][j + 1]);
//			}
//		}
//		Y_next[i][0] = next(Y[i - 1][0], Y[i][0], Y[i + 1][0], 0, Y[i][1], A, B);
//		for (int j = 1; j < M - 1; j++) {
//			/*Y_next[i][j] = A * (Y[i - 1][j] - 2 * Y[i][j] + Y[i + 1][j]) + B * (Y[i][j - 1] - 2 * Y[i][j] + Y[i][j + 1]);*/
//			Y_next[i][j] = next(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], A, B);
//		}
//		Y_next[i][M - 1] = next(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0, A, B);
//	}
//	Y_next[M - 1][0] = next(Y[M - 2][0], Y[M - 1][0], 0, 0, Y[M - 1][1], A, B);
//	for (int j = 1; j < M - 1; j++) Y_next[M - 1][j] = next(Y[M - 2][j], Y[M - 1][j], 0, Y[M - 1][j - 1], Y[M - 1][j + 1], A, B);
//	Y_next[M - 1][M - 1] = next(Y[M - 2][M - 1], Y[M - 1][M - 1], 0, Y[M - 1][M - 2], 0, A, B);
//	return Y_next;
//	
//}
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
		//Y_next[1][1] = next_max(0, Y[1][1], Y[2][1], 0, Y[1][2]);
		//for (int j = 1; j < M - 1; j++) Y_next[0][j] = next_max(0, Y[0][j], Y[1][j], Y[0][j - 1], Y[0][j + 1]);
		//Y_next[0][M - 1] = next_max(0, Y[0][M - 1], Y[1][M - 1], Y[0][M - 2], 0);
		//for (int i = 1; i < (M/* - 1*/) / 2; i++) {
		//	Y_next[i][0] = next_max(Y[i - 1][0], Y[i][0], Y[i + 1][0], 0, Y[i][1]);
		//	for (int j = 1; j < M - 1; j++) {
		//		Y_next[i][j] = next_max(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1]);
		//	}
		//	Y_next[i][M - 1] = next_max(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0);
		//}
		//Y_next[(M - 1) / 2][0] = next_max(Y[(M - 3) / 2][0], Y[(M - 1) / 2][0], 0, 0, Y[(M - 1) / 2][1]);
		//for (int j = 1; j < (M + 3) / 2; j++) Y_next[(M - 1) / 2][j] = next_max(Y[(M - 3) / 2][j], Y[(M - 1) / 2][j], 0, Y[(M - 1) / 2][j - 1], Y[(M - 1) / 2][j + 1]);
//for (int i = (M + 1) / 2; i < M - 1; i++) {
		//	Y_next[i][(M + 3) / 2] = next_max(Y[i - 1][(M + 3) / 2], Y[i][(M + 3) / 2], Y[i + 1][(M + 3) / 2], 0, Y[i][(M + 5) / 2]);
		//	for (int j = (M + 5) / 2; j < M - 1; j++) {
//		Y_next[i][j] = next_max(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1]);
		//	}
		//	Y_next[i][M - 1] = next_max(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0);
		//}
//Y_next[M - 1][(M + 3) / 2] = next_max(Y[M - 2][(M + 3) / 2], Y[M - 1][(M + 3) / 2], 0, 0, Y[M - 1][(M + 5) / 2]);
		//for (int j = (M + 5) / 2; j < M - 1; j++) Y_next[M - 1][j] = next_max(Y[M - 2][j], Y[M - 1][j], 0, Y[M - 1][j - 1], Y[M - 1][j + 1]);
		//Y_next[M - 1][M - 1] = next_max(Y[M - 2][M - 1], Y[M - 1][M - 1], 0, Y[M - 1][M - 2], 0);
//for (int i = (M + 1) / 2; i < M - 1; i++) {
		//	for (int j = 0; j < (M + 1) / 2; j++) Y_next[i][j] = Y[i][j];
		//}




		//Y_next[0][0] = next_max(0, Y[0][0], Y[1][0], 0, Y[0][1]);
		//	for (int j = 1; j < M - 1; j++) Y_next[0][j] = next_max(0, Y[0][j], Y[1][j], Y[0][j - 1], Y[0][j + 1]);
		//	Y_next[0][M - 1] = next_max(0, Y[0][M - 1], Y[1][M - 1], Y[0][M - 2], 0);
		//	for (int i = 1; i < M - 1; i++) {
		//		
		//		Y_next[i][0] = next_max(Y[i - 1][0], Y[i][0], Y[i + 1][0], 0, Y[i][1]);
		//		for (int j = 1; j < M - 1; j++) {
		//			/*Y_next[i][j] = A * (Y[i - 1][j] - 2 * Y[i][j] + Y[i + 1][j]) + B * (Y[i][j - 1] - 2 * Y[i][j] + Y[i][j + 1]);*/
		//			Y_next[i][j] = next_max(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1]);
		//		}
		//		Y_next[i][M - 1] = next_max(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0);
		//	}
		//	Y_next[M - 1][0] = next_max(Y[M - 2][0], Y[M - 1][0], 0, 0, Y[M - 1][1]);
		//	for (int j = 1; j < M - 1; j++) Y_next[M - 1][j] = next_max(Y[M - 2][j], Y[M - 1][j], 0, Y[M - 1][j - 1], Y[M - 1][j + 1]);
		//	Y_next[M - 1][M - 1] = next_max(Y[M - 2][M - 1], Y[M - 1][M - 1], 0, Y[M - 1][M - 2], 0);
		
		
		L = Lnext;
		Lnext = eigenvalue(Y, Y_next, M+1) / pow(h, 2);
		for (int i = 0; i < M+1; i++)
			for (int j = 0; j < M+1; j++)
				Y[i][j] = Y_next[i][j];
	} while ((fabs(L - Lnext) / L) > e);
	const double lmax = Lnext;

	fout <<n << ";" << lmax << "\n\n";

	//поиск наименьшего собственного числа
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

		/*Y_next[0][0] = next_min(0, Y[0][0], Y[1][0], 0, Y[0][1], lmax, h);
		for (int j = 1; j < M - 1; j++) Y_next[0][j] = next_min(0, Y[0][j], Y[1][j], Y[0][j - 1], Y[0][j + 1], lmax, h);
		Y_next[0][M - 1] = next_min(0, Y[0][M - 1], Y[1][M - 1], Y[0][M - 2], 0, lmax, h);
		for (int i = 1; i < (M - 1) / 2; i++) {
			Y_next[i][0] = next_min(Y[i - 1][0], Y[i][0], Y[i + 1][0], 0, Y[i][1], lmax, h);
			for (int j = 1; j < M - 1; j++) {
				Y_next[i][j] = next_min(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], lmax, h);
			}
			Y_next[i][M - 1] = next_min(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0, lmax, h);
		}
		Y_next[(M - 1) / 2][0] = next_min(Y[(M - 3) / 2][0], Y[(M - 1) / 2][0], 0, 0, Y[(M - 1) / 2][1], lmax, h);
		for (int j = 1; j < (M + 3) / 2; j++) Y_next[(M - 1) / 2][j] = next_min(Y[(M - 3) / 2][j], Y[(M - 1) / 2][j], 0, Y[(M - 1) / 2][j - 1], Y[(M - 1) / 2][j + 1], lmax, h);

		for (int i = (M + 1) / 2; i < M - 1; i++) {
			Y_next[i][(M + 3) / 2] = next_min(Y[i - 1][(M + 3) / 2], Y[i][(M + 3) / 2], Y[i + 1][(M + 3) / 2], 0, Y[i][(M + 5) / 2], lmax, h);
			for (int j = (M + 5) / 2; j < M - 1; j++) {

				Y_next[i][j] = next_min(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], lmax, h);
			}
			Y_next[i][M - 1] = next_min(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0, lmax, h);
		}



		Y_next[M - 1][(M + 3) / 2] = next_min(Y[M - 2][(M + 3) / 2], Y[M - 1][(M + 3) / 2], 0, 0, Y[M - 1][(M + 5) / 2], lmax, h);
		for (int j = (M + 5) / 2; j < M - 1; j++) Y_next[M - 1][j] = next_min(Y[M - 2][j], Y[M - 1][j], 0, Y[M - 1][j - 1], Y[M - 1][j + 1], lmax, h);
		Y_next[M - 1][M - 1] = next_min(Y[M - 2][M - 1], Y[M - 1][M - 1], 0, Y[M - 1][M - 2], 0, lmax, h);

		for (int i = (M + 1) / 2; i < M - 1; i++) {
			for (int j = 0; j < (M + 1) / 2; j++) Y_next[i][j] = Y[i][j];
		}

*/

		//Y_next[0][0] = next_min(0, Y[0][0], Y[1][0], 0, Y[0][1], lmax, h);
		//for (int j = 1; j < M - 1; j++) Y_next[0][j] = next_min(0, Y[0][j], Y[1][j], Y[0][j - 1], Y[0][j + 1], lmax, h);
		//Y_next[0][M - 1] = next_min(0, Y[0][M - 1], Y[1][M - 1], Y[0][M - 2], 0,lmax, h);
		//for (int i = 1; i < M - 1; i++) {
		//	Y_next[i][0] = next_min(Y[i - 1][0], Y[i][0], Y[i + 1][0], 0, Y[i][1],lmax, h);
		//	for (int j = 1; j < M - 1; j++) {
		//		/*Y_next[i][j] = A * (Y[i - 1][j] - 2 * Y[i][j] + Y[i + 1][j]) + B * (Y[i][j - 1] - 2 * Y[i][j] + Y[i][j + 1]);*/
		//		Y_next[i][j] = next_min(Y[i - 1][j], Y[i][j], Y[i + 1][j], Y[i][j - 1], Y[i][j + 1], lmax, h);
		//	}
		//	Y_next[i][M - 1] = next_min(Y[i - 1][M - 1], Y[i][M - 1], Y[i + 1][M - 1], Y[i][M - 2], 0,lmax, h);
		//}
		//Y_next[M - 1][0] = next_min(Y[M - 2][0], Y[M - 1][0], 0, 0, Y[M - 1][1],lmax, h);
		//for (int j = 1; j < M - 1; j++) 
		//	Y_next[M - 1][j] = next_min(Y[M - 2][j], Y[M - 1][j], 0, Y[M - 1][j - 1], Y[M - 1][j + 1],lmax, h);
		//Y_next[M - 1][M - 1] = next_min(Y[M - 2][M - 1], Y[M - 1][M - 1], 0, Y[M - 1][M - 2], 0,lmax, h);


		L = Lnext;
		Lnext = eigenvalue(Y, Y_next, M+1) /*/ pow(h, 2)*/;
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