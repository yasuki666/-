
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
using namespace std;


int rho = 1;
int Re =1000 ;//雷诺数
double L = 1;//方腔尺寸
const int n = 200;
double dh = L / n;//detla h
double dt = 1e-3;//时间步长
double xi[n + 1][n + 1] = { 0 };
double psi[n + 1][n + 1] = { 0 };
double u[n + 1][n + 1], v[n + 1][n + 1];
double err1, err2, temp;
double err;
int i, j;
double err1_line[50000], err2_line[50000];
double max(double a, double b) {
	if (a > b)
		return a;
	else
		return b;
}
void output1(int m)	//输出
{
	ostringstream name;
	name << "cavity_" << m << ".dat";
	ofstream out(name.str().c_str());
	out << "Title=\" Lid Driven Flow\"\n"
		<< "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"psi\",\"xi\"\n"
		<< "ZONE T= \"BOX\", I= "
		<< n + 1 << ", J=" << n + 1 << ", F=POINT" << endl;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			out << double(i) / n << " " << double(j) / n << " "
				<< u[i][j] << " " << v[i][j] << " " << psi[i][j] << " " << xi[i][j] << endl;
		}
}

void output2(double a[50000],double b[50000])	//输出
{
	ostringstream name;
	name << "err"  << ".dat";
	ofstream out(name.str().c_str());
	out << "Title=\" Lid Driven Flow\"\n"
		<< "VARIABLES=\"time\",\"err1\",\"err2\"\n"
		<< "ZONE T= \"BOX\", I= "
		<< n + 1 << ", J=" << n + 1 << ", F=POINT" << endl;
	
	for (i = 0; i <= 50000; i++)
	{
		out << double(i)/50 << " " << a[i]  << " " << b[i] << endl;
	}
}
int main() {
	for (int k = 1; k < 50000; k++) {
		err = 0;
		//边界条件
		for (int i = 2; i <= n; i++) {
			xi[i][1] = -2 * (psi[i][2] /*- psi[i][1]*/) / (dh * dh);
			xi[i][n + 1] = -2 * (psi[i][n] /*- psi[i][n + 1]*/ + dh) / (dh * dh);
		}
		for (int j = 2; j <= n; j++) {
			xi[1][j] = -2 * (psi[2][j] /*- psi[1][j]*/) / (dh * dh);
			xi[n + 1][j] = -2 * (psi[n][j] /*- psi[n + 1][j]*/) / (dh * dh);
		}
		//控制方程
		for (i = 2; i <= n; i++) {
			for (j = 2; j <= n; j++) {
				u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * dh);
				v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * dh);
				err1 = (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] + xi[i][j] * dh * dh) / 4 - psi[i][j];
				psi[i][j] = psi[i][j] + rho * err1;
				err2 = dt * (-dh / 2 * (u[i][j] * (xi[i + 1][j] - xi[i - 1][j]) + v[i][j] * (xi[i][j + 1] - xi[i][j - 1])) + (xi[i + 1][j] + xi[i - 1][j] + xi[i][j + 1] + xi[i][j - 1] - 4 * xi[i][j]) / Re) / (dh * dh);
				xi[i][j] = xi[i][j] + rho * err2;
				temp = max(abs(err1), abs(err2));
				if (err < temp)
					err = temp;
			}
		}
		if (k % 1000 == 0) {
			cout << "k=" << k << "  " << "err=" << err << "   " << err1 << "   " << err2 << endl;
			output1(k);
		}
		err1_line[k - 1] = err1;
		err2_line[k - 1] = err2;
		if (err < 1e-6)
			break;
	}

	output2(err1_line,err2_line);
}