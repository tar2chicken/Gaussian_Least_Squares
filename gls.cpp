#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;


double fun(double a, double s, double m, double x);
double squares(double a, double s, double m, string iname, int n);
vector<double> LeastSquares(int d, double amin, double amax, double smin, double smax, double mmin, double mmax, string iname, int n);


int main() {
	// Get iofile information
	string ifname;
	int number;
	cout << "Enter the name of the input file." << endl;
	cin >> ifname;
	cout << "Enter the number of data." << endl;
	cin >> number;
	string ofname = ifname;
	ofname.insert(ofname.find("."), "_gls");

	// Determine parameter ranges and division number
	double amin, amax, smin, smax, mmin, mmax;
	int divi;
	cout << "Determine the range of normalization constant." << endl << "Enter min then max." << endl;
	cin >> amin >> amax;
	cout << "Determine the range of standard deviation." << endl << "Enter min then max." << endl;
	cin >> smin >> smax;
	cout << "Determine the range of mean." << endl << "Enter min then max." << endl;
	cin >> mmin >> mmax;
	cout << "Enter the division number." << endl;
	cin >> divi;

	// Find minimum
	vector<double> minimum(4);
	minimum = LeastSquares(divi, amin, amax, smin, smax, mmin, mmax, ifname, number);

	// Write data to ofile
	ifstream ifile(ifname);
	ofstream ofile(ofname);
	ofile << "a, " << minimum.at(1) << "+-" << (amax-amin)/divi << endl;
	ofile << "s, " << minimum.at(2) << "+-" << (smax-smin)/divi << endl;
	ofile << "m, " << minimum.at(3) << "+-" << (mmax-mmin)/divi << endl;
	ofile << ", , " << endl;
	ofile << ", raw, gls" << endl;
	for (int i = 0; i < number; i++) {
		double x, y;
		char comma;
		ifile >> x >> comma >> y;
		double z = fun(minimum.at(1), minimum.at(2), minimum.at(3), x);
		if (z < 1.0e-307) {
			z = 0;
		}
		ofile << x << ", " << y << ", " << z << endl;
	}
	ofile.close();

	cout << endl << "Done!" << endl;
	return 0;
}


double fun(double a, double s, double m, double x) {
	double t, y;
	t = -pow(((x-m)/s), 2) / 2.0;
	y = a*exp(t);
	return y;
}


double squares(double a, double s, double m, string iname, int n) {
	ifstream file(iname);
	double sum = 0;
	for (int i = 0; i < n; i++) {
		double x, y;
		char comma;
		file >> x >> comma >> y;
		double dif = y - fun(a, s, m, x);
		sum += pow(dif, 2);
	}
	return sum;
}


vector<double> LeastSquares(int d, double amin, double amax, double smin, double smax, double mmin, double mmax, string iname, int n) {
	// Calculate squares
	vector<vector<vector<vector<double>>>> val(d, vector<vector<vector<double>>>(d, vector<vector<double>>(d, vector<double>(4))));
	for (int i = 0; i < d; i++) {
		double a = amin + i*(amax-amin)/d;
		for (int j = 0; j < d; j++) {
			double s = smin + j*(smax-smin)/d;
			for (int k = 0; k < d; k++) {
				double m = mmin + k*(mmax-mmin)/d;
				val.at(i).at(j).at(k).at(0) = squares(a, s, m, iname, n);
				val.at(i).at(j).at(k).at(1) = a;
				val.at(i).at(j).at(k).at(2) = s;
				val.at(i).at(j).at(k).at(3) = m;
			}
		}
	}

	// Find minimum
	vector<vector<vector<double>>> valmin2(d, vector<vector<double>>(d, vector<double>(4)));
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			valmin2.at(i).at(j).at(0) = val.at(i).at(j).at(0).at(0);
			valmin2.at(i).at(j).at(1) = val.at(i).at(j).at(0).at(1);
			valmin2.at(i).at(j).at(2) = val.at(i).at(j).at(0).at(2);
			valmin2.at(i).at(j).at(3) = val.at(i).at(j).at(0).at(3);
			for (int k = 0; k < d-1; k++) {
				double x = valmin2.at(i).at(j).at(0);
				double y = val.at(i).at(j).at(k+1).at(0);
				if (x > y) {
					valmin2.at(i).at(j).at(0) = val.at(i).at(j).at(k+1).at(0);
					valmin2.at(i).at(j).at(1) = val.at(i).at(j).at(k+1).at(1);
					valmin2.at(i).at(j).at(2) = val.at(i).at(j).at(k+1).at(2);
					valmin2.at(i).at(j).at(3) = val.at(i).at(j).at(k+1).at(3);
				}
			}
		}
	}
	vector<vector<double>> valmin1(d, vector<double>(4));
	for (int i = 0; i < d; i++) {
		valmin1.at(i).at(0) = valmin2.at(i).at(0).at(0);
		valmin1.at(i).at(1) = valmin2.at(i).at(0).at(1);
		valmin1.at(i).at(2) = valmin2.at(i).at(0).at(2);
		valmin1.at(i).at(3) = valmin2.at(i).at(0).at(3);
		for (int j = 0; j < d-1; j++) {
			double x = valmin1.at(i).at(0);
			double y = valmin2.at(i).at(j+1).at(0);
			if (x > y) {
				valmin1.at(i).at(0) = valmin2.at(i).at(j+1).at(0);
				valmin1.at(i).at(1) = valmin2.at(i).at(j+1).at(1);
				valmin1.at(i).at(2) = valmin2.at(i).at(j+1).at(2);
				valmin1.at(i).at(3) = valmin2.at(i).at(j+1).at(3);
			}
		}
	}
	vector<double> valmin(4);
	valmin.at(0) = valmin1.at(0).at(0);
	valmin.at(1) = valmin1.at(0).at(1);
	valmin.at(2) = valmin1.at(0).at(2);
	valmin.at(3) = valmin1.at(0).at(3);
	for (int i = 0; i < d-1; i++) {
		double x = valmin.at(0);
		double y = valmin1.at(i+1).at(0);
		if (x > y) {
			valmin.at(0) = valmin1.at(i+1).at(0);
			valmin.at(1) = valmin1.at(i+1).at(1);
			valmin.at(2) = valmin1.at(i+1).at(2);
			valmin.at(3) = valmin1.at(i+1).at(3);
		}
	}

	return valmin;
}

