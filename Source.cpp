#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>

using namespace std;

void solve_1(vector<vector<double>>& list1, vector<vector<double>>& list2, vector<vector<double>>& list3, vector<vector<double>>& list4,
	double dx, double dt, double a) {
	vector<vector<double>> temp_left(list1);
	vector<vector<double>> temp_right(list2);
#pragma omp parallel for
	for (int i = 1; i < list1.size() - 1; ++i) {
		for (int j = 1; j < list1[i].size() - 1; ++j) {
			temp_left[i][j] += a * dt * ((list1[i][j + 1] - 2 * list1[i][j] + list1[i][j - 1]) / dx + (list1[i + 1][j] - 2 * list1[i][j] + list1[i - 1][j]) / dx);
		}
	}
#pragma omp parallel for
	for (int i = 1; i < list2.size() - 1; ++i) {
		for (int j = 1; j < list2[i].size() - 1; ++j) {
			temp_right[i][j] += a * dt * ((list2[i][j + 1] - 2 * list2[i][j] + list2[i][j - 1]) / dx + (list2[i + 1][j] - 2 * list2[i][j] + list2[i - 1][j]) / dx);
		}
	}

	vector<double> temp_mid(list1.size(), 0);

	temp_mid[0] = temp_left[0][temp_left[0].size() - 1];
	temp_mid[temp_mid.size() - 1] = temp_left[temp_left.size() - 1][temp_left[temp_left.size() - 1].size() - 1];

	for (int i = 1; i < list1.size() - 1; ++i) {
		temp_mid[i] = list1[i][list1[i].size() - 1] + a * dt *
			((list1[i - 1][list1[i - 1].size() - 1] - 2 * list1[i][list1[i].size() - 1] + list1[i + 1][list1[i + 1].size() - 1]) / dx
				+ (list1[i][list1[i].size() - 2] - 2 * list1[i][list1[i].size() - 1] + list2[i][1]) / dx);
	}

	list1 = temp_left;
	temp_left = temp_right;
	temp_right = list3;
	for (int i = 0; i < list1.size(); ++i) {
		list1[i][list1[i].size() - 1] = temp_mid[i];
	}
	temp_mid.resize(list3.size(), 0);

#pragma omp parallel for
	for (int i = 1; i < list3.size() - 1; ++i) {
		for (int j = 1; j < list3[i].size() - 1; ++j) {
			temp_right[i][j] += a * dt * ((list3[i][j + 1] - 2 * list3[i][j] + list3[i][j - 1]) / dx + (list3[i + 1][j] - 2 * list3[i][j] + list3[i - 1][j]) / dx);
		}
	}

#pragma omp parallel for
	for (int i = 1; i < list3.size() - 1; ++i) {
		temp_mid[i] = list3[i][list3[i].size() - 1] + a * dt *
			((list3[i - 1][list3[i - 1].size() - 1] - 2 * list3[i][list3[i].size() - 1] + list3[i + 1][list3[i + 1].size() - 1]) / dx
				+ (list3[i][list3[i].size() - 2] - 2 * list3[i][list3[i].size() - 1] + list2[list2.size() - list3.size() + i][1]) / dx);
	}
	list3 = temp_right;
	list2 = temp_left;
	temp_right = temp_left;
	temp_left = list4;

	for (int i = 1; i < list1.size() - 1; ++i) {
		list2[i][0] = list1[i][list1[i].size() - 1];
	}

	for (int i = 1; i < list3.size(); ++i) {
		list3[i][list3[i].size() - 1] = temp_mid[i];
		list2[list2.size() - list3.size() + i][0] = temp_mid[i];
	}
	temp_mid.clear();
	temp_mid.resize(list4[0].size());

	for (int i = 1; i < temp_left.size() - 1; ++i) {
		for (int j = 1; j < temp_left[i].size() - 1; ++j) {
			temp_left[i][j] += a * dt * ((list4[i][j + 1] - 2 * list4[i][j] + list4[i][j - 1]) / dx + (list4[i + 1][j] - 2 * list4[i][j] + list4[i - 1][j]) / dx);
		}
	}

	for (int j = 1; j < temp_mid.size() - 1; ++j) {
		temp_mid[j] = temp_right[temp_right.size() - 1][j] +
			a * dt * ((temp_right[temp_right.size() - 2][j] - 2 * temp_right[temp_right.size() - 1][j] + temp_left[1][j]) / dx
				+ (temp_right[temp_right.size() - 1][j - 1] - 2 * temp_right[temp_right.size() - 1][j] + temp_right[temp_right.size() - 1][j + 1]) / dx);
	}

	list4 = temp_left;

	for (int j = 0; j < temp_mid.size(); ++j) {
		list4[0][j] = temp_mid[j];
		list2[list2.size() - 1][j] = temp_mid[j];
	}
}

void solve_2(vector<vector<double>>& list1, vector<vector<double>>& list2, vector<vector<double>>& list3, vector<vector<double>>& list4, double dx, double dt, double a) {
	vector<vector<double>> temp_left(list1);
	vector<vector<double>> temp_right(list2);
#pragma omp parallel for
	for (int i = 1; i < list1.size() - 1; ++i) {
		for (int j = 1; j < list1[i].size() - 1; ++j) {
			temp_left[i][j] += a * dt * ((list1[i][j + 1] - 2 * list1[i][j] + list1[i][j - 1]) / dx + (list1[i + 1][j] - 2 * list1[i][j] + list1[i - 1][j]) / dx);
		}
	}

	for (int j = 0; j < list1[0].size(); ++j) {
		temp_left[0][j] = temp_left[1][j];
		temp_left[temp_left.size() - 1][j] = temp_left[temp_left.size() - 2][j];
	}

	for (int i = 0; i < temp_left.size(); ++i) {
		temp_left[i][0] = temp_left[i][1];
	}

#pragma omp parallel for
	for (int i = 1; i < list2.size() - 1; ++i) {
		for (int j = 1; j < list2[i].size() - 1; ++j) {
			temp_right[i][j] += a * dt * ((list2[i][j + 1] - 2 * list2[i][j] + list2[i][j - 1]) / dx + (list2[i + 1][j] - 2 * list2[i][j] + list2[i - 1][j]) / dx);
		}
	}

	for (int i = 0; i < list2.size(); ++i) {
		temp_right[i][0] = temp_right[i][1];
		temp_right[i][temp_right[i].size() - 1] = temp_right[i][temp_right[i].size() - 2];
	}

	vector<double> temp_mid(list1.size(), 0);
	temp_mid[0] = temp_left[0][temp_left[0].size() - 1];
	temp_mid[temp_mid.size() - 1] = temp_left[temp_mid.size() - 1][temp_left[temp_mid.size() - 1].size() - 1];

	for (int i = 1; i < list1.size() - 1; ++i) {
		temp_mid[i] = list1[i][list1[i].size() - 1] + a * dt *
			((list1[i - 1][list1[i - 1].size() - 1] - 2 * list1[i][list1[i].size() - 1] + list1[i + 1][list1[i + 1].size() - 1]) / dx
				+ (list1[i][list1[i].size() - 2] - 2 * list1[i][list1[i].size() - 1] + list2[i][1]) / dx);
	}

	list1 = temp_left;
	temp_left = temp_right;
	temp_right = list3;
	for (int i = 0; i < list1.size(); ++i) {
		list1[i][list1[i].size() - 1] = temp_mid[i];
	}
	temp_mid.resize(list3.size(), 0);

#pragma omp parallel for
	for (int i = 1; i < list3.size() - 1; ++i) {
		for (int j = 1; j < list3[i].size() - 1; ++j) {
			temp_right[i][j] += a * dt * ((list3[i][j + 1] - 2 * list3[i][j] + list3[i][j - 1]) / dx + (list3[i + 1][j] - 2 * list3[i][j] + list3[i - 1][j]) / dx);
		}
	}

	for (int i = 0; i < list3.size(); ++i) {
		temp_right[i][0] = temp_right[i][1];
	}

	for (int j = 0; j < list3[0].size(); ++j) {
		temp_right[0][j] = temp_right[1][j];
		temp_right[temp_right.size() - 1][j] = temp_right[temp_right.size() - 2][j];
	}

	temp_mid[0] = temp_right[0][temp_right[0].size() - 1];
	temp_mid[temp_mid.size() - 1] = temp_right[temp_right.size() - 1][list2[0].size() - 1];

#pragma omp parallel for
	for (int i = 1; i < list3.size() - 1; ++i) {
		temp_mid[i] = list3[i][list3[i].size() - 1] + a * dt *
			((list3[i - 1][list3[i - 1].size() - 1] - 2 * list3[i][list3[i].size() - 1] + list3[i + 1][list3[i + 1].size() - 1]) / dx
				+ (list3[i][list3[i].size() - 2] - 2 * list3[i][list3[i].size() - 1] + list2[list2.size() - list3.size() + i][1]) / dx);
	}
	list2 = temp_left;
	list3 = temp_right;
	temp_right = temp_left;
	temp_left = list4;

	for (int i = 1; i < temp_left.size() - 1; ++i) {
		for (int j = 1; j < temp_left[i].size() - 1; ++j) {
			temp_left[i][j] += a * dt * ((list4[i][j + 1] - 2 * list4[i][j] + list4[i][j - 1]) / dx + (list4[i + 1][j] - 2 * list4[i][j] + list4[i - 1][j]) / dx);
		}
	}

	for (int i = 1; i < list3.size(); ++i) {
		list3[i][list3[i].size() - 1] = temp_mid[i];
		list2[list2.size() - list3.size() + i][0] = temp_mid[i];
	}

	temp_mid.clear();
	temp_mid.resize(list4[0].size(), 0);

	for (int j = 1; j < temp_mid.size() - 1; ++j) {
		temp_mid[j] = temp_right[temp_right.size() - 1][j] +
			a * dt * ((temp_right[temp_right.size() - 2][j] - 2 * temp_right[temp_right.size() - 1][j] + temp_left[1][j]) / dx
				+ (temp_right[temp_right.size() - 1][j - 1] - 2 * temp_right[temp_right.size() - 1][j] + temp_right[temp_right.size() - 1][j + 1]) / dx);
	}

	for (int j = 0; j < temp_mid.size(); ++j) {
		temp_right[temp_right.size() - 1][j] = temp_mid[j];
		temp_left[0][j] = temp_mid[j];
	}

	for (int i = 0; i < temp_left.size(); ++i) {
		temp_left[i][0] = temp_left[i][1];
		temp_left[i][temp_left[i].size() - 1] = temp_left[i][temp_left[i].size() - 2];
	}

	for (int j = 0; j < temp_left[0].size(); ++j) {
		temp_left[temp_left.size() - 1][j] = temp_left[temp_left.size() - 2][j];
	}

	list4 = temp_left;

	for (int i = 1; i < list1.size() - 1; ++i) {
		list2[i][0] = list1[i][list1[i].size() - 1];
	}

	for (int j = 0; j < list2[list2.size() - 1].size(); ++j) {
		list2[list2.size() - 1][j] = list4[0][j];
	}

}

int main(int argc, char** argv) {    // variant, a, dx, dt, Tmax
	setlocale(LC_ALL, "rus");

	if (argc < 6) {
		cout << "Мала параметрив" << endl;

		ofstream f("output1.txt");
		f.close();

		return -1;
	}


	cout << "ABALDET'";
	unsigned int start_time = clock();

	double t, dx = atof(argv[3]), dt = atof(argv[4]), Tmax = atof(argv[5]), a = atof(argv[2]);
	int l1 = 10, l2 = 50, l3 = 10, l4 = 20, l5 = 10, l6 = 20, var = atoi(argv[1]);
	vector<vector<double>> list1(l1), list2(l2), list3(l3), list4(l6);

	vector<double> str_(l4, 1);
	list1[0] = str_;

	for (int i = 1; i < l1; ++i) {
		vector<double> str(l4);
		for (int j = 0; j < l4; ++j) {
			str[j] = 0;
		}
		list1[i] = str;
	}

	str_.resize(l5, 1);
	list2[0] = str_;
	for (int i = 1; i < l2; ++i) {
		vector<double> str(l5 + 1, 0);
		list2[i] = str;
	}

	for (int i = 0; i < l3; ++i) {
		vector<double> str(l4, 0);
		list3[i] = str;
	}

	for (int i = 0; i < l6; ++i) {
		vector<double> str(l5 + 1, 0);
		list4[i] = str;
	}

	t = 0;

	while (t < Tmax) {
		switch (var) {
		case 1:
			solve_1(list1, list2, list3, list4, dx, dt, a);
			break;
		default:
			solve_2(list1, list2, list3, list4, dx, dt, a);
		}

		t += dt;
	}

	unsigned int finish_time = clock();

	unsigned int search_time = finish_time - start_time;

	cout << search_time << endl;

	ofstream f("output1.txt");

	for (int i = 0; i < list1.size(); ++i) {
		for (int j = 0; j < list1[i].size(); ++j) {
			f << j << "\t" << i << "\t" << list1[i][j] << "\n";
		}
	}

	for (int i = 0; i < list2.size(); ++i) {
		for (int j = 0; j < list2[i].size(); ++j) {
			f << j + l4 << "\t" << i << "\t" << list2[i][j] << "\n";
		}
	}

	for (int i = 0; i < list3.size(); ++i) {
		for (int j = 0; j < list3[i].size(); ++j) {
			f << j << "\t" << i + l2 - l3 << "\t" << list3[i][j] << "\n";
		}
	}

	for (int i = 0; i < list4.size(); ++i) {
		for (int j = 0; j < list4[i].size(); ++j) {
			f << j + l4 << "\t" << i + l2 << "\t" << list4[i][j] << "\n";
		}
	}
	f.close();

	return 0;
}