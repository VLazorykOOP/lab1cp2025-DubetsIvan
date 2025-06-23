#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>
#include <stdexcept>  
using namespace std;

struct TableRow {
    double x, T, U;
};

bool readTable(const string& filename, vector<TableRow>& table) {
    ifstream file(filename);
    if (!file.is_open()) return false;

    TableRow row;
    while (file >> row.x >> row.T >> row.U)  
        table.push_back(row);
    return true;
}

pair<double, double> interpolate(const vector<TableRow>& table, double x) {
    for (size_t i = 0; i + 1 < table.size(); ++i) {
        if (table[i].x <= x && x <= table[i + 1].x) {
            double x0 = table[i].x, x1 = table[i + 1].x;
            double T0 = table[i].T, T1 = table[i + 1].T;
            double U0 = table[i].U, U1 = table[i + 1].U;
            double T = T0 + (T1 - T0) * (x - x0) / (x1 - x0);
            double U = U0 + (U1 - U0) * (x - x0) / (x1 - x0);
            return { T, U };
        }
    }
    return { 0, 0 }; 
}

pair<double, double> get_T_U(double x, bool& ok) {
    try {
        vector<TableRow> table;
        string filename;

        if (x <= 1)
            filename = "dat_X_1_1.dat";
        else if (x < -1) {
            x = 1.0 / x;
            filename = "dat_X00_1.dat";
        }
        else {
            x = 1.0 / x;
            filename = "dat_X1_00.dat";
        }

        if (!readTable(filename, table)) {
            ok = false;
            return { 0, 0 };
        }

        ok = true;
        return interpolate(table, x);
    }
    catch (const exception& e) {
        cerr << "Exception in get_T_U: " << e.what() << endl;
        ok = false;
        return { 0, 0 };
    }
}

double T(double x, bool& ok) {
    return get_T_U(x, ok).first;
}

double U(double x, bool& ok) {
    return get_T_U(x, ok).second;
}

double Srz(double x, double y, double z);

double Srs(double x, double y, double z) {
    if (z + x * y <= 0) return Srz(x, y, z);
    if (x + z * y <= 0) return Srz(x, y, z);
    return 0;
}

double Qrz(double x, double y, double a, double b, double (*Srs_func)(double, double, double)) {
    if (y >= x) return Srs_func(y, x, y);
    return Srs_func(x, y, x);
}

double Srz(double x, double y, double z) {
    bool ok1, ok2, ok3;
    double Ux = U(x, ok1), Uy = U(y, ok2), Uz = U(z, ok3);
    double Tx = T(x, ok1), Ty = T(y, ok2);

    if (!(ok1 && ok2 && ok3))
        return 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;

    if (Uy + Uz - x <= 0)
        return Ty + Uy + Uz;
    else
        return Tx + Uz + Ty;
}

double Qrz1(double x, double y, double a, double b) {
    return Qrz(x, y, a, b, Srs);
}

double Rrz_algo2(double x, double y, double z) {
    if (x <= z)
        return Qrz1(x, y, x, y);
    return Qrz1(y, z, x, y);
}

double Rrz_algo3(double x, double y, double z) {
    if (y <= z)
        return Qrz(x, y, x, y, Srs);
    return Qrz(y, z, x, y, Srs);
}

double Rrz(double x, double y, double z) {
    if (z + x * y <= 0)
        return Rrz_algo2(x, y, z);
    else if (x + z * y <= 0)
        return Rrz_algo3(x, y, z);
    else if (y <= x)
        return Qrz(x, y, x, y, Srs);
    else
        return Qrz(y, z, x, y, Srs);
}


double Grs(double x, double y, double z) {
    return 0.1389 * Rrz(x, y, y) + 1.8389 * Rrz(x - y, z, y);
}

double fun(double x, double y, double z) {
    bool ok;
    get_T_U(x, ok);

    if (!ok)
        return 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;

    return x * Grs(x, y, z) + y * Grs(x, z, y);
}

int main() {
    try {
        double x, y, z;
        cout << "Enter x, y, z: ";
        if (!(cin >> x >> y >> z))
            throw runtime_error("Invalid input");

        double result = fun(x, y, z);
        cout << "fun(x, y, z) = " << result << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
