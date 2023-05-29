#pragma once
#include <iostream>
#include <sstream>
#include <type_traits>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

int** distMatrix;
int ins_size;

int euc_2d(double x_i, double x_j, double y_i, double y_j) {
    double x_d = x_i - x_j;
    double y_d = y_i - y_j;
    double pre_dist = sqrt(x_d * x_d + y_d * y_d);
    int dist = (int)(pre_dist + 0.5);
    return dist;
}

int ceil_2d(double x_i, double x_j, double y_i, double y_j) {
    double x_d = x_i - x_j;
    double y_d = y_i - y_j;
    double pre_dist = sqrt(x_d * x_d + y_d * y_d);
    int dist = (int)ceil(pre_dist);
    return dist;
}

int euc_3d(double x_i, double x_j, double y_i, double y_j, double z_i, double z_j) {
    double x_d = x_i - x_j;
    double y_d = y_i - y_j;
    double z_d = z_i - z_j;
    double pre_dist = sqrt(x_d * x_d + y_d * y_d + z_d * z_d);
    int dist = (int)(pre_dist + 0.5);
    return dist;
}

int man_2d(double x_i, double x_j, double y_i, double y_j) {
    double x_d = abs(x_i - x_j);
    double y_d = abs(y_i - y_j);
    double pre_dist = x_d + y_d;
    int dist = (int)(pre_dist + 0.5);
    return dist;
}

int man_3d(double x_i, double x_j, double y_i, double y_j, double z_i, double z_j) {
    double x_d = abs(x_i - x_j);
    double y_d = abs(y_i - y_j);
    double z_d = abs(z_i - z_j);
    double pre_dist = x_d + y_d + z_d;
    int dist = (int)(pre_dist + 0.5);
    return dist;
}

int max_2d(double x_i, double x_j, double y_i, double y_j) {
    double x_d = abs(x_i - x_j);
    double y_d = abs(y_i - y_j);
    int dist = max((int)(x_d + 0.5), (int)(y_d + 0.5));
    return dist;
}

int max_3d(double x_i, double x_j, double y_i, double y_j, double z_i, double z_j) {
    double x_d = abs(x_i - x_j);
    double y_d = abs(y_i - y_j);
    double z_d = abs(z_i - z_j);
    int pre_dist = max((int)(x_d + 0.5), (int)(y_d + 0.5));
    int dist = max(pre_dist, (int)(z_d + 0.5));
    return dist;
}

int att(double x_i, double x_j, double y_i, double y_j) { // psuedo-euclidean
    double x_d = x_i - x_j;
    double y_d = y_i - y_j;
    double r = sqrt((x_d * x_d + y_d * y_d) / 10.0);
    int t = (int)(r + 0.5);
    int dist = (t < r) ? t + 1 : t;
    return dist;
}

void process_explicit(vector<string> v_data, string edge_weight_type, string edge_weight_format, string dim) {
    bool has_diag = edge_weight_format.find("DIAG") != string::npos ? true : false;
    int n = stoi(dim);
    ins_size = n;
    int k = 0;
    distMatrix = new int* [n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            distMatrix[i] = new int[n];

    if (edge_weight_format.find("UPPER") != string::npos) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    distMatrix[i][j] = has_diag ? stoi(v_data[k++]) : 0;
                    if (distMatrix[i][j] > 0) {
                        cout << "Diagonal cell has an entry other than 0, that's an error...";
                        exit(-1);
                    }
                }
                else if (i > j)
                    distMatrix[i][j] = 0;
                else
                    distMatrix[i][j] = stoi(v_data[k++]);
            }
        }
    }
    else if (edge_weight_format.find("LOWER") != string::npos) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    distMatrix[i][j] = has_diag ? stoi(v_data[k++]) : 0;
                    if (distMatrix[i][j] > 0) {
                        cout << "Diagonal cell has an entry other than 0, that's an error...";
                        exit(-1);
                    }
                }
                else if (i < j)
                    distMatrix[i][j] = 0;
                else
                    distMatrix[i][j] = stoi(v_data[k++]);
            }
        }
    }
    else if (edge_weight_format.find("FULL") != string::npos) { // full matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distMatrix[i][j] = stoi(v_data[k++]);
            }
        }
    }
    else {
        cout << "No support for the edge weight format yet";
        exit(-1);
    }
}

/*
v_data := coordinates of the n cities of the instance
dim := number of cities
code := 0 (for euc_2d), 1 (for ceil_2d), 2 (for att)
*/
void process_euc_2d(vector<string> v_data, string dim, int code) {
    int n = stoi(dim);
    ins_size = n;
    int k = 0;
    distMatrix = new int* [n];
    double** coord = new double* [n];
    for (int i = 0; i < n; i++) {
        coord[i] = new double[2];
        for (int j = 0; j < 2; j++) {
            int l = j == 0 ? k + 1 : k + 2;
            coord[i][j] = stod(v_data[l]);
        }
        k += 3;
    }

    for (int i = 0; i < n; i++)
        distMatrix[i] = new int[n];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                distMatrix[i][j] = 0;
            else {
                if (code == 0)
                    distMatrix[i][j] = euc_2d(coord[i][0], coord[j][0], coord[i][1], coord[j][1]);
                else if (code == 1)
                    distMatrix[i][j] = ceil_2d(coord[i][0], coord[j][0], coord[i][1], coord[j][1]);
                else
                    distMatrix[i][j] = att(coord[i][0], coord[j][0], coord[i][1], coord[j][1]);
                distMatrix[j][i] = distMatrix[i][j];
            }
        }
    }
}

void loadInstance() {
    ifstream fin;
    string instances[] = {"berlin52", "kroA100", "att48", "rat575", "a280", "eil101", "bayg29", "eil51", "eil76"};
    string ins_name = instances[0]; // 6
    string ins_path = ins_name + ".tsp";
    cout << "TSP data: " << ins_name << endl;
    //cin >> ins_name;
    fin.open(ins_path);
    if (!fin.fail()) {
        string line = "";
        int line_number = 0;
        vector<string> raw_matrix_data;
        bool data_section = false;
        string edge_weight_type;
        string dimension;
        string edge_weight_format;
        while (!fin.eof()) {
            getline(fin, line);
            line_number++;
            if (line_number < 4)
                continue;
            if (data_section) {
                stringstream stream(line);
                //stream << line;
                string next_token;
                while (stream >> next_token) {
                    raw_matrix_data.push_back(next_token);
                }
            }
            else if (dimension == "") {
                int _start = line.find(": ") + strlen(": ");
                string mysubtr = line.substr(_start, line.size() - _start); // get the type of edge weight.
                dimension = mysubtr; // substring after removing the trailing spaces;
            }
            else if (line_number == 6 && edge_weight_type == "EXPLICIT") {
                int _start = line.find(": ") + strlen(": ");
                string mysubtr = line.substr(_start, line.size() - _start); // get the type of edge weight.
                edge_weight_format = mysubtr; // substring after removing the trailing spaces;
            }
            else if (edge_weight_type == "") {
                size_t found = line.find("EXPLICIT") || line.find("EUC_2D") || line.find("ATT") || line.find("EUC_3D") ||
                    line.find("MAX_2D") || line.find("MAX_3D") || line.find("MAN_2D") || line.find("MAN_3D") ||
                    line.find("CEIL_2D") || line.find("GEO");
                if (found != string::npos) {
                    int _start = line.find(": ") + strlen(": ");
                    string mysubtr = line.substr(_start, line.size() - _start); // get the type of edge weight.
                    edge_weight_type = mysubtr; // substring after removing the trailing spaces;
                }
            }
            else if (line.find("SECTION") != string::npos)
                data_section = true;
        }

        if (edge_weight_type == "EXPLICIT")
            process_explicit(raw_matrix_data, edge_weight_type, edge_weight_format, dimension);
        else if (edge_weight_type == "EUC_2D")
            process_euc_2d(raw_matrix_data, dimension, 0);
        else if (edge_weight_type == "ATT")
            process_euc_2d(raw_matrix_data, dimension, 2);
        else {
            cout << "No support for the data format of this instance yet..";
            exit(-1);
        }
    }
    else {
        cout << "Issues with opening instance data file, confirm the name and location of the file and try again";
        exit(-1);
    }
}
