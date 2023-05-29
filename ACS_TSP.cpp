// ACS_TSP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>
#include <random>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "read_instance.h"
using namespace std;

double* tourLength;
double bestTourLength;
vector<int> bestTour;
int bestAnt;
int optimal_tour_length;

void generateDist(int n) {
    distMatrix = new int*[n];
    int min = 2; // minimum distance between two cities
    int max = 24; // maximum distance between two cities
    int factor = max - min + 1;
    for (int i = 0; i < n; i++) {
        distMatrix[i] = new int[n];
        for (int j = 0; j < n; j++) {
            if (i == j)
                distMatrix[i][i] = -1;
            else if (i > j)
                distMatrix[i][j] = distMatrix[j][i];
            else if(i < j)
                distMatrix[i][j] = (rand() % factor) + min;
        }
    }
}

int nearest_neighbour(int n) { // implements the nearest neighbour heuristic on the TSP instance
    int current = rand() % n;
    int start_city = current;
    int nn_dist = 0; // the distance covered by the nearest neighbour tour
    vector<int> yet_to_visit;
    for (int i = 0; i < n; i++)
        yet_to_visit.push_back(i);
    yet_to_visit.erase(yet_to_visit.begin() + current);
    while (!yet_to_visit.empty()) {
        int lowest_dist = 100000000;
        int next_index = 0, city_index = 0;
        for (unsigned int iter = 0; iter < yet_to_visit.size(); iter++) {
            city_index = yet_to_visit[iter];
            if (distMatrix[current][city_index] < lowest_dist) {
                next_index = iter;
                lowest_dist = distMatrix[current][city_index];
            }
        }

        nn_dist += distMatrix[current][city_index];
        yet_to_visit.erase(yet_to_visit.begin() + next_index);
        current = city_index;
    }

    nn_dist += distMatrix[current][start_city]; // back to the initial city
    return nn_dist;
}

int complete_enumeration(int n) {
    int tourLength = 100000000;
    vector<int> perm;
    for (int i = 0; i < n; i++)
        perm.push_back(i);
    do {
        int temp = 0;
        for (int i = 1; i < n; i++) {
            int city1 = perm[i - 1];
            int city2 = perm[i];
            temp += distMatrix[city1][city2];
        }
        temp += distMatrix[perm[n-1]][perm[0]];
        if (tourLength > temp)
            tourLength = temp;
    } while (next_permutation(perm.begin(), perm.end()));

    return tourLength;
}

int roulette_wheel(vector<double> arr, double sum_arr) {
    int i = 0;
    double q = (float)rand() / RAND_MAX;
    double threshold = q * sum_arr;
    double sum_slice = arr[i];
    while (sum_slice < threshold) {
        i++;
        sum_slice += arr[i];
    }

    return i;
}

int main()
{
    /**
    * A C++ implementation of Ant Colony System (ACS) -- basic version for TSPs.
    * @Since 07:45; 03-05-2023
    * @Ref Dorigo, M., & Gambardella, L. M. (1997). Ant colony system: a cooperative learning approach to the traveling 
       salesman problem. IEEE Transactions on evolutionary computation, 1(1), 53-66.
    **/
    int n; // number of cities in the TSP
    double alpha = 0.10, rho = 0.10, beta = 2.5, q_0 = 0.9; // parameters
    cout << "***************************************************\n";
    cout << " Ant Colony System (ACS) implementation on TSPs!\n";
    cout << "***************************************************\n";

    bool random = false;

    if (random) { // instances to be generated at random
        n = 400;
        generateDist(n); // generate the distance between the n cities randomly
    }
    else { // read from a dataset, e.g. from a TSPLIB dataset
        loadInstance();
        n = ins_size;
    }

    int trials = 0;
    int n_trials = 10;
    int* fitness = new int[n_trials];
    while (trials < n_trials) {
        srand((unsigned)time(NULL));
        if (n <= 10) {
            optimal_tour_length = complete_enumeration(n);
            cout << "Optimal Tour Length = " << optimal_tour_length << " units";
        }

        int nntour_length = nearest_neighbour(n);
        cout << "\nNearest Neighbour Heuristic Tour Length = " << nntour_length << " units";

        //initialization
        int iter = 0; // number of iterations
        int m = 20; // number of ants // 25
        int* init_loc_ants = new int[m]; // an array storing intial location of the m ants
        bestTourLength = 100000000.0;
        tourLength = new double[m];
        double init_pher = 1.0 / ((double)n * (double)nntour_length);
        double** pher;
        double** heu_info; // to be computed once
        pher = new double* [n];
        heu_info = new double* [n];

        for (int i = 0; i < n; i++) {
            pher[i] = new double[n];
            heu_info[i] = new double[n];
        }

        for (int i = 0; i < n; i++) {
            pher[i][i] = 0.0;
            heu_info[i][i] = 0.0;
            for (int j = i + 1; j < n; j++) {
                pher[i][j] = init_pher;
                pher[j][i] = pher[i][j];
                heu_info[i][j] = 1.0 / ((double)distMatrix[i][j]);
                heu_info[j][i] = heu_info[i][j];
            }
        }

        vector<int>* memory; // a vector, indicating the cities yet to be visited by the current ant k'
        vector<int>* tour; // a vector, stores the tour of an ant k
        tour = new vector<int>[m];
        memory = new vector<int>[m];

        while (iter < 1250) {
            int i;
            bool bestTourFound = false;

            for (int k = 0; k < m; k++) // place k ants into initial positions
                init_loc_ants[k] = (int)rand() % n;

            for (int k = 0; k < m; k++) {
                memory[k].clear();
                for (int c = 0; c < n; c++)
                    memory[k].push_back(c);
                tourLength[k] = 0.0;
            }

            for (i = 0; i < n; i++) {
                int cur_city_k;
                int choice_index, selected_city;
                for (int k = 0; k < m; k++) { //iteration of ants' activities
                    if (i == 0) {
                        memory[k].erase(memory[k].begin() + init_loc_ants[k]);
                        cur_city_k = init_loc_ants[k]; // current location of ant k
                        tour[k].clear();
                        tour[k].push_back(cur_city_k); // update the tour of ant i
                    }
                    else {
                        cur_city_k = tour[k][i];
                    }

                    if (i < n - 1) {
                        vector<double> utility_values_k; // utility values of the cities yet to be visited by k when in city r
                        int max_index = 0;
                        double sum_of_utilities = 0.0;
                        double max_utility = -100000000.0;
                        for (size_t s = 0; s < memory[k].size(); s++) {
                            int cand_city_k = memory[k][s];
                            double h_info = heu_info[cur_city_k][cand_city_k];
                            double utility = pher[cur_city_k][cand_city_k] * pow(h_info, beta);
                            if (utility > max_utility) {
                                max_utility = utility;
                                max_index = s;
                            }
                            sum_of_utilities += utility;
                            utility_values_k.push_back(utility);
                        }

                        double q = (float)rand() / RAND_MAX; // generate a random real number between 0 and 1
                        if (q < q_0) {
                            //perform exploitation
                            choice_index = max_index;
                            selected_city = memory[k][choice_index];
                        }
                        else {
                            // perform exploration
                            choice_index = roulette_wheel(utility_values_k, sum_of_utilities);
                            selected_city = memory[k][choice_index];
                        }

                        tour[k].push_back(selected_city); // update the tour of ant i
                        tourLength[k] += distMatrix[cur_city_k][selected_city];
                        memory[k].erase(memory[k].begin() + choice_index);
                    }
                    else {
                        tour[k].push_back(init_loc_ants[k]); // go back to the starting city
                        tourLength[k] += distMatrix[cur_city_k][init_loc_ants[k]]; // update the distance travelled.
                        if (bestTourLength > tourLength[k]) {
                            if (!bestTourFound)
                                bestTourFound = true;
                            bestTourLength = tourLength[k];
                            bestTour.clear();
                            // update the best tour
                            for (size_t t = 0; t < tour[k].size(); t++)
                                bestTour.push_back(tour[k][t]);
                            bestAnt = k;
                        }
                    }
                } // loops runs from 1 to m; where k E {1, 2, ..., m} is the index of an ant

                // local pheromone update on the edges ants just travelled through.
                for (int l = 0; l < m; l++) {
                    int city1 = tour[l][i];
                    int city2 = tour[l][i + 1];
                    pher[city1][city2] = (1 - rho) * pher[city1][city2] + rho * init_pher;
                    pher[city2][city1] = pher[city1][city2]; // comment out if TSP is Asymmetric
                }
            } // end of major for loop, that runs from 1 to n.

            // time for global pheromone update
            double inv_bestTour = 1.0 / (double)bestTourLength;
            for (size_t x = 1; x < bestTour.size(); x++) {
                int city1 = bestTour[x - 1];
                int city2 = bestTour[x];
                pher[city1][city2] = (1 - alpha) * pher[city1][city2] + alpha * inv_bestTour;
                pher[city2][city1] = pher[city1][city2]; // comment out if TSP is Asymmetric
            }

            if (bestTourFound) {
                cout << "\nNew best tour with a length of " << bestTourLength << " units found in iteration " << iter + 1;
                //cout << " Init Location of Ant 0 = " << init_loc_ants[0];
            }
            iter++;
        } // end of while (iter < ?) {...} to control the number of epochs

        cout << "\nBest tour after ACS: " << bestTourLength;
        cout << "\nProgram ends...\n";
        fitness[trials] = bestTourLength;
        trials++;
    }

    cout << "\n\nFitness values after trials include\n";
    double sum_fitness = 0.0;
    for (int i = 0; i < n_trials; i++) {
        cout << fitness[i] << endl;
        sum_fitness += fitness[i];
    }

    double av_fitness = (double)sum_fitness / n_trials;
    cout << "\nMean Fitness value after trials is: " << av_fitness << endl;
}
