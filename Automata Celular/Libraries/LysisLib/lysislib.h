#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <time.h>
#include <cmath>
#include <ctime>
#include <fstream> 
#include <cmath> 
#include <vector>
#include "../../constants.hpp"
#include "../ToolsLib/toolslib.h"

using namespace constants;

// Alphabetiacally ordered
void effectorCellPlacement(int T[], int E[]);
void Emigration(int T[], int E[], int Ecount[], int H[], int D[], int node, std::mt19937 generator);
void first_quad(int E[], int T[]);
void fourth_quad(int E[], int T[]);
float get_lysis_ratio(int T[], int T0);
void inactivation(int T[], int E[], int Ecount[], int H[], int node, std::mt19937 generator);
void lysis(int T[], int E[], int Ecount[], int D[], int H[], int node, std::mt19937 generator);
bool no_cells(int mat[]);
bool no_lysis(int T[], int T0);
void recruitment(int T[], int E[], int D[], int H[], int node, std::mt19937 generator);
void second_quad(int E[], int T[]);
void sector(int E[], int T[], int quadrant);
int summation(int mat[], std::vector<int> neighbours);
void third_quad(int E[], int T[]);
void tumor_lysis(int T[], int E[], int Ecount[], int D[], int H[]);









