#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <time.h>
#include <vector>
#include "../../constants.hpp"
#include "../ToolsLib/toolslib.h"

using namespace constants; //quitar???

void necrosis(double M[], int T[], int D[], int node, std::mt19937 generator);
void division(double N[], int T[],  int D[], int H[], float DIV_mat[], int node, int xsize, int ysize, std::mt19937 generator);
void migracion(double M[], int T[], int D[], int H[], float DIV_mat[], int node, int xsize, int ysize, std::mt19937 generator);
void grow(double M[], double N[], int T[], int D[], int H[], float DIV_mat[], int xsize, int ysize);
void changeNegativeValue(double &value);
void mutation(int T[], float DIV_mat[]);
