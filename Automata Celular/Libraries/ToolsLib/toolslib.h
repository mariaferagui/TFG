#include <iostream>
#include <fstream> //getMat
#include <cmath> //is_integer
#include <vector> // neighbours
#include <array>
#include "../../constants.hpp"
using namespace constants;

int* get_mat(std::string filename, int matlen);
void node_to_coordinates(int node, int &x, int &y);
void coordinates_to_node(int &node, int x, int y);
std::vector<int> get_neighbours(int *mat, int node);
std::vector<int> get_specific_neighbours(int *mat, int node, int value, char mode);
void create_vec(int node_num, int mat[], int value);
void create_vec(int node_num, float mat[], float value);
void save_mat(int node_num, int mat[], std::string filename);
void save_mat(int node_num, float mat[], std::string filename);
void save_mat(int node_num, double mat[], std::string filename);
void save_vec(std::vector<int> mat, std::string filename);
int* get_random_nodes(int xsize, int ysize);
bool metastasis(int *mat, int xsize, int ysize);
void changeNegativeValue(double &value);
double* int_2_double(int mat[], int matlen);
int cell_counter(int mat[]);
void get_occupied_nodes(int mat[], int mat_nodes[]);
bool is_negative(int value);
bool out_of_mat(int value);
void filename_inc ( std::string *filename );

