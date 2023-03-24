# include <iostream>
# include "Libraries/LysisLib/lysislib.h"
# include "Libraries/ToolsLib/toolslib.h"
# include "constants.hpp"
using namespace constants;

float lysis_ratio;
int *T, *Ecount, *D, *H, *E, *T_cells; 
int xsize = 2*NX-1, ysize = 2*NY-1, Tdi, Tdf, diff = 0, Tdead = 0;
int N_T_cells, d = 0, no_lysis_counter = 0, no_lysis_lim = 3;
std::string folder = "Results/Lysis/";
std::string T_filename = folder + "T/0000.txt";
std::string E_filename = folder + "E/0000.txt";

int main(){
    T = new int[NODE_NUM];
    D = new int[NODE_NUM];
    E = new int[NODE_NUM];
    H = new int[NODE_NUM];
    Ecount = new int[NODE_NUM];
    create_vec(NODE_NUM, E, 0);
    create_vec(NODE_NUM, D, 0);
    create_vec(NODE_NUM, H, 1);
    create_vec(NODE_NUM, Ecount, 0);
    std::vector<int> Tcells_num;
    T = get_mat("Results/Growth/T/0275.txt", NODE_NUM); 
    
    // Create H mat
    for(int node = 0;node<NODE_NUM;node++){
        if(T[node]==1) H[node]=0;
    }

    // Place effector cells and save matrix
    effectorCellPlacement(T, E);
    save_mat(NODE_NUM, E, E_filename);
    save_mat(NODE_NUM, T, T_filename);

    int T0 = cell_counter(T);

    for(int i=1; i<DESTRUCTION_IT; i++){
        
        // Lysis
        tumor_lysis(T, E, Ecount, D, H);

        // Increase filename
        filename_inc ( &T_filename );
        filename_inc ( &E_filename );

        // Save each 5 time units
        if(i%5 == 0){
            save_mat(NODE_NUM, T,T_filename);
            save_mat(NODE_NUM, E,E_filename);
        }
        
        // Stop if no presence of E or T cells
        if((no_cells(E)) || (no_cells(T))){
            break;
        };

        // Stop if 0 lysis in 3 iterations
        if(no_lysis(T, T0)){
            no_lysis_counter++;
            if (no_lysis_counter>=no_lysis_lim){
                break;
            }
        }
        T0 = cell_counter(T);
        Tcells_num.push_back(T0);
        std::cout<<"ITERATION: "<<i<<std::endl;
    }

    // Print lysis ratio
    lysis_ratio = get_lysis_ratio(T, T0);
    std::cout<<"Lisis ratio: "<<lysis_ratio<<std::endl;
    save_vec(Tcells_num, "Results/Lysis/T_cell_count.txt");

    return 0;
}