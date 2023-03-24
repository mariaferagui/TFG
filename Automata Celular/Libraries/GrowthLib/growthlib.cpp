#include "growthlib.h"


void grow(double M[], double N[], int T[], int D[], int H[], float DIV_mat[], int xsize, int ysize){  
    /* Main function. Call all the posible actions that a tumor cell can take. */
    
    int var, dice, *random_nodes, node;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> dice_distribution(1,3);
    
    if(MUTATED_CELLS>0)mutation(T, DIV_mat);
    random_nodes = new int[xsize*ysize];
    random_nodes = get_random_nodes(xsize, ysize);
    for(int i = 0; i < xsize*ysize; i++)
    {   
        node = random_nodes[i];
        if(T[node] >= 1){
            dice = dice_distribution(generator);
            switch (dice)
            {
            case 1:
                necrosis(M, T, D, node, generator);
                break;
            case 2: 
                migracion(M, T, D, H, DIV_mat, node, xsize, ysize, generator);
                break;
            case 3:
                division(N, T, D, H, DIV_mat, node, xsize, ysize, generator);
                break;
            } 
        }
    } 
}

void necrosis(double M[], int T[], int D[], int node, std::mt19937 generator){
    /* Decides whether a tumor cell dies or not. If so, computes P probability 
    and if its greated than a random number modifies D matrix (Dead cells).*/

    std::normal_distribution<double> distribution(0,1);
    double rnd_n = distribution(generator), P; 
    P = exp(- pow(M[node]/(T[node]*NEC) ,2));
    
    if(P>fabs(rnd_n)){ 
        T[node]--;
        if(T[node] == 0){
            D[node] = 1;
        }
    }
    return;
}

void migracion(double M[], int T[], int D[], int H[], float DIV_mat[], int node, int xsize, int ysize, std::mt19937 generator){
    /* Determines if a tumor cell migrates. Exchanges two pixels and modifies T, H and D matrices.*/

    std::normal_distribution<double> distribution(0,1);
    std::uniform_int_distribution<int> dice_distr(1,50);
    std::vector<int> neighbour_nodes, free_nodes;
    int neighbour_node, index;
    double rnd_n = distribution(generator);
    double P = 1 - exp( - (T[node] * pow(M[node]/MIG,2)));
    
    if(P>fabs(rnd_n)){
        neighbour_nodes = get_neighbours(T, node);
        free_nodes = get_specific_neighbours(T, node, 0, '=');
        if( !free_nodes.empty() ){ // two neighbour nodes
            index = dice_distr(generator) % free_nodes.size();
            neighbour_node = free_nodes[index];
            if(T[node] > 1){ // unmodified if there is more than one cell
                T[node]--;
                T[neighbour_node] = 1;
                DIV_mat[neighbour_node] = DIV_mat[node];
                if (H[neighbour_node] == 1){
                    H[neighbour_node] = 0;
                }
                else if (D[neighbour_node] == 1) {
                    D[neighbour_node] = 0;
                }
            }
            else{
                T[node]--;
                T[neighbour_node]++;
                DIV_mat[neighbour_node] = DIV_mat[node];
                DIV_mat[node] = DIV;
                if (H[neighbour_node] == 1){
                    H[neighbour_node] = 0;
                    H[node] = 1;
                }
                else if (D[neighbour_node] == 1) {
                    D[neighbour_node] = 0;
                    D[node] = 1;
                }
            }
        }
        else{
            index = dice_distr(generator) % neighbour_nodes.size();
            neighbour_node = neighbour_nodes[index];
            T[node]--;
            T[neighbour_node]++;
            DIV_mat[neighbour_node] = DIV_mat[node];
            if(T[node]==0){
                H[node] = 1;
                DIV_mat[node] = DIV;
            }
        }
    }
};

void mutation(int T[], float DIV_mat[]){
    /* Asigns mutations to tumor cells.*/ 

    int var, *T_nodes, number_T_cells, node, index;
    static int mutated = 0;
    T_nodes = new int[NODE_NUM];
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
    static_cast<long long>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now().time_since_epoch()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var)))};
    std::knuth_b generator(seed);
    std::uniform_int_distribution<int> dice_distribution(1,4);
    int rnd = dice_distribution(generator);
    std::normal_distribution<float> distribution(2*DIV/3,0.1);

    number_T_cells = cell_counter(T);
    if( (number_T_cells>2000) & (mutated<MUTATED_CELLS) & (rnd == 2) ){
        get_occupied_nodes(T, T_nodes);
        index = dice_distribution(generator) % number_T_cells;
        node = T_nodes[index];
        DIV_mat[node] = fabs(distribution(generator));
        mutated++;
    }
};

void division(double N[], int T[], int D[], int H[], float DIV_mat[], int node, int xsize, int ysize, std::mt19937 generator){
    /* Determines the division of the cell. If affirmative, it carries it out modifying T, DIV_mat, H and D.*/

    std::normal_distribution<double> distribution(0,1);
    std::uniform_int_distribution<int> dice_dist(1,50);
    std::vector<int> free_nodes;
    int neighbour_node, index;
    double rnd_n = distribution(generator), P;
    
    P = 1 - exp( - pow(N[node]/(T[node]*DIV_mat[node]), 2));
    if(P>fabs(rnd_n)){
        free_nodes = get_specific_neighbours(T,node,0,'=');
        if( free_nodes.empty() ){
            T[node]++;
        }
        else{
            index = dice_dist(generator) % free_nodes.size();
            neighbour_node = free_nodes[index];
            T[neighbour_node] = 1;
            DIV_mat[neighbour_node] = DIV_mat[node];
            if(H[neighbour_node] == 1){
                H[neighbour_node] = 0;
            }
            else{
                D[neighbour_node] = 0;
            }
        }
    }
    return;
};
