#include "lysislib.h"

// Functionally ordered

void tumor_lysis(int T[], int E[], int Ecount[], int D[], int H[]){
    int var;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> dice_distribution(1,3);
    std::vector<int> Tneighbours;
    int dice, node = 0;
    for(int node = 0; node < NODE_NUM; node++){
        if(E[node]!=0){
            Tneighbours = get_specific_neighbours(T, node, 0,'>');
            if(Tneighbours.size()!=0){
                lysis(T, E, Ecount, D, H, node, generator);
            }
            else{
                dice = dice_distribution(generator);
                switch (dice)
                {
                case 1:
                    inactivation(T, E, Ecount, H, node, generator);
                    break;
                case 2:
                    recruitment(T, E, D, H, node, generator);
                    break;
                case 3:
                    Emigration(T, E, Ecount, H, D, node, generator);
                    break;
                }
            } 
        }
    }
};

void lysis(int T[], int E[], int Ecount[], int D[], int H[], int node, std::mt19937 generator){
    std::vector<int> Eneighbours, Tneighbours, Hneighbours, Dneighbours;
    std::normal_distribution<double> distribution(0,1);
    std::uniform_int_distribution<> u_distrib(1,50);
    int neignode, index;
    double rnd_n = distribution(generator), P;
    Eneighbours = get_specific_neighbours(E, node, 0, '>');
    P = 1-exp(-pow(Eneighbours.size()/LYS,2));
    if(P>fabs(rnd_n)){
        Tneighbours = get_specific_neighbours(T, node, 0, '>');
        index = u_distrib(generator) % Tneighbours.size();
        neignode = Tneighbours[index];
        T[neignode]--; D[neignode]++; Ecount[node]++;

        recruitment(T, E, D, H, node, generator);

        if(Ecount[node] == 3){
            E[node] = 0; H[node] = 1; Ecount[node] = 0;
        }
    }
};

void recruitment(int T[], int E[], int D[], int H[], int node, std::mt19937 generator){
    std::normal_distribution<double> distribution(0,1);
    std::uniform_int_distribution<> u_distrib(1,50);
    std::vector<int> Tneighbours, Hneighbours, Dneighbours;
    int index, neignode;
    double rnd_n = distribution(generator), P;
    

    Tneighbours = get_specific_neighbours(T, node, 0, '>');

    P = exp(-1/pow((summation(T, Tneighbours)*REC),2));


    Hneighbours = get_specific_neighbours(H, node, 0, '>');
    Dneighbours = get_specific_neighbours(D, node, 0, '>');
    std::vector<int> HDneighbours = Hneighbours; // set union of both vectors
    HDneighbours.insert(HDneighbours.end(), Dneighbours.begin(), Dneighbours.end());
    for(int j = 0; j<HDneighbours.size(); j++){
        rnd_n = distribution(generator);
        if(P>fabs(rnd_n)){
            index = u_distrib(generator) % HDneighbours.size();
            neignode = HDneighbours[index];
            if(T[neignode] == 0){
                D[neignode] = 0; H[neignode] = 0; E[neignode] = 30;
                HDneighbours.erase(HDneighbours.begin()+index);
            }   
        }
    }
};

void inactivation(int T[], int E[], int Ecount[], int H[], int node, std::mt19937 generator){
    std::normal_distribution<double> distribution(0,1); 
    std::vector<int> Tneighbours;
    double rnd_n = distribution(generator), P;
    Tneighbours = get_specific_neighbours(T, node, 0, '>');
    P = 1 - exp(- pow(1/(summation(T, Tneighbours)*INC),2));
    if(P>fabs(rnd_n)){
        E[node] = 0; Ecount[node] = 0;
        H[node] = 1;
    }
};

int summation(int mat[], std::vector<int> neighbours){
    int node, result = 0;
    for(int i = 0; i<neighbours.size(); i++){
        node = neighbours[i];
        result += mat[node];
    }
    return result;
};

void Emigration(int T[], int E[], int Ecount[], int H[], int D[], int node, std::mt19937 generator){
    std::uniform_int_distribution<int> u_distrib(1,8);
    std::vector<int> neighbours;
    int index, neignode, temp;

    neighbours = get_neighbours(E, node);
    index = u_distrib(generator) % neighbours.size();
    neignode = neighbours[index];

    if(E[neignode] != 0){
        temp = Ecount[neignode];
        Ecount[neignode] = Ecount[node];
        Ecount[node] = temp;

        temp = E[neignode];
        E[neignode] = E[node]; 
        E[node] = temp;
    }

    else if( (H[neignode] == 1 ) ){
        E[neignode] = E[node]; E[node] = 0;
        Ecount[neignode] = Ecount[node]; Ecount[node] = 0; 
        H[node] = H[neignode]; H[neignode] = 0;
    }

    else if(((D[neignode] > 0) & (T[neignode] == 0))){
        E[neignode] = E[node];
        Ecount[neignode] = Ecount[node]; Ecount[node] = 0;
        D[neignode] = 0;
    }
};

bool no_cells(int mat[]){
    bool result = true;
    for(int node=0; node<NODE_NUM; node++){
        if(mat[node]!=0){
            result = false;
            break;
        }
    }
    return result;
}

void effectorCellPlacement(int T[], int E[]){
    std::vector<int> Tneighbours, Eneighbours;
    int ECells;
    sector(E, T, QUADRANT);
}

void sector(int E[], int T[], int quadrant){

    switch (quadrant)
    {
    case 1:
        first_quad(E,T);
        break;
    case 2:
        second_quad(E,T);
        break;
    case 3:
        third_quad(E,T);
        break;
    case 4:
        fourth_quad(E,T);
        break;
    }

};

void first_quad(int E[], int T[]){
    int node, n_cells, var, x, y;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    n_cells = int(E_PERCENTAGE*cell_counter(T));
    std::uniform_int_distribution<int> x_distribution(0, X_SIZE/2);
    std::uniform_int_distribution<int> y_distribution(Y_SIZE/2, Y_SIZE);
    while(n_cells>0){
        x = x_distribution(generator);
        y = y_distribution(generator);
        coordinates_to_node(node, x,y);
        if(T[node] == 0){
            E[node] = 1;
            n_cells--;
        }
    }
};

void second_quad(int E[], int T[]){
    int node, n_cells, var, x, y;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> x_distribution(0, X_SIZE);
    std::uniform_int_distribution<int> y_distribution(Y_SIZE/2, Y_SIZE); 
    n_cells = int(E_PERCENTAGE*cell_counter(T));
    while(n_cells>0){
        x = x_distribution(generator);
        y = y_distribution(generator);
        coordinates_to_node(node, x,y);
        if(T[node] == 0){
            E[node] = 1;
            n_cells--;
        }
    }
};

void third_quad(int E[], int T[]){
    int node, n_cells, var, x, y;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> x_distribution(0, X_SIZE);
    std::uniform_int_distribution<int> y_distribution(0, Y_SIZE); //0, Y_SIZE/2
    n_cells = int(E_PERCENTAGE*cell_counter(T));
    while(n_cells>0){
        x = x_distribution(generator);
        y = y_distribution(generator);
        coordinates_to_node(node, x,y);
        if((T[node] == 0)){
            if((x<(X_SIZE/2)) && (y<(Y_SIZE/2))){
                E[node] = 0;
            }
            else{
                E[node] = 1;
                n_cells--;
            }

        }
    }
};

void fourth_quad(int E[], int T[]){
    int node, n_cells, var, x, y;
    std::seed_seq seed{static_cast<long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
                        static_cast<long long>(reinterpret_cast<intptr_t>(&var))};
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> x_distribution(0, X_SIZE);
    std::uniform_int_distribution<int> y_distribution(0, Y_SIZE);
    n_cells = int(E_PERCENTAGE*cell_counter(T));
    while(n_cells>0){
        x = x_distribution(generator);
        y = y_distribution(generator);
        coordinates_to_node(node, x,y);
        if(T[node] == 0){
            E[node] = 1;
            n_cells--;
        }
    }
};

float get_lysis_ratio(int T[], int T0){
    int Tactual = cell_counter(T);
    float lysis;
    lysis = (float(T0)-float(Tactual))/float(T0);
    return lysis;
};

bool no_lysis(int T[], int T0){
    int Tactual = cell_counter(T);
    bool no_lysis_var = false;
    if((float(T0)-float(Tactual)) == 0){
        no_lysis_var = true;
    };
    return no_lysis_var;

};