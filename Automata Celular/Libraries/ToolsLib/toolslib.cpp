#include "toolslib.h"
#include <random>
#include <chrono>
#include <array>
#include <iostream>
#include <algorithm>

int * get_mat(std::string filename, int matlen){
    int * mat;
    int loop = 0;
    mat = new int[matlen];
    std::ifstream file;
    file.open(filename);
    if(file.is_open()){
        while(file >> mat[loop]){
            loop++;
        }
        file.close();
    }
    else{
        std::cout<<"No se pudo abrir el fichero.";
    }
    return mat;
};

void node_to_coordinates(int node, int &x, int &y){
    if(node>NODE_NUM){
        std::cout<<"La matriz es demasiado pequena para albergar ese nodo.";
    }
    y = int(node/X_SIZE);
    x = node - y*X_SIZE;

};

void coordinates_to_node(int &node, int x, int y){
    node = ((2*NX-1)-1)*y+(y+x);
}

std::vector<int> get_neighbours(int *mat, int node){
    // Devuelve los vecinos de un nodo
    int x,y,right_node, left_node, up_node,down_node, diag_ur_node,diag_ul_node, diag_dr_node, diag_dl_node;
    node_to_coordinates(node, x, y);
    up_node = node-1;
    down_node = node+1;
    coordinates_to_node(right_node,x,y-1);
    coordinates_to_node(left_node,x,y+1);
    coordinates_to_node(diag_ur_node,x-1,y+1);
    coordinates_to_node(diag_ul_node,x-1,y-1);
    coordinates_to_node(diag_dr_node,x+1,y+1);
    coordinates_to_node(diag_dl_node,x+1,y-1);
    std::vector<int> neighbour_nodes{right_node,left_node,
                                    up_node,down_node,
                                    diag_ur_node,diag_ul_node, 
                                    diag_dr_node, diag_dl_node};
    neighbour_nodes.erase(std::remove_if(neighbour_nodes.begin(), neighbour_nodes.end(), is_negative), neighbour_nodes.end());
    neighbour_nodes.erase(std::remove_if(neighbour_nodes.begin(), neighbour_nodes.end(), out_of_mat), neighbour_nodes.end());

    return neighbour_nodes;
}

bool out_of_mat(int value){
    return (value>NODE_NUM);
}

bool is_negative(int value){
    return (value<=0);
}

std::vector<int> get_specific_neighbours(int *mat, int node, int value, char mode){
    // gets the n neighbours with a specific value and a operator specified
    std::vector<int> neighbour_nodes, specific_neighbour_nodes;
    neighbour_nodes = get_neighbours(mat, node);

    for(int i = 0; i<neighbour_nodes.size(); i++){
        switch (mode)
        {
        case '=':
            if((mat[neighbour_nodes[i]] == value)){
                specific_neighbour_nodes.push_back(neighbour_nodes[i]);
            }
            break;
        
        case '>':
            if(mat[neighbour_nodes[i]] > value){
                specific_neighbour_nodes.push_back(neighbour_nodes[i]);
            }
            break;
        case '<':
            if(mat[neighbour_nodes[i]] < value){
                specific_neighbour_nodes.push_back(neighbour_nodes[i]);
            }
        }
    }

    return specific_neighbour_nodes;
};

void create_vec(int node_num, int mat[], int value){ //poner un mensaje si node_num no es correcto
    for (int i = 0; i<node_num; i++) 
    {
        mat[i] = value;
    }
}

void create_vec(int node_num, float mat[], float value){
    for (int i = 0; i<node_num; i++) 
    {
        mat[i] = value;
    }
};

void save_mat(int node_num, int mat[], std::string filename){
    std::ofstream File(filename); 
    for(int node = 0; node<node_num; node++){
        File << mat[node] << std::endl;
    }
    File.close();
}

void save_mat(int node_num, float mat[], std::string filename){ //quitar y usar solo double
    std::ofstream File(filename); 
    for(int node = 0; node<node_num; node++){
        File << mat[node] << std::endl;
    }
    File.close();
}

void save_mat(int node_num, double mat[], std::string filename){
    std::ofstream File(filename); 
    for(int node = 0; node<node_num; node++){
        File << mat[node] << std::endl;
    }
    File.close();
}

void save_vec(std::vector<int> mat, std::string filename){
    std::ofstream File(filename); 
    for(int node = 0; node<mat.size(); node++){
        File << mat[node] << std::endl;
    }
    File.close();
}

void changeNegativeValue(double &value){
    if(value<0){
        value = 0;
    }
    return;
}

int* get_random_nodes(int xsize, int ysize){
    int* random_nodes;
    random_nodes = new int[xsize*ysize];

    for(int node = 0; node < (xsize*ysize); node++){
        random_nodes[node] = node;
    }
    std::random_shuffle(&random_nodes[0], &random_nodes[xsize*ysize]);
    return random_nodes;
}

bool metastasis(int *mat, int xsize, int ysize){ //cambiar nombre
    bool result=false, loop = true;
    int node = 0;
    while(loop){
        for(int x = 0; x<xsize; x++){
            for(int y = 0; y<ysize; y++){
                if(mat[node]>=MAX_PILED_CELL){
                    result = true;
                    loop = false;
                }
                if( (x==0) || (y==0) || (x==(xsize-1)) || (y==(ysize-1))){
                    if(mat[node] != 0){
                        result = true;
                        loop = false;
                    }
                }
                node++;
            }
        }
        loop = false;
    }
    return result;
};

double* int_2_double(int mat[], int matlen){
    double* matd;
    matd = new double[NODE_NUM];
    for(int node=0; node<matlen; node++){
        matd[node] = (double)mat[node];
    }
    return matd;
};

int cell_counter(int mat[]){
    int n_cells = 0;
    for(int node = 0;node<NODE_NUM; node++){
        n_cells = n_cells + mat[node];
        if(n_cells<0){
            std::cout<<"Hay mas nodos de los que deberia."<<std::endl;
        }
    };
    return n_cells;
};

void get_occupied_nodes(int mat[], int mat_nodes[]){
    int i=0;
    for(int node = 0; node<NODE_NUM; node++){
        if(mat[node]>0){
            mat_nodes[i] = node;
            i++;
        }
    }
};

void filename_inc ( std::string *filename )

//  Purpose: increments a partially numeric file name.

{
  char c;
  int change;
  int i;
  int lens;
  std::string numstr;
  int num;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    std::cout << "\n";
    std::cout << "FILENAME_INC - Fatal error!\n";
    std::cout << "  The input std::string is empty.\n";
    exit ( 1 );
  }

  change = 0;
  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];
    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}

