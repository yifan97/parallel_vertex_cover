#ifndef __MPIVERTEXCOVER_H__
#define __MPIVERTEXCOVER_H__

#include <omp.h>
#include <unordered_set>
#include <map>  
#include <vector>

using namespace std;

struct Vertex {
    int id;
    int weight;
    float score;
    bool isLeaf;
    omp_lock_t lock;
};

struct Edge {
   Vertex* v1;
   Vertex* v2;
   bool isStar;
};

void program();
void main_logic(unordered_set<Edge> available_edges);
void remove_related_edges(Vertex v);
void check_finish(Vertex* v, vector<Vertex*> &toBeDeleted);

void delete_finish(Vertex* v, vector<Vertex*> &toBeDeleted, vector<Vertex*> &tmp);
void run_leaf(Vertex* v);
void run_root(Vertex* v);
vector<Edge*> get_star_edge(Vertex* v);
void remove_related_edges(Vertex* v);
bool istep(Vertex* v, Edge* edge);
bool step(Vertex* v, Edge* edge);
void write_cover_to_file();
bool root_leaf_edge(Edge* edge);
void check_correctness();
void assign_role(Vertex* v);
void print_info(Vertex* v);

bool compare_float(float x, float y);


void take_share(ifstream file, int procID, int nproc);
#endif