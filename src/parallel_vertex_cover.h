#ifndef __PARALLELVERTEXCOVER_H__
#define __PARALLELVERTEXCOVER_H__

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
void step(Edge* edge);
void check_finish(Vertex* v, int i);
void run_leaf(Vertex* v, int i);
void run_root(Vertex* v, int i);
vector<Edge*> get_star_edge(Vertex* v);
void remove_related_edges(Vertex* v);
bool istep(Vertex* v, Edge* edge);
void write_cover_to_file();
void check_correctness();

#endif