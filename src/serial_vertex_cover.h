#ifndef __SERIALVERTEXCOVER_H__
#define __SERIALVERTEXCOVER_H__

#include <unordered_set>
#include <map>  

using namespace std;

struct Vertex {
    int id;
    int weight;
    float score;
};

struct Edge {
   Vertex *v1;
   Vertex *v2;
};

void main_logic(unordered_set<Edge*> available_edges);
void remove_related_edges(Vertex *v);
void step(Edge* edge);
void write_cover_to_file();
void check_correctness();

#endif