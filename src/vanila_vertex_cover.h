#ifndef __VANILAVERTEXCOVER_H__
#define __VANILAVERTEXCOVER_H__

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

void program();
void main_logic(unordered_set<Edge*> available_edges);
void remove_related_edges(Vertex *v);
void write_cover_to_file();

#endif