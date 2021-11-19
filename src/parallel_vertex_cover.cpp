#include <unordered_set>
#include <map>  
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <unordered_set>
#include <algorithm>
#include <vector>   
#include "serial_vertex_cover.h"
#include <omp.h>

using namespace std;

unordered_set<Edge> available_edges;
vector<Vertex> available_vertex;
unordered_set<Vertex> covered_vertex;
map<int, Vertex> id_to_vertex;
map<Vertex, unordered_set<Edge>> vertex_to_edges;

void program(){
    
    int already_handled = 0;
    std::ifstream file("../data/graph_data");
    if (file.is_open()) {
        string line;
        getline(file, line);
        int i = line.find(" ");
        int num_vertex = stoi(line.substr(0, i).c_str());
        int num_edges = stoi(line.substr(i+1).c_str());
        while (getline(file, line)) {
            i = line.find(" ");
            if (already_handled < num_vertex) {
                int id = stoi(line.substr(0, i).c_str());
                int weight = stoi(line.substr(i+1).c_str());
                Vertex v;
                v.id = id;
                v.weight = weight;
                v.score = 0.0;
                id_to_vertex[id] = v;
            } else{
                int id1 = stoi(line.substr(0, i).c_str());
                int id2 = stoi(line.substr(i+1).c_str());
                Edge edge;
                edge.v1 = id_to_vertex[id1];
                edge.v2 = id_to_vertex[id2];
                available_edges.insert(edge);
                vertex_to_edges[edge.v1].insert(edge);
                vertex_to_edges[edge.v2].insert(edge);
            }
            already_handled++;
        }
        file.close();
    }

    while(available_vertex.size() >= 0){
#pragma omp parallel for shared(wires) private(i) schedule(dynamic)
        for(int i = 0; i < available_vertex.size(); i++){
            main_logic(available_vertex[i]);
        }
    }
}

// While a vertex is not finished:
//    if all edges are covered, finish
//    else flip a coin, heads->leaf node, tails->root node
//    if leaf v:
//      label edge （v,w）active if w is root and step(x,(v,w)) would add v
//      randomly choose one active edge as star edge 
//         For (v, w) in available_edges:
//    if root w:
//      flip a coin
//      if heads:
//          for each star edge(fixed order):
//              if w not yet in the cover, do step(x,(v,w))
//      if tails:
//          do step for only the last edge in above order          


// 		step(v ,w)
//      beta = min((1-Xv)*Cv, (1-Xw)Cw)
//      Xv += beta/Cv, Xw += beta/Cw
//   	If v == 1:
// 		    cover v, and all v's edges
//      else if w == 1:
// 		    cover w, and all w's edges

void main_logic(Vertex v, int i){
    if(vertex_to_edges[v].size() == 0){
        omp_set_lock(&v.lock);
        vector<Vertex>::iterator it = find(available_vertex.begin(), available_vertex.end(), v);
        if(it != available_vertex.end()){
            available_vertex.erase(it);
            covered_vertex.insert(v);
        }
        return;
        omp_unset_lock(&v.lock);
    }
    int head = rand() % 2;
    unordered_set<Edge>::iterator it = available_edges.begin();
    for(int i = 0; i < random_idx; i++){
        it++;
    }
    Vertex v1 = (*it).v1;
    Vertex v2 = (*it).v2;
    step(v1, v2);
}

void remove_related_edges(Vertex v){
    unordered_set<Edge> neighbors = vertex_to_edges[v];
    for(unordered_set<Edge>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
        Vertex v1 = (*it).v1;
        Vertex v2 = (*it).v2;
        if(v1.id != v.id) {
            vertex_to_edges[v1].erase(it);
        } else{
            vertex_to_edges[v2].erase(it);
        }
        available_edges.erase(it);
    }
    vertex_to_edges[v].clear();
}

void step(Vertex v1, Vertex v2){
    float beta = min((1 - v1.score) * v1.weight, (1 - v2.score) * v2.weight);
    v1.score += beta /  v1.weight;
    v2.score += beta /  v2.weight;
    if(v1.score == 1.0) {
        covered_vertex.insert(v1);
        remove_related_edges(v1);
    } else if(v2.score == 1.0) {
        covered_vertex.insert(v2);
        remove_related_edges(v2);
    }
}
