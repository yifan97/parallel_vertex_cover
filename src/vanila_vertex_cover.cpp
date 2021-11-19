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
#include "serial_vertex_cover.h"

using namespace std;

unordered_set<Edge> available_edges;
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

    while(available_edges.size() >= 0){
        main_logic(available_edges);
    }
}

// While a edge is not covered:
//         For (v, w) in available_edges:
// 		step(v ,w)
//   		If v == 1:
//             		add_to_cover(v)
// 			remove all ages v connects
//         	else if w == 1:
//             		add_to_cover(w)
// 			remove all ages w connects

void main_logic(unordered_set<Edge> available_edges){
    int random_idx = rand() % available_edges.size();
    unordered_set<Edge>::iterator it = available_edges.begin();
    for(int i = 0; i < random_idx; i++){
        it++;
    }
    Vertex v1 = (*it).v1;
    Vertex v2 = (*it).v2;
    covered_vertex.insert(v1);
    remove_related_edges(v1);
    covered_vertex.insert(v2);
    remove_related_edges(v2);
}

void remove_related_edges(Vertex v){
    unordered_set<Edge> neighbors = vertex_to_edges[v];
    for(auto it = neighbors.begin(); it != neighbors.end(); it++){
        Vertex v1 = (*it).v1;
        Vertex v2 = (*it).v2;

        if(v1.id != v.id) {
            vertex_to_edges[v1].erase(it);
        } else{
            vertex_to_edges[v2].erase(it);
        }
    }
    vertex_to_edges[v].clear();
}
