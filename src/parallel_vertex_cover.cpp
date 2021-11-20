#include "parallel_vertex_cover.h"

#include <unordered_set>
#include <map>  
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <algorithm>
#include <vector>
#include <set>
#include <omp.h>

using namespace std;

vector<Vertex*> available_vertex;
set<Vertex*> covered_vertex;
map<int, Vertex*> id_to_vertex;
map<Vertex*, set<Edge*>> vertex_to_edges;
unordered_set<Edge*> all_edges;

int main(int argc, const char *argv[]){
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
                Vertex* v = (Vertex*)malloc(sizeof(struct Vertex));
                v->id = id;
                v->weight = weight;
                v->score = 0.0;
                v->isLeaf =false;
                id_to_vertex[id] = v;
                available_vertex.push_back(v);
            } else{
                int id1 = stoi(line.substr(0, i).c_str());
                int id2 = stoi(line.substr(i+1).c_str());
                Edge* edge = (Edge*)malloc(sizeof(struct Edge));
                edge->v1 = id_to_vertex[id1];
                edge->v2 = id_to_vertex[id2];
                edge->isStar = false;
                all_edges.insert(edge);
                vertex_to_edges[edge->v1].insert(edge);
                vertex_to_edges[edge->v2].insert(edge);
            }
            already_handled++;
        }
        file.close();
    }

    omp_set_num_threads(1);

    while(available_vertex.size() > 0){
        cout<< "available vertex size " << available_vertex.size() <<endl;
        int i;
        vector<Vertex*> toBeDeleted;
#pragma omp parallel for shared(available_vertex, toBeDeleted) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            check_finish(available_vertex[i], toBeDeleted);
        }
        if(available_vertex.size() == toBeDeleted.size()){
            cout<< "break" <<endl;
            break;
        }
        cout << "after check finish" << endl;
        vector<Vertex*> tmp;
#pragma omp parallel for shared(available_vertex, toBeDeleted) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            delete_finish(available_vertex[i], toBeDeleted, tmp);
        }
        toBeDeleted.clear();
        available_vertex = tmp;
        
        cout<< "ok 1" <<endl;

        //  need to set all star edges first
#pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            run_leaf(available_vertex[i]);
        }

        cout<< "ok 2" <<endl;

#pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            run_root(available_vertex[i]);
        }
    }

    write_cover_to_file();
    check_correctness();
    return 0;
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

void delete_finish(Vertex* v, vector<Vertex*> toBeDeleted, vector<Vertex*> tmp){
    vector<Vertex*>::iterator it = find(toBeDeleted.begin(), toBeDeleted.end(), v);
    if(it != toBeDeleted.end()){
        tmp.push_back(*it);
    }
}

void check_finish(Vertex* v, vector<Vertex*> toBeDeleted){
    cout << "in check" <<endl;
    if(vertex_to_edges[v].size() == 0){
        cout << "in if" <<endl;
        // omp_set_lock(&v->lock);
        vector<Vertex*>::iterator it = find(available_vertex.begin(), available_vertex.end(), v);
        if(it != available_vertex.end()){
            toBeDeleted.push_back(*it);
            covered_vertex.insert(v);
            cout<<"insert finish"<<endl;
        }
        // omp_unset_lock(&v->lock);
    }
    cout << "out check" << endl;
    return;
}

void run_leaf(Vertex* v){
    int head = rand() % 2;
    if(head == 0){  //  leaf node
        v->isLeaf = true;
        vector<Edge*> active_edge;
        for(set<Edge*>::iterator it = vertex_to_edges[v].begin(); it != vertex_to_edges[v].end(); it++){
            if(istep(v, *it)){  // active
                active_edge.push_back(*it);
            }
        }
        cout << "active edge size " << active_edge.size() << endl;
        if (active_edge.size() == 0) {
            return;
        }

        int starIdx = rand() % active_edge.size();
        Edge* starEdge = active_edge[starIdx];
        starEdge->isStar = true;
    }else{
        v->isLeaf = false;
        return;
    }
}

void run_root(Vertex* v){
    if(!v->isLeaf){  //  root node
        int runLast = rand()%2;
        vector<Edge*> starEdges = get_star_edge(v);
        if(starEdges.size() == 0) {
            return;
        }
        if(runLast){  //  only run last edge
            step(starEdges[starEdges.size()-1]);
        }else{
            for(int i = 0; i < starEdges.size(); i++){
                step(starEdges[i]);
            }
        }
    }else{
        v->isLeaf = false;
    }
}

vector<Edge*> get_star_edge(Vertex* v){
    vector<Edge*> starEdges;
    set<Edge*> neighbors = vertex_to_edges[v];
    for(set<Edge*>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
        if((*it)->isStar){
            starEdges.push_back(*it);
        }
        (*it)->isStar = false;
    }
    return starEdges;
}

void remove_related_edges(Vertex* v){
    set<Edge*> neighbors = vertex_to_edges[v];
    for(set<Edge*>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
        Vertex* v1 = (*it)->v1;
        Vertex* v2 = (*it)->v2;
        if(v1->id != v->id) {
            vertex_to_edges[v1].erase(*it);
        } else{
            vertex_to_edges[v2].erase(*it);
        }
    }
    vertex_to_edges[v].clear();
    vector<Vertex*>::iterator idx = find(available_vertex.begin(), available_vertex.end(), v);
    if (idx != available_vertex.end()) {
        available_vertex.erase(idx);
    }
}

bool istep(Vertex* v, Edge* edge){
    Vertex* v1 = edge->v1;
    Vertex* v2 = edge->v2;
    float beta = min((1 - v1->score) * v1->weight, (1 - v2->score) * v2->weight);
    float score = 0.0;

    if(v1->id == v->id){
        score = v1->score + beta /  v1->weight;
    }else if(v2->id == v->id) {
        score = v2->score + beta /  v2->weight;
    }

    if(score == 1.0){
        return true;
    }else{
        return false;
    }
}

void step(Edge* edge){
    Vertex* v1 = edge->v1;
    Vertex* v2 = edge->v2;
    float beta = min((1 - v1->score) * v1->weight, (1 - v2->score) * v2->weight);
    v1->score += beta /  v1->weight;
    v2->score += beta /  v2->weight;
    if(v1->score == 1.0) {
        covered_vertex.insert(v1);
        remove_related_edges(v1);
    } else if(v2->score == 1.0) {
        covered_vertex.insert(v2);
        remove_related_edges(v2);
    }
}

void write_cover_to_file(){
    cout << "total number of covered vertex is " << covered_vertex.size() << endl;
    for(set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        cout << "id is " << (*it)->id << endl;
    }
}

void check_correctness(){
    unordered_set<int> covered_ids;
    for(set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        covered_ids.insert((*it)->id);
    }
    int count = 0;
    for(unordered_set<Edge*>:: iterator it = all_edges.begin(); it != all_edges.end(); it++){
        int id1 = (*it)->v1->id;
        int id2 = (*it)->v2->id;

        if(covered_ids.find(id1) == covered_ids.end() && covered_ids.find(id2) == covered_ids.end()) {
            // cout << "You missed some vertex!!!!" << endl;
            count++;
        }
    }
    cout << "You missed " << count << endl;
}