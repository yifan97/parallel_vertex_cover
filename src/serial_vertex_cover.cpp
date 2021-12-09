#include <unordered_set>
#include <map>  
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <string.h> 
#include <unordered_set>
#include <algorithm>   
#include <chrono>
#include "serial_vertex_cover.h"

using namespace std;

static int _argc;
static const char **_argv;

unordered_set<Edge*> available_edges;
unordered_set<Vertex*> covered_vertex;
map<int, Vertex*> id_to_vertex;
map<Vertex*, unordered_set<Edge*> > vertex_to_edges;
unordered_set<Edge*> all_edges;

const char *get_option_string(const char *option_name,
                              const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

static void show_help(const char *program_path)
{
  printf("Usage: %s OPTIONS\n", program_path);
  printf("\n");
  printf("OPTIONS:\n");
  printf("\t-f <input_filename> (required)\n");
}

int num_vertex;

int main(int argc, const char *argv[]){
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;

    auto init_start = Clock::now();
    double init_time = 0;

    _argc = argc - 1;
    _argv = argv + 1;

    const char *input_filename = get_option_string("-f", "../data/graph_data");
    if (input_filename == NULL) {
        printf("Input file is required! Check the instruction and run it again!\n");
        exit(-1);
    } 

    int already_handled = 0;
    std::ifstream file(input_filename);
    if (file.is_open()) {
        string line;
        getline(file, line);
        int i = line.find(" ");
        num_vertex = stoi(line.substr(0, i).c_str());
        int num_edges = stoi(line.substr(i+1).c_str());
        while (getline(file, line)) {
            i = line.find(" ");
            if (already_handled < num_vertex) {
                int id = stoi(line.substr(0, i).c_str());
                int weight = stoi(line.substr(i+1).c_str());
                Vertex *v = (Vertex*) malloc(sizeof(struct Vertex));
                v->id = id;
                v->weight = weight;
                v->score = 0.0;
                id_to_vertex[id] = v;
            } else{
                int id1 = stoi(line.substr(0, i).c_str());
                int id2 = stoi(line.substr(i+1).c_str());
                Edge *edge = (Edge*) malloc(sizeof(Vertex) * 2);
                edge->v1 = id_to_vertex[id1];
                edge->v2 = id_to_vertex[id2];
                available_edges.insert(edge);
                all_edges.insert(edge);
                vertex_to_edges[edge->v1].insert(edge);
                vertex_to_edges[edge->v2].insert(edge);
            }
            already_handled++;
        }
        file.close();
    }

    init_time += duration_cast<dsec>(Clock::now() - init_start).count();
    printf("Initialization Time: %lf.\n", init_time);

    auto compute_start = Clock::now();
    double compute_time = 0;

    while(available_edges.size() > 0){
        main_logic(available_edges);
    }

    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", compute_time);

    // write_cover_to_file();
    check_correctness();
    return 0;
}

void main_logic(unordered_set<Edge*> available_edges){
    int random_idx = rand() % available_edges.size();
    unordered_set<Edge*>::iterator it = available_edges.begin();;
    for(int i = 0; i < random_idx; i++){
        it++;
    }
    
    step(*it);
}

void remove_related_edges(Vertex *v){
    unordered_set<Edge*> neighbors = vertex_to_edges[v];
    for(unordered_set<Edge*>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
        Vertex *v1 = (*it)->v1;
        Vertex *v2 = (*it)->v2;
        
        if(v1->id != v->id) {
            vertex_to_edges[v1].erase(*it);
        } else{
            vertex_to_edges[v2].erase(*it);
        }
        available_edges.erase(*it);
    }

    vertex_to_edges[v].clear();
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

    // remove the edge 
    // available_edges.erase(edge);
}

void write_cover_to_file(){
    cout << "total number of covered vertex is " << covered_vertex.size() << endl;
    for(unordered_set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        cout << "id is " << (*it)->id << endl;
    }
}

void check_correctness(){
    unordered_set<int> covered_ids;
    for(unordered_set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        covered_ids.insert((*it)->id);
    }
    int count = 0;
    for(unordered_set<Edge*>:: iterator it = all_edges.begin(); it != all_edges.end(); it++){
        int id1 = (*it)->v1->id;
        int id2 = (*it)->v2->id;

        if(covered_ids.find(id1) == covered_ids.end() && covered_ids.find(id2) == covered_ids.end()) {
            count++;
        }
    }
    cout << "Covered " << covered_vertex.size() << " vertices out of " << num_vertex <<endl;
    cout << "You missed " << count << " vertex!!!!" << endl;
}