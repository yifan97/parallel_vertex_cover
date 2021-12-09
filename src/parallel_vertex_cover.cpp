#include "parallel_vertex_cover.h"

#include <unordered_set>
#include <map>  
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>      
#include <stdlib.h>     
#include <time.h>
#include <string.h> 
#include <algorithm>
#include <vector>
#include <set>
#include <chrono>
#include <omp.h>
#include <atomic>

using namespace std;

static int _argc;
static const char **_argv;

// no lock needed
map<int, Vertex*> id_to_vertex;
set<Edge*> all_edges;

// lock needed
vector<Vertex*> available_vertex;
set<Vertex*> covered_vertex;
map<Vertex*, set<Edge*>> vertex_to_edges;
vector<Vertex*> toBeDeleted;
vector<Vertex*> tmp;


// for correctness check only
map<int, set<Edge*>> final_neighbors;


omp_lock_t available_lock;
omp_lock_t cover_lock;
atomic<int> atomic_count{0};

const char *get_option_string(const char *option_name,
                              const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

int get_option_int(const char *option_name, int default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

static void show_help(const char *program_path)
{
  printf("Usage: %s OPTIONS\n", program_path);
  printf("\n");
  printf("OPTIONS:\n");
  printf("\t-f <input_filename> (required)\n");
  printf("\t-n <num_of_threads> (required)\n");
}

int num_vertex;
int total_weight = 0;

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
    
    int num_of_threads = get_option_int("-n", 1);

    int already_handled = 0;
    std::ifstream file(input_filename);
    if (file.is_open()) {
        string line;
        getline(file, line);
        int i = line.find(" ");
        num_vertex = stoi(line.substr(0, i).c_str());
        int num_edges = stoi(line.substr(i+1).c_str());
        omp_init_lock(&available_lock);
        omp_init_lock(&cover_lock);
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
                total_weight += weight;
                omp_init_lock(&(v->lock));
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
                final_neighbors[edge->v1->id].insert(edge);
                final_neighbors[edge->v2->id].insert(edge);
            }
            already_handled++;
        }
        file.close();
    }

    omp_set_num_threads(num_of_threads);
    cout << "!!!number of threads " << num_of_threads <<"\n";
    init_time += duration_cast<dsec>(Clock::now() - init_start).count();
    printf("Initialization Time: %lf.\n", init_time);

    auto compute_start = Clock::now();
    double compute_time = 0;

    // printf("num of core %d\n",num_of_threads);
    int prev_size = -1;
    while(available_vertex.size() > 0){
        if(available_vertex.size() == prev_size){
            cout << "here !" << endl;
            // break;
        }
        prev_size = available_vertex.size();
        // cout << "???"<<available_vertex.size() << endl;
        int i;
#pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            check_finish(available_vertex[i], toBeDeleted);
        }
        // cout << "after check finish\n";
#pragma omp parallel for shared(available_vertex, toBeDeleted) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            delete_finish(available_vertex[i], toBeDeleted, tmp);
        }
        toBeDeleted.clear();
        available_vertex = tmp;
        tmp.clear();
//  set leaf/root role
// #pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
//         for(i = 0; i < available_vertex.size(); i++){
//             assign_role(available_vertex[i]);
//         }
        //  need to set all star edges first
#pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            run_leaf(available_vertex[i]);
        }
#pragma omp parallel for shared(available_vertex) private(i) schedule(dynamic)
        for(i = 0; i < available_vertex.size(); i++){
            run_root(available_vertex[i]);
        }
    }

    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", compute_time);

    // write_cover_to_file();
    check_correctness();
    return 0;
}


void delete_finish(Vertex* v, vector<Vertex*> &toBeDeleted, vector<Vertex*> &tmp){
    vector<Vertex*>::iterator it = find(toBeDeleted.begin(), toBeDeleted.end(), v);
    if(it == toBeDeleted.end()){  //  put it back if not to be deleted
        omp_set_lock(&available_lock);
        tmp.push_back(v);
        omp_unset_lock(&available_lock);
    }
}

void print_info(Vertex* v){
    cout <<"current id: " << v->id << "\nprinting available vertex\n";
    for(Vertex* i: available_vertex){
        cout << "id: " << i->id;
    }
    cout <<"end printing\n";
}

void check_finish(Vertex* v, vector<Vertex*> &toBeDeleted){
    // cout<< "before if\n";
    if(vertex_to_edges[v].size() == 0){
        omp_set_lock(&available_lock);
        // vector<Vertex*>::iterator it = find(available_vertex.begin(), available_vertex.end(), v);
        // if(it != available_vertex.end()){
        //     available_vertex.erase(it);
        // }else{
        //     cout << "didn't find finished vertex\n";
        //     print_info(v);
        // }
        // cout << "before push\n";
        toBeDeleted.push_back(v);
        // cout << "after push\n";
        omp_unset_lock(&available_lock);
    }
    return;
}

// void assign_role(Vertex* v){
//     int head = rand() % 2;
//     if(head == 0){  //  leaf node
//         v->isLeaf = true;
//     }else{
//         v->isLeaf = false;
//     }
// }

void run_leaf(Vertex* v){
    int head = rand() % 2;
    if(head == 0){  //  leaf node
        v->isLeaf = true;
    }else{
        v->isLeaf = false;
    }
    if(v->isLeaf){  //  leaf node
        vector<Edge*> active_edge;
        for(set<Edge*>::iterator it = vertex_to_edges[v].begin(); it != vertex_to_edges[v].end(); it++){
            if((*it)->isStar || !root_leaf_edge(*it)) {
                continue;
            }
            // if((*it)->isStar || istep(v, *it)){  // active
            //     active_edge.push_back(*it);
            // }
            if(istep(v, *it)){  // active
                active_edge.push_back(*it);
            }
        }
        if (active_edge.size() == 0) {
            return;
        }
        int starIdx = rand() % active_edge.size();
        Edge* starEdge = active_edge[starIdx];
        starEdge->isStar = true;
    }else{
        return;
    }
}

bool root_leaf_edge(Edge* edge){
    return (edge->v1->isLeaf && !edge->v2->isLeaf)||(edge->v2->isLeaf && !edge->v1->isLeaf);
}

void run_root(Vertex* v){
    if(!v->isLeaf){  //  root node
        int runLast = rand()%2;
        vector<Edge*> starEdges = get_star_edge(v);
        if(starEdges.size() == 0) {
            return;
        }
        if(runLast){  //  only run last edge
            int stepIdx = starEdges.size()-1;
            for(int i = 0; i < starEdges.size(); i++){
                if(istep(v, starEdges[i])){
                    stepIdx = i;
                    break;
                }
            }
            step(v, starEdges[stepIdx]);
        }else{
            for(int i = 0; i < starEdges.size(); i++){
                if(step(v, starEdges[i])){
                    break;
                }
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
        // (*it)->isStar = false;
    }
    return starEdges;
}

bool compare_float(float x, float y){
    float epsilon = 0.01f;
    if(fabs(x - y) < epsilon){
        return true; //they are same
    }
    return false; //they are not same
}

bool istep(Vertex* v, Edge* edge){
    Vertex* v1 = edge->v1;
    Vertex* v2 = edge->v2;
    float beta = min((1.0 - v1->score) * v1->weight, (1.0 - v2->score) * v2->weight);
    float score = 0.0;

    if(v1->id == v->id){
        score = v1->score + beta /  v1->weight;
    }else if(v2->id == v->id) {
        score = v2->score + beta /  v2->weight;
    }

    if(compare_float(score, 1.0)){
        return true;
    }else{
        return false;
    }
}

bool step(Vertex* v, Edge* edge){
    Vertex* v1 = edge->v1;
    Vertex* v2 = edge->v2;
    float beta = min((1 - v1->score) * v1->weight, (1 - v2->score) * v2->weight);
    v1->score += beta /  v1->weight;
    v2->score += beta /  v2->weight;
    if(compare_float(v1->score,1.0)) {
        omp_set_lock(&cover_lock);
        covered_vertex.insert(v1);
        remove_related_edges(v1);
        omp_unset_lock(&cover_lock);
        // omp_set_lock(&v1->lock);
        // remove_related_edges(v1);
        // omp_unset_lock(&v1->lock);
        return v1==v;
    }else if(compare_float(v2->score,1.0)) {
        omp_set_lock(&cover_lock);
        covered_vertex.insert(v2);
        remove_related_edges(v2);
        omp_unset_lock(&cover_lock);
        // omp_set_lock(&v2->lock);
        // remove_related_edges(v2);
        // omp_unset_lock(&v2->lock);
        return v2==v;
    }
}

void remove_related_edges(Vertex* v){
    set<Edge*> neighbors = vertex_to_edges[v];
    if(neighbors.size() == 0) return;
    for(set<Edge*>::iterator it = neighbors.begin(); it != neighbors.end(); it++){
        Vertex* v1 = (*it)->v1;
        Vertex* v2 = (*it)->v2;
        set<Edge*>::iterator itr;
        if(v1->id != v->id) {
            omp_set_lock(&v1->lock);
            itr = vertex_to_edges[v1].find(*it);
            if(itr != vertex_to_edges[v1].end()){
                vertex_to_edges[v1].erase(itr);
            }
            omp_unset_lock(&v1->lock);
        } else{
            omp_set_lock(&v2->lock);
            itr = vertex_to_edges[v2].find(*it);
            if(itr != vertex_to_edges[v2].end()){
                vertex_to_edges[v2].erase(itr);
            }
            omp_unset_lock(&v2->lock);
        }
    }
    // omp_set_lock(&v->lock);
    // vertex_to_edges[v].clear();
    // omp_unset_lock(&v->lock);

    vertex_to_edges[v].clear();
    
    omp_set_lock(&available_lock);
    vector<Vertex*>::iterator idx = find(available_vertex.begin(), available_vertex.end(), v);
    if (idx != available_vertex.end()) {
        available_vertex.erase(idx);
    }else{
        cout << "in remove no remove vertex\n";
    }
    omp_unset_lock(&available_lock);
}

void write_cover_to_file(){
    cout << "total number of covered vertex is " << covered_vertex.size() << endl;
    for(set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        cout << "id is " << (*it)->id << endl;
    }
}

void check_correctness(){
    cout << "Now checking correctness!"<< endl;
    set<int> final_cover;
    for(int i = 0; i < num_vertex; i++){
        if(covered_vertex.find(id_to_vertex[i]) == covered_vertex.end()){  // take complement
            final_cover.insert(i);
        }
    }
    cout << "before for1\n";
    // for(set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
    //     // cout << "before id\n";
    //     // cout << "id: " << (*it)->id << endl;
    //     // cout << "after id\n";
    //     final_cover.insert((*it)->id);
    //     // cout << "after insert\n";
    // }
    cout << "after for1\n";
    set<int> correctly_covered;
    int covered_weight = 0;
        for(set<int>::iterator it = final_cover.begin(); it != final_cover.end(); it++){
        int id = *it;
        covered_weight += id_to_vertex[id]->weight;
        map<int, set<Edge*>>::iterator t = final_neighbors.find(id);
        if(t == final_neighbors.end()) {
            continue;
        }
        set<Edge*> neighbors = t->second;
        // cout << "before if\n";
        for(set<Edge*>::iterator itr = neighbors.begin(); itr != neighbors.end(); itr++) {
            correctly_covered.insert((*itr)->v1->id);
            correctly_covered.insert((*itr)->v2->id);
        }
    }
    cout << "Covered " << final_cover.size() << " vertices out of " << num_vertex <<endl;
    cout << "Covered weight " << covered_weight << " out of " << total_weight <<endl;
    cout << "Average weight per node " << (double)total_weight/num_vertex;
    cout << " vs. ours average weight per node " << (double)covered_weight/final_cover.size() << endl;
    cout << "Average total covered weight for "<<final_cover.size() <<" vertices is " << (total_weight/num_vertex)*final_cover.size() << endl;
    cout << "You missed " << num_vertex- correctly_covered.size() << " vertex!!!!" << endl;
}