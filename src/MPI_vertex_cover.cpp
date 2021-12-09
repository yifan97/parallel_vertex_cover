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
#include <cmath>

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


void take_share(ifstream file, int procID, int nproc){
    if (file.is_open()) {
        string line;
        getline(file, line);
        int i = line.find(" ");
        num_vertex = stoi(line.substr(0, i).c_str());
        int width = sqrt(num_vertex);
        int num_edges = stoi(line.substr(i+1).c_str());
        int a = sqrt(nproc);
        int already_handled = 0;
        int rowID, colID, left, right, up, down;
        if(a*a == nproc){  //  1, 4 or 16
            int unitA = width/a;
            rowID = procID/a;
            colID = procID%a;
            left = colID*unitA;
            right = left+unitA;
            up = rowID*unitA;
            down = up+unitA;
        }else{  //  8->2*4
            int unitA = width/4;
            rowID = procID/4;
            colID = procID%4;
            left = colID*unitA;
            right = left+unitA;
            up = rowID*2*unitA;
            down = up+2*unitA;
        }
        while (getline(file, line)) {
            i = line.find(" ");
            if (already_handled < num_vertex) {  //  store vertex
                int id = stoi(line.substr(0, i).c_str());
                int weight = stoi(line.substr(i+1).c_str());
                int row = id/a;
                int col = id%a;
                if(row >= up && row < down && col >= left && col < right){
                    Vertex* v = (Vertex*)malloc(sizeof(struct Vertex));
                    v->id = id;
                    v->weight = weight;
                    v->score = 0.0;
                    v->isLeaf =false;
                    id_to_vertex[id] = v;
                    available_vertex.push_back(v);
                }
            }else{  //  store edges
                int id1 = stoi(line.substr(0, i).c_str());
                int id2 = stoi(line.substr(i+1).c_str());
                int row1 = id1/a;
                int col1 = id1%a;
                int row2 = id2/a;
                int col2 = id2%a;
                if(row1 >= up && row1 < down && col1 >= left && col1 < right){
                    if(row2 >= up && row2 < down && col2 >= left && col2 < right){
                        Edge* edge = (Edge*)malloc(sizeof(struct Edge));
                        edge->v1 = id_to_vertex[id1];
                        edge->v2 = id_to_vertex[id2];
                        edge->isStar = false;
                        all_edges.insert(edge);
                        vertex_to_edges[edge->v1].insert(edge);
                        vertex_to_edges[edge->v2].insert(edge);       
                    }
                }
            }
            already_handled++;
        }
        file.close();
    }
}

int compute(int procID, int nproc, char *inputFilename){
    
    std::ifstream file(input_filename);
    take_share(file, procID, nproc);  //  each processor take corresponding share of data


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
        for(i = 0; i < available_vertex.size(); i++){
            check_finish(available_vertex[i], toBeDeleted);
        }
        // cout << "after check finish\n";
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
        for(i = 0; i < available_vertex.size(); i++){
            run_leaf(available_vertex[i]);
        }
        for(i = 0; i < available_vertex.size(); i++){
            run_root(available_vertex[i]);
        }
    }

    communicate_cover();
    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", compute_time);

    // write_cover_to_file();
    check_correctness();
    return 0;
}

void communicate_cover(int procID, int nproc){
    if(nproc == 1) return;
    MPI_Request send_request1;
    //  first send number of covered vertices
    //  then send the actual Ids and dedup
    if(procID == 0){
        int total = covered_vertex.size();
        int* vertex_count = (int*)calloc(nproc, sizeof(int));
        for(int i = 1; i < nproc; i++){
            receive_data(vertex_count+i, 1, i);
            total += vertex_count[i];
        }
        int* buffer = (int*)calloc(total, sizeof(int));
        int offset = covered_vertex.size();
        for(int i = 1; i < nproc; i++){
            receive_data(buffer+offset, vertex_count[i], i);
            offset += vertex_count[i];
        }
    }else{
        int* count = (int*)calloc(1, sizeof(int));
        count[0] = covered_vertex.size();
        start_send_data(count, 1, 0, &send_request1);
        finish_data(&send_request1);
        int* data = (int*)calloc(count[0], sizeof(int));
        int idx = 0;
        for(Vertex* v: covered_vertex){
            data[idx++] = v->id;
        }
        start_send_data(data, count[0], 0, &send_request1);
        finish_data(&send_request1);
    }
    
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
    vector<int> covered_ids;
    cout << "before for1\n";
    for(set<Vertex*>:: iterator it = covered_vertex.begin(); it != covered_vertex.end(); it++){
        // cout << "before id\n";
        // cout << "id: " << (*it)->id << endl;
        // cout << "after id\n";
        covered_ids.push_back((*it)->id);
        // cout << "after insert\n";
    }
    int count = 0;
    cout << "after for1\n";
    for(set<Edge*>:: iterator it = all_edges.begin(); it != all_edges.end(); it++){
        int id1 = (*it)->v1->id;
        int id2 = (*it)->v2->id;
        // cout << "before if\n";
        if(find(covered_ids.begin(), covered_ids.end() ,id1) == covered_ids.end() && find(covered_ids.begin(), covered_ids.end() ,id2) == covered_ids.end()) {
            cout << "Miss edge" << id1 << " , " << id2 <<endl;
            count++;
        }
    }
    cout << "Covered " << covered_vertex.size() << " vertices out of " << num_vertex <<endl;
    cout << "You missed " << count << " vertex!!!!" << endl;
}

void receive_data(int *data, int count, int process_id)
{
  MPI_Recv(data, count, MPI_INT, process_id,
           0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void start_send_data(int *data, int count, int process_id, MPI_Request *request)
{
  MPI_Isend(data, count, MPI_INT, process_id,
            0, MPI_COMM_WORLD, request);
}

void finish_data(MPI_Request *request)
{
  MPI_Wait(request, MPI_STATUS_IGNORE);
}