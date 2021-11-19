#include <unordered_set>
#include <map>  

using namespace std;

struct Vertex {
    int id;
    int weight;
    float score;
};

struct Edge {
   Vertex v1;
   Vertex v2;
};

void program();
void main_logic(unordered_set<Edge> available_edges);
void remove_related_edges(Vertex v);
void step(Vertex v1, Vertex v2);