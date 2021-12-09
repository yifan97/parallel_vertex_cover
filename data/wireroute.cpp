/**
 * Parallel VLSI Wire Routing via OpenMP
 * Yifan Xu(yifanxu2), Zhelong Li(zhelongl)
 */

#include "wireroute.h"

#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <algorithm>
#include <string.h>
#include <iostream>
#include <stdio.h>

using namespace std;

static int _argc;
static const char **_argv;

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

float get_option_float(const char *option_name, float default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

static void show_help(const char *program_path)
{
  printf("Usage: %s OPTIONS\n", program_path);
  printf("\n");
  printf("OPTIONS:\n");
  printf("\t-f <input_filename> (required)\n");
  printf("\t-n <num_of_threads> (required)\n");
  printf("\t-p <SA_prob>\n");
  printf("\t-i <SA_iters>\n");
}

int main(int argc, const char *argv[])
{
  using namespace std::chrono;
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> dsec;

  auto init_start = Clock::now();
  double init_time = 0;

  _argc = argc - 1;
  _argv = argv + 1;

  const char *input_filename = get_option_string("-f", NULL);
  int num_of_threads = get_option_int("-n", 1);
  double SA_prob = get_option_float("-p", 0.1f);
  int SA_iters = get_option_int("-i", 5);

  int error = 0;

  if (input_filename == NULL)
  {
    printf("Error: You need to specify -f.\n");
    error = 1;
  }

  if (error)
  {
    show_help(argv[0]);
    return 1;
  }

  printf("Number of threads: %d\n", num_of_threads);
  printf("Probability parameter for simulated annealing: %lf.\n", SA_prob);
  printf("Number of simulated annealing iterations: %d\n", SA_iters);
  printf("Input file: %s\n", input_filename);

  FILE *input = fopen(input_filename, "r");

  if (!input)
  {
    printf("Unable to open file: %s.\n", input_filename);
    return 1;
  }

  int dim_x, dim_y;
  int num_of_wires;
  int i;

  fscanf(input, "%d %d\n", &dim_x, &dim_y);
  fscanf(input, "%d\n", &num_of_wires);

  wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));
  /* Read the grid dimension and wire information from file */

  cost_t *costs = (cost_t *)calloc(1, sizeof(cost_t));
  /* Initialize cost matrix */

  for (i = 0; i < num_of_wires; i++)
  {
    int x1, y1, x2, y2;
    fscanf(input, "%d %d %d %d\n", &x1, &y1, &x2, &y2);
    wires[i].bendCounts = 0;
    wires[i].Endpoints[0] = x1;
    wires[i].Endpoints[1] = y1;
    wires[i].Endpoints[2] = x2;
    wires[i].Endpoints[3] = y2;
  }

  costs->grid = (cell_t *)calloc(dim_x * dim_y, sizeof(cell_t));
  costs->dim_x = dim_x;
  costs->dim_y = dim_y;

  //  initialize all locks
  for (int k = 0; k < dim_x * dim_y; k++)
  {
    omp_init_lock(&(costs->grid[k].lock));
  }

  /* Initailize additional data structures needed in the algorithm */

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  auto compute_start = Clock::now();
  double compute_time = 0;

  /**
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm.
   */

  omp_set_num_threads(num_of_threads);

// first use simple path to start off the iteration
#pragma omp parallel for shared(wires) private(i) schedule(dynamic)
  for (i = 0; i < num_of_wires; i++)
  {
    simple_path(&wires[i]);
    place_wire(&wires[i], costs, i);
  }
  /* In each itration, we need to do the following things
  1. Calculate cost of current path, if not known. This is the current minimum path.
  2. Consider all paths which first travel horizontally. If any cost less than the current minimum
  path, that is the new minimum path.
  3. Consider all paths which first travel vertically. If any cost less than the current minimum
  path, that is the new minimum path.
  4. With probability 1 − P, choose the current minimum path. Otherwise, choose a path
  uniformly at random from the space of ∆x + ∆y possible routes.
  */

  cell_t *grid = costs->grid;

  int y, x;
  for (int itr = 0; itr < SA_iters; itr++)
  {
    //  clean the board
// #pragma omp parallel for private(y, x) shared(grid) schedule(dynamic)
//     for (y = 0; y < dim_y; y++)
//     {
//       for (x = 0; x < dim_x; x++)
//       {
//         grid[y * dim_y + x].cost = 0;      // clean up cost
//         grid[y * dim_y + x].wireCount = 0; // clean up wireCount
//       }
//     }
//     //  place all wires' best path
// #pragma omp parallel for shared(wires) private(i) schedule(dynamic)
//     for (i = 0; i < num_of_wires; i++)
//     {
//       place_wire(&wires[i], costs, i);
//     }

    //  get best path's cost for each wire
// #pragma omp parallel for shared(wires) private(i) schedule(dynamic)
//     for (i = 0; i < num_of_wires; i++)
//     {
//       collect_wire_cost(&wires[i], costs, i);
//     }

    // run algorithm
#pragma omp parallel for shared(wires, costs) private(i) schedule(dynamic)
    for (i = 0; i < num_of_wires; i++)
    {
      srand(time(NULL));
      if ((SA_prob * 100) <= (rand() % 100))
      {
        // find all paths
        delete_wire(&wires[i], costs, i);
        all_paths(&wires[i], costs, i);
        place_wire(&wires[i], costs, i);
      }
      else
      {
        // use random path
        delete_wire(&wires[i], costs, i);
        random_path(&wires[i]);
        place_wire(&wires[i], costs, i);
      }
    }
  }

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);

  /* Write wires and cost_t to files */
  string inputFile;
  inputFile.append(input_filename);
  size_t start = inputFile.find_last_of("/\\");
  size_t end = inputFile.find_last_of(".");

  char routes_file_name[100];
  memset(routes_file_name, 0, sizeof(routes_file_name));
  strcat(routes_file_name, "output_");
  strcat(routes_file_name, inputFile.substr(start + 1, end - start - 1).c_str());
  strcat(routes_file_name, "_");
  strcat(routes_file_name, std::to_string(num_of_threads).c_str());
  strcat(routes_file_name, ".txt");

  char costs_file_name[100];
  memset(costs_file_name, 0, sizeof(costs_file_name));
  strcat(costs_file_name, "costs_");
  strcat(costs_file_name, inputFile.substr(start + 1, end - start - 1).c_str());
  strcat(costs_file_name, "_");
  strcat(costs_file_name, std::to_string(num_of_threads).c_str());
  strcat(costs_file_name, ".txt");

  printf("routes file name: %s.\n", routes_file_name);
  printf("costs file name: %s.\n", costs_file_name);

  FILE *routes_file, *costs_file;

  int maxCost = 0;
  // write cost_t
  costs_file = fopen(costs_file_name, "w");
  fprintf(costs_file, "%d %d\n", dim_x, dim_y);
  for (int r = 0; r < dim_y; r++)
  {
    for (int c = 0; c < dim_x; c++)
    {
      maxCost = max(maxCost, costs->grid[r * dim_x + c].cost);
      fprintf(costs_file, "%d ", costs->grid[r * dim_x + c].cost);
    }
    fprintf(costs_file, "\n");
  }
  printf("Cost: %d\n", maxCost);

  // write routes
  routes_file = fopen(routes_file_name, "w");
  fprintf(routes_file, "%d %d\n", dim_x, dim_y);
  fprintf(routes_file, "%d\n", num_of_wires);
  for (int w = 0; w < num_of_wires; w++)
  {
    fprintf(routes_file, "%d %d ", wires[w].Endpoints[0],
            wires[w].Endpoints[1]);

    if (wires[w].bendCounts == 1)
    {
      fprintf(routes_file, "%d %d ", wires[w].bendPos[0],
              wires[w].bendPos[1]);
    }

    if (wires[w].bendCounts == 2)
    {
      fprintf(routes_file, "%d %d ", wires[w].bendPos[0],
              wires[w].bendPos[1]);
      fprintf(routes_file, "%d %d ", wires[w].bendPos[2],
              wires[w].bendPos[3]);
    }

    fprintf(routes_file, "%d %d\n", wires[w].Endpoints[2],
            wires[w].Endpoints[3]);
  }

  fclose(costs_file);
  fclose(routes_file);

  return 0;
}

void all_paths(wire_t *wire, cost_t *costs, int wireId)
{
  int x1 = (int)wire->Endpoints[0];
  int y1 = (int)wire->Endpoints[1];
  int x2 = (int)wire->Endpoints[2];
  int y2 = (int)wire->Endpoints[3];
  int bendPos[4];
  int best[4];
  int bestCount = wire->bendCounts;
  // copy(begin(wire->bendPos), end(wire->bendPos), begin(bendPos));
  memcpy(bendPos, wire->bendPos, sizeof(wire->bendPos));
  memcpy(best, wire->bendPos, sizeof(wire->bendPos));

  int minX = min(x1, x2);
  int minY = min(y1, y2);
  int maxX = max(x1, x2);
  int maxY = max(y1, y2);

  int width = maxX - minX;
  int height = maxY - minY;
  // printGrid(local_cost);

  int minCost = 1e7;
  int oldCost = minCost;

  int dirH = x2 > x1 ? 1 : -1;
  int dirV = y2 > y1 ? 1 : -1;

  // delete_wire(wire, grid, wireId);

  if (x1 != x2 && y1 != y2)
  { //  dx+dy paths
    //  vertical first | --
    //  1 bend vertical first
    bendPos[0] = x1;
    bendPos[1] = y2;
    // int cur = calculate_path_cost(wire, costs, bendPos, 1, wireId);
    int x1v[height + 1];
    int y2h[width + 1];
    int x2v[height + 1];
    int y1h[width + 1];

    // printPrefix(x1v, height+1);
    // printPrefix(x2v, height+1);
    // printPrefix(y1h, width+1);
    // printPrefix(y2h, width+1);

    calculateVerticalBorder(x1, y1, x1, y2, costs, wireId, x1v, dirV);
    calculateHorizontalBorder(x1, y2, x2, y2, costs, wireId, y2h, dirH);

    int cur = x1v[height] + y2h[width] - y2h[0]; //  v inclusive h exclusive
    if (cur < minCost)
    {
      minCost = cur;
      memcpy(best, bendPos, sizeof(bendPos));
      bestCount = 1;
      // printf("bestCount: %d, bends: %d %d cost: %d\n", bestCount, bendPos[0], bendPos[1], cur);
    }

    //  1 bend horizontal first
    bendPos[0] = x2;
    bendPos[1] = y1;
    calculateHorizontalBorder(x1, y1, x2, y1, costs, wireId, y1h, dirH);
    calculateVerticalBorder(x2, y1, x2, y2, costs, wireId, x2v, dirV);
    cur = x2v[height] + y1h[width] - x2v[0]; //  v inclusive h exclusive
    if (cur < minCost)
    {
      minCost = cur;
      // copy(begin(bendPos), end(bendPos), begin(best));
      memcpy(best, bendPos, sizeof(bendPos));
      bestCount = 1;
      // printf("bestCount: %d, bends: %d %d cost: %d\n", bestCount, bendPos[0], bendPos[1], cur);
    }
    //  vertical first  竖横竖
    bendPos[0] = x1;
    bendPos[2] = x2;
    int yCur = y1 + dirV;
    int index = 1;
    while (yCur != y2) //  2 bends
    {
      bendPos[1] = yCur;
      bendPos[3] = yCur;
      //  竖*2
      int cur = x1v[index] + x2v[height] - x2v[index - 1];
      //  横
      cur += calculateHorizontal(x1 + dirH, yCur, x2 - dirH, yCur, costs, wireId, dirH); //  both end inclusive
      index++;
      yCur += dirV;
      if (cur < minCost)
      {
        minCost = cur;
        // copy(begin(bendPos), end(bendPos), begin(best));
        memcpy(best, bendPos, sizeof(bendPos));
        bestCount = 2;
        // printf("bestCount: %d, bends: %d %d %d %dcost: %d\n", bestCount, bendPos[0], bendPos[1], bendPos[2], bendPos[3], cur);
      }
    }

    // printPrefix(x1v, height+1);
    // printPrefix(x2v, height+1);
    // printPrefix(y1h, width+1);
    // printPrefix(y2h, width+1);

    //  horizontal first 横竖横
    bendPos[1] = y1;
    bendPos[3] = y2;
    int xCur = x1 + dirH;
    index = 1;
    while (xCur != x2) //  2 bends
    {
      bendPos[0] = xCur;
      bendPos[2] = xCur;
      //  2*横
      int cur = y1h[index] + y2h[width] - y2h[index - 1];
      //  竖  -> cache miss， calculate together
      cur += calculateVertical(xCur, y1 + dirV, xCur, y2 - dirV, costs, wireId, dirV); //  both end inclusive
      index++;
      xCur += dirH;
      if (cur < minCost)
      {
        minCost = cur;
        // copy(begin(bendPos), end(bendPos), begin(best));
        memcpy(best, bendPos, sizeof(bendPos));
        bestCount = 2;
        // printf("bestCount: %d, bends: %d %d %d %dcost: %d\n", bestCount, bendPos[0], bendPos[1], bendPos[2], bendPos[3], cur);
      }
    }
    // copy(begin(best), end(best), begin(wire->bendPos));
  }
  // copy(begin(best), end(best), begin(wire->bendPos));
  // place_wire(wire, grid, wireId);
  if (minCost != 1e7)
  {
    memcpy(wire->bendPos, best, sizeof(best));
    wire->bendCounts = bestCount;
    // printf("!!!!!new cost %d best: %d %d %d %d\n", minCost, wire->bendPos[0], wire->bendPos[1], wire->bendPos[2], wire->bendPos[3]);
  }
  else
  {
    // printf("@@@@@old cost %d best: %d %d %d %d\n", minCost, wire->bendPos[0], wire->bendPos[1], wire->bendPos[2], wire->bendPos[3]);
  }
}

void printPrefix(int *A, int length)
{
  for (int i = 0; i < length; i++)
  {
    printf("%d ", A[i]);
  }
  printf("End prefix\n");
}

void simple_path(wire_t *wire)
{
  int x1 = wire->Endpoints[0];
  int y1 = wire->Endpoints[1];
  int x2 = wire->Endpoints[2];
  int y2 = wire->Endpoints[3];
  if (x1 == x2 || y1 == y2)
  {
    wire->bendCounts = 0;
    return;
  }

  // simply do vertical and then horizontal, only one bend
  wire->bendCounts = 1;
  wire->bendPos[0] = x1;
  wire->bendPos[1] = y2;
}

void random_path(wire_t *wire)
{
  /*
  50% start vertically, 50% start horizontally
  once start direction is set, multiply the probability to find the path
  */
  srand(time(NULL));
  int x1, y1, x2, y2, width, height;

  x1 = wire->Endpoints[0];
  y1 = wire->Endpoints[1];
  x2 = wire->Endpoints[2];
  y2 = wire->Endpoints[3];
  int minX = min(x1, x2);
  int minY = min(y1, y2);
  width = abs(x1 - x2);
  height = abs(y1 - y2);
  //  lie on a striaght line
  if (x1 == x2 || y1 == y2)
    return;

  // odd: start vertically
  // even : start horizontally
  if (rand() % 2 == 1)
  {
    //  vertically first | --
    int i = rand() % (height + 1);
    if (i == 0)
    {
      wire->bendPos[0] = x2;
      wire->bendPos[1] = y1;
      wire->bendCounts = 1;
    }
    else if (i == height)
    {
      wire->bendPos[0] = x1;
      wire->bendPos[1] = y2;
      wire->bendCounts = 1;
    }
    else
    {
      wire->bendPos[0] = x1;
      wire->bendPos[1] = minY + i;
      wire->bendPos[2] = x2;
      wire->bendPos[3] = minY + i;
      wire->bendCounts = 2; //  # bends
    }
  }
  else
  {
    //  horizontal first -- |
    int i = rand() % (width + 1);
    if (i == 0)
    {
      wire->bendPos[0] = x1;
      wire->bendPos[1] = y2;
      wire->bendCounts = 1;
    }
    else if (i == width)
    {
      wire->bendPos[0] = x2;
      wire->bendPos[1] = y1;
      wire->bendCounts = 1;
    }
    else
    {
      wire->bendPos[0] = minX + i;
      wire->bendPos[1] = y1;
      wire->bendPos[2] = minX + i;
      wire->bendPos[3] = y2;
      wire->bendCounts = 2; //  # bends
    }
  }
}

void calculateVerticalBorder(int xStart, int yStart, int xEnd, int yEnd, cost_t *costs, int wireId, int *xv, int dirV)
{
  int dim_x = costs->dim_x;
  int base = yStart * dim_x + xStart;
  int dirBase = dim_x * dirV;
  xv[0] = getCellCost(wireId, costs, base);
  base += dirBase;
  for (int idx = 1; idx <= abs(yStart - yEnd); idx++)
  {
    // printf("in x!: idx-1:%d\n", xv[idx-1]);
    xv[idx] = xv[idx - 1] + getCellCost(wireId, costs, base);
    base += dirBase;
  }
}

void calculateHorizontalBorder(int xStart, int yStart, int xEnd, int yEnd, cost_t *costs, int wireId, int *yh, int dirH)
{
  int dim_x = costs->dim_x;
  int base = yStart * dim_x + xStart;
  yh[0] = getCellCost(wireId, costs, base);
  base += dirH;
  for (int idx = 1; idx <= abs(xStart - xEnd); idx++)
  {
    // printf("int y!: idx-1:%d\n", yh[idx-1]);
    yh[idx] = yh[idx - 1] + getCellCost(wireId, costs, base);
    base += dirH;
  }
}

int calculateVertical(int xStart, int yStart, int xEnd, int yEnd, cost_t *costs, int wireId, int dirV)
{
  if ((yEnd - yStart) * dirV < 0)
    return 0;
  int res = 0;
  int dim_x = costs->dim_x;
  int base = yStart * dim_x + xStart;
  int dirBase = dim_x * dirV;
  for (int i = 0; i <= abs(yStart - yEnd); i++)
  {
    res += getCellCost(wireId, costs, base);
    base += dirBase;
  }
  return res;
}

int calculateHorizontal(int xStart, int yStart, int xEnd, int yEnd, cost_t *costs, int wireId, int dirH)
{
  if ((xEnd - xStart) * dirH < 0)
    return 0;
  int res = 0;
  int dim_x = costs->dim_x;
  int base = yStart * dim_x + xStart;
  for (int i = 0; i <= abs(xStart - xEnd); i++)
  {
    res += getCellCost(wireId, costs, base);
    base += dirH;
  }
  return res;
}

void printGrid(cost_t *costs)
{
  printf("Start printing grid\n");
  int width = costs->dim_x;
  int height = costs->dim_y;
  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      printf("%d ", costs->grid[i * width + j].cost);
    }
    printf("\n");
  }
}

//  place each wire after bends are set
void place_wire_line(int x1, int y1, int x2, int y2, cost_t *costs, int wireId)
{
  int dim_x = costs->dim_x;
  if (x1 == x2)
  {
    int minY = min(y1, y2);
    int maxY = max(y1, y2);
    int base = (minY + 1) * dim_x + x1;
    for (int i = minY + 1; i < maxY; i++)
    {
      incrCell(wireId, costs, base);
      base += dim_x;
    }
  }
  else if (y1 == y2)
  {
    int minX = min(x1, x2);
    int maxX = max(x1, x2);
    int base = y1 * dim_x + minX + 1;
    for (int i = minX + 1; i < maxX; i++)
    {
      incrCell(wireId, costs, base);
      base++;
    }
  }
}

void incrCell(int wireId, cost_t *costs, int offset)
{
  cell_t *cur = &costs->grid[offset];
  omp_set_lock(&cur->lock);
  cur->cost++;
  if (cur->wireCount < MAX_COUNT)
  {
    cur->record[cur->wireCount] = wireId;
    cur->wireCount++;
  }
  omp_unset_lock(&cur->lock);
}

//  place each wire after bends are set
void place_wire(wire_t *wire, cost_t *costs, int wireId)
{
  int xStart = wire->Endpoints[0];
  int yStart = wire->Endpoints[1];
  incrCell(wireId, costs, yStart * costs->dim_x + xStart);
  int count = wire->bendCounts;
  int xEnd;
  int yEnd;
  for (int i = 0; i < count; i++)
  {
    xEnd = wire->bendPos[i * 2];
    yEnd = wire->bendPos[i * 2 + 1];
    place_wire_line(xStart, yStart, xEnd, yEnd, costs, wireId);
    xStart = xEnd;
    yStart = yEnd;
    incrCell(wireId, costs, yStart * costs->dim_x + xStart);
  }
  xEnd = wire->Endpoints[2];
  yEnd = wire->Endpoints[3];
  place_wire_line(xStart, yStart, xEnd, yEnd, costs, wireId);
  incrCell(wireId, costs, yEnd * costs->dim_x + xEnd);
}

//  after all wires are placed, calculate cost for each wire segment
int collect_wire_line_cost(int x1, int y1, int x2, int y2, cost_t *costs, int wireId)
{
  int dim_x = costs->dim_x;
  int res = 0;
  if (x1 == x2)
  {
    int minY = min(y1, y2);
    int maxY = max(y1, y2);
    int base = (minY + 1) * dim_x + x1;
    for (int i = minY + 1; i < maxY; i++)
    {
      res += costs->grid[base].cost;
      base += dim_x;
    }
  }
  else if (y1 == y2)
  {
    int minX = min(x1, x2);
    int maxX = max(x1, x2);
    int base = y1 * dim_x + minX + 1;
    for (int i = minX + 1; i < maxX; i++)
    {
      res += costs->grid[base].cost;
      base++;
    }
  }
  return res;
}

//  after all wires are placed, calculate cost for each wire
void collect_wire_cost(wire_t *wire, cost_t *costs, int wireId)
{
  int xStart = wire->Endpoints[0];
  int yStart = wire->Endpoints[1];
  int count = wire->bendCounts;
  int xEnd;
  int yEnd;
  wire->cost += costs->grid[yStart * costs->dim_x + xStart].cost;
  for (int i = 0; i < count; i++)
  {
    xEnd = wire->bendPos[i * 2];
    yEnd = wire->bendPos[i * 2 + 1];
    wire->cost += collect_wire_line_cost(xStart, yStart, xEnd, yEnd, costs, wireId);
    xStart = xEnd;
    yStart = yEnd;
    wire->cost += costs->grid[yStart * costs->dim_x + xStart].cost;
  }
  xEnd = wire->Endpoints[2];
  yEnd = wire->Endpoints[3];
  wire->cost += collect_wire_line_cost(xStart, yStart, xEnd, yEnd, costs, wireId);
  wire->cost += costs->grid[yEnd * costs->dim_x + xEnd].cost;
}

//  after all wires are placed, assume put one line segment in grid and estimate cost
int calculate_path_line_cost(int x1, int y1, int x2, int y2, cost_t *costs, int wireId)
{
  int dim_x = costs->dim_x;
  int res = 0;
  if (x1 == x2)
  {
    int minY = min(y1, y2);
    int maxY = max(y1, y2);
    int base = (minY + 1) * dim_x + x1;
    for (int i = minY + 1; i < maxY; i++)
    {
      res += getCellCost(wireId, costs, base) + 1;
      base += dim_x;
    }
  }
  else if (y1 == y2)
  {
    int minX = min(x1, x2);
    int maxX = max(x1, x2);
    int base = y1 * dim_x + minX + 1;
    for (int i = minX + 1; i < maxX; i++)
    {
      res += getCellCost(wireId, costs, base) + 1;
      base++;
    }
  }
  return res;
}

// after all wires are placed, assume put one wire route in grid and estimate cost
int calculate_path_cost(wire_t *wire, cost_t *costs, int *bendPos, int count, int wireId)
{
  int xStart = wire->Endpoints[0];
  int yStart = wire->Endpoints[1];
  int res = 0;
  int xEnd;
  int yEnd;
  res += getCellCost(wireId, costs, yStart * costs->dim_x + xStart) + 1;
  for (int i = 0; i < count; i++)
  {
    xEnd = bendPos[i * 2];
    yEnd = bendPos[i * 2 + 1];
    res += calculate_path_line_cost(xStart, yStart, xEnd, yEnd, costs, wireId);
    xStart = xEnd;
    yStart = yEnd;
    res += getCellCost(wireId, costs, yStart * costs->dim_x + xStart) + 1;
  }
  xEnd = wire->Endpoints[2];
  yEnd = wire->Endpoints[3];
  res += calculate_path_line_cost(xStart, yStart, xEnd, yEnd, costs, wireId);
  res += getCellCost(wireId, costs, yEnd * costs->dim_x + xEnd) + 1;
  return res;
}

// int getCellCost(int wireId, cost_t *costs, int offset)
// {
//   cell_t *cur = &costs->grid[offset];
//   int res = 0;
//   for (int i = 0; i < cur->wireCount; i++)
//   {
//     if (wireId == cur->record[i])
//     {
//       res = cur->cost - 1;
//       break;
//     }
//   }
//   res = res == 0 ? cur->cost: res;
//   return res > 20 ? 200 : res;
// }

// int getCellCost(int wireId, cost_t *costs, int offset)
// {
//   cell_t *cur = &costs->grid[offset];
//   for (int i = 0; i < cur->wireCount; i++)
//   {
//     if (wireId == cur->record[i])
//     {
//       return cur->cost - 1;
//     }
//   }
//   return cur->cost;
// }

int getCellCost(int wireId, cost_t *costs, int offset)
{
  cell_t *cur = &costs->grid[offset];
  // for (int i = 0; i < cur->wireCount; i++)
  // {
  //   if (wireId == cur->record[i])
  //   {
  //     return cur->cost - 1;
  //   }
  // }
  return cur->cost;
}

void delete_wire_line(int x1, int y1, int x2, int y2, cost_t *costs, int wireId)
{
  if (x1 == x2)
  {
    int minY = min(y1, y2);
    int maxY = max(y1, y2);
    int base = (minY + 1) * costs->dim_x + x1;
    for (int i = minY + 1; i < maxY; i++)
    {
      decrCell(wireId, costs, base);
      base += costs->dim_x;
    }
  }
  else if (y1 == y2)
  {
    int minX = min(x1, x2);
    int maxX = max(x1, x2);
    int base = y1 * costs->dim_x + minX + 1;
    for (int i = minX + 1; i < maxX; i++)
    {
      decrCell(wireId, costs, base);
      base++;
    }
  }
}

void delete_wire(wire_t *wire, cost_t *costs, int wireId)
{
  int xStart = wire->Endpoints[0];
  int yStart = wire->Endpoints[1];
  decrCell(wireId, costs, yStart * costs->dim_x + xStart);
  int count = wire->bendCounts;
  int xEnd;
  int yEnd;
  for (int i = 0; i < count; i++)
  {
    xEnd = wire->bendPos[i * 2];
    yEnd = wire->bendPos[i * 2 + 1];
    delete_wire_line(xStart, yStart, xEnd, yEnd, costs, wireId);
    xStart = xEnd;
    yStart = yEnd;
    decrCell(wireId, costs, yStart * costs->dim_x + xStart);
  }
  xEnd = wire->Endpoints[2];
  yEnd = wire->Endpoints[3];
  delete_wire_line(xStart, yStart, xEnd, yEnd, costs, wireId);
  decrCell(wireId, costs, yEnd * costs->dim_x + xEnd);
}

void decrCell(int wireId, cost_t *costs, int offset)
{
  cell_t *cur = &costs->grid[offset];
  omp_set_lock(&cur->lock);
  cur->cost--;
  for (int i = 0; i < cur->wireCount; i++)
  {
    if (wireId == cur->record[i])
    {
      cur->record[i] = cur->record[cur->wireCount-1];
      cur->wireCount--;
      omp_unset_lock(&cur->lock);
      return;
    }
  }
  omp_unset_lock(&cur->lock);
}