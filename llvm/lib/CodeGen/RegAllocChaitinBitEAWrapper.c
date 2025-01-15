#include "RegAllocChaitinBitEAWrapper.h"
#include "RegAllocChaitinBitEA.h"
#include "RegAllocChaitinBitEAStdgraph.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct BitEA_Graph {
  unsigned nodeCount;
  block_t arr[];
};

struct BitEA_Weights {
  unsigned count;
  int arr[];
};

struct BitEA_Colors {
  unsigned nodeCount;
  unsigned colorCount;
  block_t arr[];
};

static size_t calcTotalBlockNum(size_t vertexNum) {
  return TOTAL_BLOCK_NUM(vertexNum);
}

struct BitEA_Graph *BitEA_Graph_create(unsigned nodeCount) {
  struct BitEA_Graph *graph =
      calloc(1, sizeof(struct BitEA_Graph) + sizeof(graph->arr[0]) *
                                                 (size_t)nodeCount *
                                                 calcTotalBlockNum(nodeCount));
  graph->nodeCount = nodeCount;
  return graph;
}

void BitEA_Graph_destroy(struct BitEA_Graph *graph) { free(graph); }

void BitEA_Graph_setEdge(struct BitEA_Graph *graph, unsigned node1,
                         unsigned node2) {
  size_t totalBlockNum = calcTotalBlockNum(graph->nodeCount);
  block_t(*edgesP)[][totalBlockNum] = (block_t(*)[][totalBlockNum])graph->arr;
  SET_EDGE(node1, node2, (*edgesP));
}

struct BitEA_Weights *BitEA_Weights_create(unsigned weightCount) {
  struct BitEA_Weights *weights = calloc(
      1, sizeof(struct BitEA_Weights) + sizeof(weights->arr[0]) * weightCount);
  weights->count = weightCount;
  return weights;
}

void BitEA_Weights_destroy(struct BitEA_Weights *weights) { free(weights); }

void BitEA_Weights_setWeight(struct BitEA_Weights *weights, unsigned index,
                             double weight, bool spillable) {
  weights->arr[index] =
      (spillable && isfinite(weight)) ? (weight * 100000.0) : INT_MAX;
}

struct BitEA_Colors *BitEA_Colors_create(unsigned nodeCount,
                                         unsigned colorCount) {
  struct BitEA_Colors *colors =
      calloc(1, sizeof(struct BitEA_Colors) + colorCount *
                                                  calcTotalBlockNum(nodeCount) *
                                                  sizeof(colors->arr[0]));
  colors->nodeCount = nodeCount;
  colors->colorCount = colorCount;
  return colors;
}

void BitEA_Colors_destroy(struct BitEA_Colors *colors) { free(colors); }

bool BitEA_Colors_getColor(const struct BitEA_Colors *colors, unsigned node,
                           unsigned *color) {
  size_t totalBlockNum = calcTotalBlockNum(colors->nodeCount);
  const block_t(*colorsP)[][totalBlockNum] =
      (const block_t(*)[][totalBlockNum])colors->arr;
  for (unsigned i = 0; i < colors->colorCount; ++i) {
    if (CHECK_COLOR((*colorsP)[i], node)) {
      *color = i;
      return true;
    }
  }
  return false;
}

static void finalize_bitea(block_t *graph, int graph_size, int *weights,
                           block_t *colors, int num_colors, block_t *pool) {
  block_t(*edge_mat)[][TOTAL_BLOCK_NUM(graph_size)] =
      (block_t(*)[][TOTAL_BLOCK_NUM(graph_size)])graph;

  memset(pool, 0, TOTAL_BLOCK_NUM(graph_size) * sizeof(block_t));

  block_t(*color_mat)[][TOTAL_BLOCK_NUM(graph_size)] =
      (block_t(*)[][TOTAL_BLOCK_NUM(graph_size)])colors;

  int conflict_count[graph_size];
  int total_conflicts, pool_total = 0;

  for (int i = 0; i < num_colors; i++) {
    total_conflicts = 0;
    count_conflicts(graph_size, (*color_mat)[i], graph, conflict_count);

    fix_conflicts(graph_size, graph, weights, conflict_count, &total_conflicts,
                  (*color_mat)[i], pool, &pool_total);
  }
}

struct BitEA_Colors *BitEA_solve(struct BitEA_Graph *graph,
                                 struct BitEA_Weights *weights,
                                 unsigned colorCount, unsigned populationSize,
                                 unsigned iterationCount) {
  int bestFitness;
  float bestSolutionTime;
  int uncoloredNum;
  struct BitEA_Colors *colors =
      BitEA_Colors_create(graph->nodeCount, colorCount);
  BitEA(graph->nodeCount, graph->arr, weights->arr, populationSize, colorCount,
        iterationCount, colors->arr, &bestFitness, &bestSolutionTime,
        &uncoloredNum);
  block_t *pool = malloc(sizeof(*pool) * calcTotalBlockNum(graph->nodeCount));
  finalize_bitea(graph->arr, graph->nodeCount, weights->arr, colors->arr,
                 colorCount, pool);
  free(pool);
  return colors;
}
