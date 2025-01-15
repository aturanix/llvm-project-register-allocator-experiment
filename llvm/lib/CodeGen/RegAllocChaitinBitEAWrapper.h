#ifndef LLVM_LIB_CODEGEN_REGALLOCCHAITINBITEAWRAPPER_H
#define LLVM_LIB_CODEGEN_REGALLOCCHAITINBITEAWRAPPER_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stdbool.h>

struct BitEA_Graph;
struct BitEA_Weights;
struct BitEA_Colors;

struct BitEA_Graph *BitEA_Graph_create(unsigned nodeCount);
void BitEA_Graph_destroy(struct BitEA_Graph *graph);
void BitEA_Graph_setEdge(struct BitEA_Graph *graph, unsigned node1, unsigned node2);

struct BitEA_Weights *BitEA_Weights_create(unsigned weightCount);
void BitEA_Weights_destroy(struct BitEA_Weights *weights);
void BitEA_Weights_setWeight(struct BitEA_Weights *weights, unsigned index, double weight, bool spillable);

struct BitEA_Colors *BitEA_Colors_create(unsigned nodeCount, unsigned colorCount);
void BitEA_Colors_destroy(struct BitEA_Colors *colors);
bool BitEA_Colors_getColor(const struct BitEA_Colors *colors, unsigned node, unsigned *color);

struct BitEA_Colors *BitEA_solve(struct BitEA_Graph *graph, struct BitEA_Weights *weights, unsigned colorCount, unsigned populationSize, unsigned iterationCount);

#ifdef __cplusplus
}
#endif
#endif
