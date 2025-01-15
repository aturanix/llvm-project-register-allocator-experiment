#include "RegAllocChaitinSolverBitEA.h"
#include "RegAllocChaitinBitEAWrapper.h"
#include "RegAllocChaitinGraph.h"

namespace alihan {
auto solveBitEA(const InterferenceGraph &graph,
                  unsigned numberOfColors) -> SolutionMap {
  BitEA_Graph *biteaGraph = BitEA_Graph_create(graph.getSize());
  for (unsigned node : graph.getNodeRange()) {
    auto edgeRange = graph.getEdgeRange(node).value();
    for (unsigned edgeNode : edgeRange) {
      BitEA_Graph_setEdge(biteaGraph, node, edgeNode);
    }
  }

  BitEA_Weights *biteaWeights = BitEA_Weights_create(graph.getSize());
  for (unsigned node : graph.getNodeRange()) {
    BitEA_Weights_setWeight(biteaWeights, node, graph.getWeight(node).value(), graph.getSpillable(node).value());
  }

  BitEA_Colors *biteaColors = BitEA_solve(biteaGraph, biteaWeights, numberOfColors, 100, 160000);

  SolutionMap solution;
  for (unsigned node : graph.getNodeRange()) {
    unsigned color;
    if (BitEA_Colors_getColor(biteaColors, node, &color)) {
      solution.insert({node, color});
    }
  }

  BitEA_Colors_destroy(biteaColors);
  BitEA_Weights_destroy(biteaWeights);
  BitEA_Graph_destroy(biteaGraph);
  return solution;
}
} // namespace alihan
