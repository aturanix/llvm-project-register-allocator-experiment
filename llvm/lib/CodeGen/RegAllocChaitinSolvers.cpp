#include "RegAllocChaitinSolvers.h"
#include "RegAllocChaitinGraph.h"
#include "RegAllocChaitinRegisters.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <queue>
#include <vector>

namespace {
auto findUnusedColor(const alihan::InterferenceGraph &graph,
                     unsigned numberOfColors,
                     const alihan::SolutionMap &solution,
                     unsigned node) -> std::optional<unsigned> {
  std::vector<char> colorUsage(numberOfColors);
  if (auto edgeRangeOpt = graph.getEdgeRange(node)) {
    for (unsigned edge : *edgeRangeOpt) {
      auto it = solution.find(edge);
      if (it != solution.end()) {
        colorUsage[it->second] = true;
      }
    }
  } else {
    return {};
  }

  for (unsigned color{0}; color != numberOfColors; ++color) {
    if (!colorUsage[color]) {
      return color;
    }
  }
  return {};
}
} // namespace

namespace alihan {
auto solveGreedy(const InterferenceGraph &graph,
                 unsigned numberOfColors) -> SolutionMap {
  auto comp = [&](unsigned const &virt1, unsigned const &virt2) {
    return graph.isNodeLessThan(virt1, virt2).value();
  };

  std::priority_queue<unsigned, std::vector<unsigned>, decltype(comp)> virts(
      comp);

  for (unsigned node : graph.getNodeRange()) {
    virts.push(node);
  }

  SolutionMap solution;

  while (!virts.empty()) {
    unsigned virt = virts.top();
    virts.pop();

    std::optional<unsigned> color =
        findUnusedColor(graph, numberOfColors, solution, virt);
    if (color.has_value()) {
      solution.insert({virt, color.value()});
    }
  }

  return solution;
}

auto solveChaitin(const InterferenceGraph &graph,
                  unsigned numberOfColors) -> SolutionMap {
  InterferenceGraph tempGraph = graph;

  auto isLess = [&](unsigned node1, unsigned node2) {
    double w1{tempGraph.getWeight(node1).value()};
    double w2{tempGraph.getWeight(node2).value()};
    std::size_t e1{tempGraph.getEdgeCount(node1).value()};
    std::size_t e2{tempGraph.getEdgeCount(node2).value()};
    bool s1{tempGraph.getSpillable(node1).value()};
    bool s2{tempGraph.getSpillable(node2).value()};
    if (s1 && s2) {
      return w1 / e1 < w2 / e2;
    } else if (!s1 && !s2) {
      return e1 > e2;
    } else {
      return s1;
    }
  };

  std::vector<unsigned> stack;
  while (!tempGraph.isEmpty()) {
    std::optional<unsigned> max;
    for (unsigned node : tempGraph.getNodeRange()) {
      if (tempGraph.getEdgeCount(node).value() >= numberOfColors) {
        continue;
      }

      if (!max || (max && isLess(*max, node))) {
        max = node;
      }
    }

    if (max) {
      stack.push_back(*max);
      tempGraph.removeNode(*max);
    } else {
      auto nodeRange = tempGraph.getNodeRange();
      unsigned min =
          *std::min_element(nodeRange.begin(), nodeRange.end(), isLess);
      tempGraph.removeNode(min);
    }
  }

  SolutionMap solution;

  while (!stack.empty()) {
    unsigned node = stack.back();
    stack.pop_back();
    tempGraph.addNode(node, 0.0, true);
    auto range = graph.getEdgeRange(node);
    for (unsigned edge : *range) {
      if (tempGraph.hasNode(edge)) {
        tempGraph.addEdge(node, edge);
      }
    }
    unsigned color =
        findUnusedColor(tempGraph, numberOfColors, solution, node).value();
    solution.insert({node, color});
  }
  return solution;
}
} // namespace alihan
