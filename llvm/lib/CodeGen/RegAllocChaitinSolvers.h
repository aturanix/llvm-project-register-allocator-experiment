#pragma once

#include "RegAllocChaitinRegisters.h"

namespace alihan {
class InterferenceGraph;
[[nodiscard]] auto solveGreedy(const InterferenceGraph &graph,
                               unsigned numberOfColors) -> SolutionMap;
[[nodiscard]] auto solveChaitin(const InterferenceGraph &graph,
                                unsigned numberOfColors) -> SolutionMap;
} // namespace alihan
