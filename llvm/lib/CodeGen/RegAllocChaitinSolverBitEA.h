#pragma once

#include "RegAllocChaitinRegisters.h"

namespace alihan {
class InterferenceGraph;
[[nodiscard]] auto solveBitEA(const InterferenceGraph &graph,
                               unsigned numberOfColors) -> SolutionMap;
} // namespace alihan
