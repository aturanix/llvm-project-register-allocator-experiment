#pragma once

#include "RegAllocChaitinGraph.h"

#include <optional>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace alihan {
using SolutionMap = std::unordered_map<unsigned, unsigned>;
using SolutionMapLLVM = std::unordered_map<unsigned, unsigned>;

class Registers {
public:
  struct VirtualRegister {
    double weight;
    bool spillable;
    std::unordered_set<unsigned> interferences;
    std::unordered_set<unsigned> candidatePhysRegs;
  };

  void addVirt(unsigned id, std::unordered_set<unsigned> candidatePhysIds,
               double weight, bool spillable);
  [[nodiscard]] auto getVirtCount() const -> unsigned;
  [[nodiscard]] auto getVirtReg(unsigned virtId) const -> const VirtualRegister *;
  [[nodiscard]] auto getVirtOrdinalId(unsigned virtId) const -> std::optional<unsigned>;
  [[nodiscard]] auto getVirtId(unsigned virtOrdinalId) const -> std::optional<unsigned>;
  [[nodiscard]] auto getGroupIdFirst() const -> unsigned;
  [[nodiscard]] auto getGroupIdLast() const -> unsigned;
  [[nodiscard]] auto getVirtOrdinalIdFirst() const -> unsigned;
  [[nodiscard]] auto getVirtOrdinalIdLast() const -> unsigned;
  auto addVirtInterference(unsigned virtId1, unsigned virtId2) -> bool;
  auto addPhys(unsigned id, std::vector<unsigned> const &subregIds) -> unsigned;
  [[nodiscard]] auto getGroupCount() const -> unsigned;
  [[nodiscard]] auto getPhysGroupId(unsigned physId) const -> std::optional<unsigned>;
  [[nodiscard]] auto getVirtCandPhysInGroup(unsigned virtId, unsigned groupId) const -> std::optional<unsigned>;
  [[nodiscard]] auto createInterferenceGraph() const -> InterferenceGraph;
  std::ostream &print(std::ostream &os) const;

private:
  [[nodiscard]] auto getVirtReg(unsigned virtId) -> VirtualRegister *;
  std::ostream &printVirtOrdinal(std::ostream &os) const;
  std::ostream &printVirt(std::ostream &os) const;
  std::ostream &printPhys(std::ostream &os) const;

  std::unordered_map<unsigned, VirtualRegister> mVirtRegs;
  std::unordered_map<unsigned, unsigned> mVirtToVirtOrdinal;
  std::vector<unsigned> mVirtOrdinalToVirt;
  std::unordered_map<unsigned, unsigned> mPhysToGroupidx;
  std::vector<std::unordered_set<unsigned>> mGroups;
};

enum class SolutionError {
  Success,
  HasEdgeNodesOfSameColor,
  HasUncoloredUnspillableNodes,
  HasNonexistentColor
};

[[nodiscard]] auto convertSolutionMapToSolutionMapLLVM(
    const InterferenceGraph &graph,
    const Registers &registers,
    const SolutionMap &solution,
    SolutionMapLLVM &solutionLLVM) -> SolutionError;

auto printSolutionMap(std::ostream &os, const SolutionMap &solution) -> std::ostream &;
} // namespace alihan
