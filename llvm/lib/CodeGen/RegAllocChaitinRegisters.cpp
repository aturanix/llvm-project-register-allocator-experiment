#include "RegAllocChaitinRegisters.h"
#include "RegAllocChaitinGraph.h"

#include <algorithm>
#include <limits>
#include <optional>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace alihan {
void Registers::addVirt(unsigned id,
                        std::unordered_set<unsigned> candidatePhysIds,
                        double weight, bool spillable) {
  VirtualRegister reg;
  reg.weight = weight;
  reg.spillable = spillable;
  reg.candidatePhysRegs = std::move(candidatePhysIds);
  mVirtRegs.insert({id, std::move(reg)});
  mVirtToVirtOrdinal[id] = mVirtOrdinalToVirt.size();
  mVirtOrdinalToVirt.push_back(id);
}

auto Registers::getVirtCount() const -> unsigned { return mVirtRegs.size(); }

auto Registers::getVirtReg(unsigned virtId) const -> const VirtualRegister * {
  auto it = mVirtRegs.find(virtId);
  return (it == mVirtRegs.cend()) ? nullptr : &it->second;
}

auto Registers::getVirtOrdinalId(unsigned virtId) const
    -> std::optional<unsigned> {
  auto it = mVirtToVirtOrdinal.find(virtId);
  if (it == mVirtToVirtOrdinal.cend()) {
    return {};
  }
  return it->second + mGroups.size();
}

auto Registers::getVirtId(unsigned virtOrdinalId) const
    -> std::optional<unsigned> {
  unsigned i{virtOrdinalId - getVirtOrdinalIdFirst()};
  if (i >= mVirtOrdinalToVirt.size()) {
    return {};
  }
  return mVirtOrdinalToVirt[i];
}

auto Registers::getGroupIdFirst() const -> unsigned { return 0; }

auto Registers::getGroupIdLast() const -> unsigned { return getGroupCount(); }

auto Registers::getVirtOrdinalIdFirst() const -> unsigned {
  return getGroupIdLast();
}

auto Registers::getVirtOrdinalIdLast() const -> unsigned {
  return getVirtOrdinalIdFirst() + getVirtCount();
}

auto Registers::addVirtInterference(unsigned virtId1,
                                    unsigned virtId2) -> bool {
  auto itEnd = mVirtRegs.end();
  if (auto it1 = mVirtRegs.find(virtId1); it1 != itEnd) {
    if (auto it2 = mVirtRegs.find(virtId2); it2 != itEnd) {
      it1->second.interferences.insert(virtId2);
      it2->second.interferences.insert(virtId1);
      return true;
    }
  }
  return false;
}

auto Registers::addPhys(unsigned id,
                        const std::vector<unsigned> &subregIds) -> unsigned {
  auto endIt = mPhysToGroupidx.end();
  auto it = mPhysToGroupidx.find(id);
  if (it != endIt) {
    return it->second;
  }

  for (unsigned subregId : subregIds) {
    it = mPhysToGroupidx.find(subregId);
    if (it != endIt) {
      break;
    }
  }

  unsigned groupId;
  if (it == endIt) {
    groupId = mGroups.size();
    mGroups.emplace_back();
  } else {
    groupId = it->second;
  }

  std::unordered_set<unsigned> &group = mGroups[groupId];
  group.insert(id);
  mPhysToGroupidx.insert({id, groupId});
  for (unsigned subregId : subregIds) {
    group.insert(subregId);
    mPhysToGroupidx.insert({subregId, groupId});
  }
  return groupId;
}

auto Registers::getGroupCount() const -> unsigned { return mGroups.size(); }

auto Registers::getPhysGroupId(unsigned physId) const
    -> std::optional<unsigned> {
  auto it = mPhysToGroupidx.find(physId);
  if (it == mPhysToGroupidx.cend()) {
    return {};
  }
  return it->second;
}

auto Registers::getVirtCandPhysInGroup(unsigned virtId, unsigned groupId) const
    -> std::optional<unsigned> {
  VirtualRegister const *virtReg = getVirtReg(virtId);
  if (!virtReg) {
    return {};
  }

  for (unsigned phys_id : virtReg->candidatePhysRegs) {
    if (getPhysGroupId(phys_id) == groupId) {
      return phys_id;
    }
  }
  return {};
}

auto Registers::createInterferenceGraph() const -> InterferenceGraph {
  InterferenceGraph graph;
  for (unsigned group{getGroupIdFirst()}, e{getGroupIdLast()}; group != e;
       ++group) {
    graph.addNode(group, std::numeric_limits<double>::infinity(), false);
    for (unsigned j{0}; j != group; ++j) {
      graph.addEdge(group, j);
    }
  }

  for (unsigned virt{getVirtOrdinalIdFirst()}, e{getVirtOrdinalIdLast()};
       virt != e; ++virt) {
    const VirtualRegister *virtReg = getVirtReg(getVirtId(virt).value());
    graph.addNode(virt, virtReg->weight, virtReg->spillable);

    std::unordered_set<unsigned> candidateGroups;
    for (unsigned candPhysId : virtReg->candidatePhysRegs) {
      candidateGroups.insert(getPhysGroupId(candPhysId).value());
    }

    for (unsigned group{getGroupIdFirst()}, e{getGroupIdLast()}; group != e;
         ++group) {
      if (!candidateGroups.count(group)) {
        graph.addEdge(virt, group);
      }
    }
  }

  for (unsigned virt{getVirtOrdinalIdFirst()}, e{getVirtOrdinalIdLast()};
       virt != e; ++virt) {
    const VirtualRegister *virtReg = getVirtReg(getVirtId(virt).value());
    for (unsigned interference : virtReg->interferences) {
      graph.addEdge(virt, getVirtOrdinalId(interference).value());
    }
  }
  return graph;
}

auto Registers::print(std::ostream &os) const -> std::ostream & {
  printVirt(os) << '\n';
  return printPhys(os);
}

auto Registers::getVirtReg(unsigned int virtId) -> VirtualRegister * {
  auto it = mVirtRegs.find(virtId);
  return (it == mVirtRegs.end()) ? nullptr : &it->second;
}

auto Registers::printVirtOrdinal(std::ostream &os) const -> std::ostream & {
  os << "ORTOVI: {";
  bool first_ordinal{true};
  for (unsigned i{0}; i != mVirtOrdinalToVirt.size(); ++i) {
    if (!first_ordinal) {
      os << ", ";
    }
    unsigned ordinal = i + mGroups.size();
    os << ordinal << " -> " << *getVirtId(ordinal);
    first_ordinal = false;
  }
  return os << '}';
}

auto Registers::printVirt(std::ostream &os) const -> std::ostream & {
  bool firstVirt{true};
  os << '{';
  for (auto [virtId, virtReg] : mVirtRegs) {
    if (!firstVirt) {
      os << ",\n";
    }
    os << virtId << ": {" << "WE: " << virtReg.weight
       << ", IS: " << (virtReg.spillable ? "true" : "false") << ", PH: {";
    bool firstPhys{true};
    for (unsigned physId : virtReg.candidatePhysRegs) {
      if (!firstPhys) {
        os << ", ";
      }
      os << physId;
      firstPhys = false;
    }
    os << "}, IN: {";
    bool firstInterference{true};
    for (unsigned interferenceVirtId : virtReg.interferences) {
      if (!firstInterference) {
        os << ", ";
      }
      os << interferenceVirtId;
      firstInterference = false;
    }
    os << "}}";
    firstVirt = false;
  }
  os << "}\n";
  return printVirtOrdinal(os);
}

auto Registers::printPhys(std::ostream &os) const -> std::ostream & {
  os << "PHGRPS [";
  bool firstGroup{true};
  for (auto const &group : mGroups) {
    if (!firstGroup) {
      os << ", ";
    }
    os << '{';
    bool firstPhys{true};
    for (unsigned phys : group) {
      if (!firstPhys) {
        os << ", ";
      }
      os << phys;
      firstPhys = false;
    }
    os << '}';
    firstGroup = false;
  }
  os << ']';
  return os;
}

auto checkSolutionMapEdges(const InterferenceGraph &graph, const SolutionMap &solution) -> bool {
  auto solutionEnd = solution.end();
  for (unsigned n1 : graph.getNodeRange()) {
    auto it1 = solution.find(n1);
    if (it1 == solutionEnd) {
      continue;
    }

    auto edgeRange = graph.getEdgeRange(n1).value();
    for (unsigned n2 : edgeRange) {
      auto it2 = solution.find(n2);
      if (it2 != solutionEnd && it1->second == it2->second) {
        return false;
      }
    }  
  }
  return true;
}

auto checkSolutionMapUnspillable(const InterferenceGraph &graph, const SolutionMap &solution) -> bool {
  for (unsigned n : graph.getNodeRange()) {
    if (!graph.getSpillable(n).value() && !solution.count(n)) {
      return false;
    }
  }
  return true;
}

auto convertSolutionMapToSolutionMapLLVM(const InterferenceGraph &graph,
                                         const Registers &registers,
                                         const SolutionMap &solution,
                                         SolutionMapLLVM &solutionLLVM) -> SolutionError {
  if (!checkSolutionMapEdges(graph, solution)) {
    return SolutionError::HasEdgeNodesOfSameColor;
  }

  if (!checkSolutionMapUnspillable(graph, solution)) {
    return SolutionError::HasUncoloredUnspillableNodes;
  }

  std::unordered_map<std::size_t, unsigned> colorToGroupMapping;
  for (unsigned group{registers.getGroupIdFirst()},
       e{registers.getGroupIdLast()};
       group != e; ++group) {
    auto it = solution.find(group);
    colorToGroupMapping.insert({it->second, group});
  }

  unsigned numberOfColors{registers.getGroupCount()};

  SolutionMapLLVM tSolutionLLVM;
  unsigned virtOrdinalIdFirst{registers.getVirtOrdinalIdFirst()};
  unsigned virtOrdinalIdLast{registers.getVirtOrdinalIdLast()};
  for (auto [reg, color] : solution) {
    if (color >= numberOfColors) {
      return SolutionError::HasNonexistentColor;
    }

    if (virtOrdinalIdFirst <= reg && reg < virtOrdinalIdLast) {
      unsigned virtId{registers.getVirtId(reg).value()};
      unsigned physId{registers.getVirtCandPhysInGroup(virtId, colorToGroupMapping[color]).value()};
      tSolutionLLVM.insert({virtId, physId});
    }
  }
  solutionLLVM = std::move(tSolutionLLVM);
  return SolutionError::Success;
}

auto printSolutionMap(std::ostream &os, const SolutionMap &solution) -> std::ostream & {
  std::vector<std::pair<unsigned, unsigned>> sortedNodes(solution.begin(), solution.end());
  auto compFunc = [](const std::pair<unsigned, unsigned> &n1, const std::pair<unsigned, unsigned> &n2) {
    return n1.first < n2.first;
  };
  std::sort(sortedNodes.begin(), sortedNodes.end(), compFunc);

  bool first{true};
  for (auto [node, color] : sortedNodes) {
    if (!first) {
      os << '\n';
    }
    os << node << ' ' << color;
    first = false;
  }
  return os;
}
} // namespace alihan
