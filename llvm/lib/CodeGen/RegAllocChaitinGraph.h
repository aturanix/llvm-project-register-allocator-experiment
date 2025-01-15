#pragma once

#include "RegAllocChaitinRange.h"

#include <cassert>
#include <cstddef>
#include <optional>
#include <ostream>
#include <unordered_map>
#include <unordered_set>

namespace alihan {
class InterferenceGraph {
private:
  class Node {
  public:
    using EdgeIterator = std::unordered_set<unsigned>::const_iterator;

    Node() = delete;
    Node(double weight, bool spillable);

    [[nodiscard]] auto getWeight() const -> double;
    [[nodiscard]] auto getSpillable() const -> bool;
    [[nodiscard]] auto getEdgeCount() const -> std::size_t;

    [[nodiscard]] auto hasEdge(unsigned node) const -> bool;
    void addEdge(unsigned node);
    void removeEdge(unsigned node);

    [[nodiscard]] auto isLessThan(const Node &node) const -> bool;

    [[nodiscard]] auto begin() const -> EdgeIterator;
    [[nodiscard]] auto end() const -> EdgeIterator;

  private:
    double mWeight;
    bool mSpillable;
    std::unordered_set<unsigned> mEdges;
  };

public:
  using EdgeIterator = Node::EdgeIterator;

  class NodeIterator {
  private:
    using iterator = std::unordered_map<unsigned, Node>::const_iterator;

  public:
    using value_type = unsigned;

    NodeIterator(iterator node);
    auto operator++() -> NodeIterator &;
    [[nodiscard]] auto operator*() const -> unsigned;
    [[nodiscard]] auto operator==(NodeIterator const &it) const -> bool;
    [[nodiscard]] auto operator!=(NodeIterator const &it) const -> bool;

  private:
    iterator mNode;
  };

  [[nodiscard]] auto isEmpty() const -> bool;
  [[nodiscard]] auto getSize() const -> std::size_t;
  [[nodiscard]] auto getWeight(unsigned node) const -> std::optional<double>;
  [[nodiscard]] auto getSpillable(unsigned node) const -> std::optional<bool>;
  [[nodiscard]] auto getEdgeCount(unsigned node) const -> std::optional<std::size_t>;

  [[nodiscard]] auto hasNode(unsigned node) const -> bool;
  [[nodiscard]] auto hasEdge(unsigned node1, unsigned node2) const -> bool;
  void addNode(unsigned id, double weight, bool spillable);
  auto addEdge(unsigned node1, unsigned node2) -> bool;
  void removeNode(unsigned node);
  auto removeEdge(unsigned node1, unsigned node2) -> bool;

  [[nodiscard]] auto isNodeLessThan(unsigned node1, unsigned node2) const -> std::optional<bool>;

  [[nodiscard]] auto getNodeRange() const -> Range<NodeIterator>;
  [[nodiscard]] auto getEdgeRange(unsigned node) const -> std::optional<Range<EdgeIterator>>;

  auto print(std::ostream &os) const -> std::ostream &;

private:
  [[nodiscard]] auto getNode(unsigned node) const -> const Node *;
  [[nodiscard]] auto getNode(unsigned node) -> Node *;

  std::unordered_map<unsigned, Node> mGraph;
};
} // namespace alihan
