#pragma once

namespace alihan {
template <typename Iterator> class Range {
public:
  Range(Iterator beginIt, Iterator endIt) : mBeginIt(beginIt), mEndIt(endIt) {}
  auto begin() -> Iterator { return mBeginIt; }
  auto end() -> Iterator { return mEndIt; }

private:
  Iterator mBeginIt;
  Iterator mEndIt;
};
} // namespace alihan