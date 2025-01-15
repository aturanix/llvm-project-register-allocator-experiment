//===-- RegAllocChaitin.cpp - Chaitin Register Allocator ----------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

#include "RegAllocChaitinRegisters.h"
#include "RegAllocChaitinSolverBitEA.h"
#include "RegAllocChaitinSolvers.h"
#include "RegAllocChaitinGraph.h"

#include "AllocationOrder.h"
#include "LiveDebugVariables.h"
#include "RegAllocBase.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/CodeGen/CalcSpillWeights.h"
#include "llvm/CodeGen/LiveIntervals.h"
#include "llvm/CodeGen/LiveRangeEdit.h"
#include "llvm/CodeGen/LiveRegMatrix.h"
#include "llvm/CodeGen/LiveStacks.h"
#include "llvm/CodeGen/MachineBlockFrequencyInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineLoopInfo.h"
#include "llvm/CodeGen/Passes.h"
#include "llvm/CodeGen/RegAllocRegistry.h"
#include "llvm/CodeGen/Spiller.h"
#include "llvm/CodeGen/TargetRegisterInfo.h"
#include "llvm/CodeGen/VirtRegMap.h"
#include "llvm/Pass.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include <queue>
#include <sstream>

using namespace llvm;

#define DEBUG_TYPE "regalloc"

static RegisterRegAlloc chaitinRegAlloc("chaitin", "chaitin register allocator",
                                      createChaitinRegisterAllocator);

namespace {
  struct CompSpillWeight {
    bool operator()(const LiveInterval *A, const LiveInterval *B) const {
      return A->weight() < B->weight();
    }
  };
}

namespace {
class RAChaitin : public MachineFunctionPass,
                public RegAllocBase,
                private LiveRangeEdit::Delegate {
  // context
  MachineFunction *MF = nullptr;

  // state
  std::unique_ptr<Spiller> SpillerInstance;
  std::priority_queue<const LiveInterval *, std::vector<const LiveInterval *>,
                      CompSpillWeight>
      Queue;

  // Scratch space.  Allocated here to avoid repeated malloc calls in
  // selectOrSplit().
  BitVector UsableRegs;

  bool LRE_CanEraseVirtReg(Register) override;
  void LRE_WillShrinkVirtReg(Register) override;

public:
  RAChaitin(const RegClassFilterFunc F = allocateAllRegClasses);

  /// Return the pass name.
  StringRef getPassName() const override { return "Chaitin Register Allocator"; }

  /// RAChaitin analysis usage.
  void getAnalysisUsage(AnalysisUsage &AU) const override;

  void releaseMemory() override;

  Spiller &spiller() override { return *SpillerInstance; }

  void enqueueImpl(const LiveInterval *LI) override { Queue.push(LI); }

  const LiveInterval *dequeue() override {
    if (Queue.empty())
      return nullptr;
    const LiveInterval *LI = Queue.top();
    Queue.pop();
    return LI;
  }

  MCRegister selectOrSplit(const LiveInterval &VirtReg,
                           SmallVectorImpl<Register> &SplitVRegs) override;

  /// Perform register allocation.
  bool runOnMachineFunction(MachineFunction &mf) override;

  MachineFunctionProperties getRequiredProperties() const override {
    return MachineFunctionProperties().set(
        MachineFunctionProperties::Property::NoPHIs);
  }

  MachineFunctionProperties getClearedProperties() const override {
    return MachineFunctionProperties().set(
      MachineFunctionProperties::Property::IsSSA);
  }

  // Helper for spilling all live virtual registers currently unified under preg
  // that interfere with the most recently queried lvr.  Return true if spilling
  // was successful, and append any new spilled/split intervals to splitLVRs.
  bool spillInterferences(const LiveInterval &VirtReg, MCRegister PhysReg,
                          SmallVectorImpl<Register> &SplitVRegs);

  static char ID;

private:
  unsigned assignRemainingIntervals(
      std::function<alihan::SolutionMap(const alihan::InterferenceGraph &, unsigned)>
          Solver);
};

char RAChaitin::ID = 0;

} // end anonymous namespace

char &llvm::RAChaitinID = RAChaitin::ID;

INITIALIZE_PASS_BEGIN(RAChaitin, "regallocchaitin", "Chaitin Register Allocator",
                      false, false)
INITIALIZE_PASS_DEPENDENCY(LiveDebugVariables)
INITIALIZE_PASS_DEPENDENCY(SlotIndexes)
INITIALIZE_PASS_DEPENDENCY(LiveIntervals)
INITIALIZE_PASS_DEPENDENCY(RegisterCoalescer)
INITIALIZE_PASS_DEPENDENCY(MachineScheduler)
INITIALIZE_PASS_DEPENDENCY(LiveStacks)
INITIALIZE_PASS_DEPENDENCY(AAResultsWrapperPass)
INITIALIZE_PASS_DEPENDENCY(MachineDominatorTree)
INITIALIZE_PASS_DEPENDENCY(MachineLoopInfo)
INITIALIZE_PASS_DEPENDENCY(VirtRegMap)
INITIALIZE_PASS_DEPENDENCY(LiveRegMatrix)
INITIALIZE_PASS_END(RAChaitin, "regallocchaitin", "Chaitin Register Allocator", false,
                    false)

bool RAChaitin::LRE_CanEraseVirtReg(Register VirtReg) {
  LiveInterval &LI = LIS->getInterval(VirtReg);
  if (VRM->hasPhys(VirtReg)) {
    Matrix->unassign(LI);
    aboutToRemoveInterval(LI);
    return true;
  }
  // Unassigned virtreg is probably in the priority queue.
  // RegAllocBase will erase it after dequeueing.
  // Nonetheless, clear the live-range so that the debug
  // dump will show the right state for that VirtReg.
  LI.clear();
  return false;
}

void RAChaitin::LRE_WillShrinkVirtReg(Register VirtReg) {
  if (!VRM->hasPhys(VirtReg))
    return;

  // Register is assigned, put it back on the queue for reassignment.
  LiveInterval &LI = LIS->getInterval(VirtReg);
  Matrix->unassign(LI);
  enqueue(&LI);
}

RAChaitin::RAChaitin(RegClassFilterFunc F):
  MachineFunctionPass(ID),
  RegAllocBase(F) {
}

void RAChaitin::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesCFG();
  AU.addRequired<AAResultsWrapperPass>();
  AU.addPreserved<AAResultsWrapperPass>();
  AU.addRequired<LiveIntervals>();
  AU.addPreserved<LiveIntervals>();
  AU.addPreserved<SlotIndexes>();
  AU.addRequired<LiveDebugVariables>();
  AU.addPreserved<LiveDebugVariables>();
  AU.addRequired<LiveStacks>();
  AU.addPreserved<LiveStacks>();
  AU.addRequired<MachineBlockFrequencyInfo>();
  AU.addPreserved<MachineBlockFrequencyInfo>();
  AU.addRequiredID(MachineDominatorsID);
  AU.addPreservedID(MachineDominatorsID);
  AU.addRequired<MachineLoopInfo>();
  AU.addPreserved<MachineLoopInfo>();
  AU.addRequired<VirtRegMap>();
  AU.addPreserved<VirtRegMap>();
  AU.addRequired<LiveRegMatrix>();
  AU.addPreserved<LiveRegMatrix>();
  MachineFunctionPass::getAnalysisUsage(AU);
}

void RAChaitin::releaseMemory() {
  SpillerInstance.reset();
}


// Spill or split all live virtual registers currently unified under PhysReg
// that interfere with VirtReg. The newly spilled or split live intervals are
// returned by appending them to SplitVRegs.
bool RAChaitin::spillInterferences(const LiveInterval &VirtReg,
                                 MCRegister PhysReg,
                                 SmallVectorImpl<Register> &SplitVRegs) {
  // Record each interference and determine if all are spillable before mutating
  // either the union or live intervals.
  SmallVector<const LiveInterval *, 8> Intfs;

  // Collect interferences assigned to any alias of the physical register.
  for (MCRegUnit Unit : TRI->regunits(PhysReg)) {
    LiveIntervalUnion::Query &Q = Matrix->query(VirtReg, Unit);
    for (const auto *Intf : reverse(Q.interferingVRegs())) {
      if (!Intf->isSpillable() || Intf->weight() > VirtReg.weight())
        return false;
      Intfs.push_back(Intf);
    }
  }
  LLVM_DEBUG(dbgs() << "spilling " << printReg(PhysReg, TRI)
                    << " interferences with " << VirtReg << "\n");
  assert(!Intfs.empty() && "expected interference");

  // Spill each interfering vreg allocated to PhysReg or an alias.
  for (unsigned i = 0, e = Intfs.size(); i != e; ++i) {
    const LiveInterval &Spill = *Intfs[i];

    // Skip duplicates.
    if (!VRM->hasPhys(Spill.reg()))
      continue;

    // Deallocate the interfering vreg by removing it from the union.
    // A LiveInterval instance may not be in a union during modification!
    Matrix->unassign(Spill);

    // Spill the extracted interval.
    LiveRangeEdit LRE(&Spill, SplitVRegs, *MF, *LIS, VRM, this, &DeadRemats);
    spiller().spill(LRE);
  }
  return true;
}

// Driver for the register assignment and splitting heuristics.
// Manages iteration over the LiveIntervalUnions.
//
// This is a minimal implementation of register assignment and splitting that
// spills whenever we run out of registers.
//
// selectOrSplit can only be called once per live virtual register. We then do a
// single interference test for each register the correct class until we find an
// available register. So, the number of interference tests in the worst case is
// |vregs| * |machineregs|. And since the number of interference tests is
// minimal, there is no value in caching them outside the scope of
// selectOrSplit().
MCRegister RAChaitin::selectOrSplit(const LiveInterval &VirtReg,
                                  SmallVectorImpl<Register> &SplitVRegs) {
  // Populate a list of physical register spill candidates.
  SmallVector<MCRegister, 8> PhysRegSpillCands;

  // Check for an available register in this class.
  auto Order =
      AllocationOrder::create(VirtReg.reg(), *VRM, RegClassInfo, Matrix);
  for (MCRegister PhysReg : Order) {
    assert(PhysReg.isValid());
    // Check for interference in PhysReg
    switch (Matrix->checkInterference(VirtReg, PhysReg)) {
    case LiveRegMatrix::IK_Free:
      // PhysReg is available, allocate it.
      return PhysReg;

    case LiveRegMatrix::IK_VirtReg:
      // Only virtual registers in the way, we may be able to spill them.
      PhysRegSpillCands.push_back(PhysReg);
      continue;

    default:
      // RegMask or RegUnit interference.
      continue;
    }
  }

  // Try to spill another interfering reg with less spill weight.
  for (MCRegister &PhysReg : PhysRegSpillCands) {
    if (!spillInterferences(VirtReg, PhysReg, SplitVRegs))
      continue;

    assert(!Matrix->checkInterference(VirtReg, PhysReg) &&
           "Interference after spill.");
    // Tell the caller to allocate to this newly freed physical register.
    return PhysReg;
  }

  // No other spill candidates were found, so spill the current VirtReg.
  LLVM_DEBUG(dbgs() << "spilling: " << VirtReg << '\n');
  if (!VirtReg.isSpillable())
    return ~0u;
  LiveRangeEdit LRE(&VirtReg, SplitVRegs, *MF, *LIS, VRM, this, &DeadRemats);
  spiller().spill(LRE);

  // The live virtual register requesting allocation was spilled, so tell
  // the caller not to allocate anything during this round.
  return 0;
}


unsigned RAChaitin::assignRemainingIntervals(
    std::function<alihan::SolutionMap(const alihan::InterferenceGraph &, unsigned)>
        Solver) {
  std::vector<const LiveInterval *> Intervals;
  for (unsigned I{0u}, E = MRI->getNumVirtRegs(); I != E; ++I) {
    Register Reg = Register::index2VirtReg(I);
    if (MRI->reg_nodbg_empty(Reg) || VRM->hasPhys(Reg)) {
      continue;
    }
    Intervals.push_back(&LIS->getInterval(Reg));
  }

  if (Intervals.empty()) {
    return 0;
  }

  alihan::Registers RegsData;
  for (const LiveInterval *VirtReg : Intervals) {
    auto Order =
        AllocationOrder::create(VirtReg->reg(), *VRM, RegClassInfo, Matrix);
    LLVM_DEBUG(dbgs() << *VirtReg << '\n');

    std::unordered_set<unsigned> CandidatePhys;
    for (MCRegister PhysReg : Order) {
      assert(PhysReg.isValid());
      if (Matrix->checkInterference(*VirtReg, PhysReg) ==
          LiveRegMatrix::IK_Free) {
        CandidatePhys.insert(PhysReg);
      }

      auto SubregsRange = TRI->subregs(PhysReg);
      std::vector<unsigned> Subregs(SubregsRange.begin(), SubregsRange.end());
      RegsData.addPhys(PhysReg, Subregs);
    }
    RegsData.addVirt(VirtReg->reg(), std::move(CandidatePhys),
                     VirtReg->weight(), VirtReg->isSpillable());
  }

  for (unsigned I{0u}; I != Intervals.size(); ++I) {
    for (unsigned J{I + 1}; J != Intervals.size(); ++J) {
      bool Interference{!Intervals[I]->empty() &&
                        Intervals[I]->overlaps(*Intervals[J])};
      if (Interference) {
        RegsData.addVirtInterference(Intervals[I]->reg(), Intervals[J]->reg());
      }
    }
  }

  alihan::InterferenceGraph Graph = RegsData.createInterferenceGraph();

  std::ostringstream Oss;
  Oss << "GRAPH\n";
  Graph.print(Oss) << '\n';

  alihan::SolutionMap Solution = Solver(Graph, RegsData.getGroupCount());

  Oss << "NUMBEROFCOLORS " << RegsData.getGroupCount() << '\n';
  Oss << "SOLUTION\n";
  alihan::printSolutionMap(Oss, Solution) << '\n';
  LLVM_DEBUG(dbgs() << Oss.str());

  alihan::SolutionMapLLVM SolutionLLVM;
  alihan::SolutionError Error = alihan::convertSolutionMapToSolutionMapLLVM(Graph, RegsData, Solution, SolutionLLVM);
  switch (Error) {
  case alihan::SolutionError::Success:
    break;
  case alihan::SolutionError::HasEdgeNodesOfSameColor:
    LLVM_DEBUG(dbgs() << "Generated solution has edge nodes of same color\n");
    return 0;
  case alihan::SolutionError::HasNonexistentColor:
    LLVM_DEBUG(dbgs() << "Generated solution has a nonexistent color\n");
    return 0;
  case alihan::SolutionError::HasUncoloredUnspillableNodes:
    LLVM_DEBUG(dbgs() << "Generated solution has uncolored unspillable nodes\n");
    return 0;
  }

  LLVM_DEBUG(dbgs() << "Generated solution has " << SolutionLLVM.size() << " assignments\n");

  for (auto [VirtId, PhysId] : SolutionLLVM) {
    Matrix->assign(LIS->getInterval(VirtId), PhysId);
  }

  Matrix->invalidateVirtRegs();
  return SolutionLLVM.size();
}


bool RAChaitin::runOnMachineFunction(MachineFunction &mf) {
  LLVM_DEBUG(dbgs() << "********** CHAITIN REGISTER ALLOCATION **********\n"
                    << "********** Function: " << mf.getName() << '\n');

  MF = &mf;
  RegAllocBase::init(getAnalysis<VirtRegMap>(),
                     getAnalysis<LiveIntervals>(),
                     getAnalysis<LiveRegMatrix>());
  VirtRegAuxInfo VRAI(*MF, *LIS, *VRM, getAnalysis<MachineLoopInfo>(),
                      getAnalysis<MachineBlockFrequencyInfo>());
  VRAI.calculateSpillWeightsAndHints();

  SpillerInstance.reset(createInlineSpiller(*this, *MF, *VRM, VRAI));

  unsigned N = assignRemainingIntervals(alihan::solveBitEA);
  LLVM_DEBUG(dbgs() << "Assigned " << N << " intervals\n");

  allocatePhysRegs();
  postOptimization();

  // Diagnostic output before rewriting
  LLVM_DEBUG(dbgs() << "Post alloc VirtRegMap:\n" << *VRM << "\n");

  releaseMemory();
  return true;
}

FunctionPass* llvm::createChaitinRegisterAllocator() {
  return new RAChaitin();
}

FunctionPass* llvm::createChaitinRegisterAllocator(RegClassFilterFunc F) {
  return new RAChaitin(F);
}
