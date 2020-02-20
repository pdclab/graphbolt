// NOTE: This scheduler is taken from GraphChi's scheduler.
// Link:
// https://github.com/GraphChi/graphchi-cpp/blob/master/src/engine/bitset_scheduler.hpp
// Be careful before using it in production and ensure that it adheres to the
// requried license.

#ifndef __BITSET_SCHEDULER_HPP__
#define __BITSET_SCHEDULER_HPP__

#include "densebitset.h"

class BitsetScheduler {
private:
  DenseBitset *currBitset;
  DenseBitset *nextBitset;

public:
  bool scheduledTasks;

  BitsetScheduler(IdType nvertices) {
    currBitset = new DenseBitset(nvertices);
    nextBitset = new DenseBitset(nvertices);

    currBitset->clear();
    nextBitset->clear();

    scheduledTasks = false;
  }

  void newIteration() {
    DenseBitset *tmp = currBitset;
    currBitset = nextBitset;
    nextBitset = tmp;
    nextBitset->clear();
    scheduledTasks = false;
  }

  void reset() {
    currBitset->clear();
    nextBitset->clear();
  }

  virtual ~BitsetScheduler() {
    delete nextBitset;
    delete currBitset;
  }

  inline void schedule(IdType vertex, bool curr = true) {
    nextBitset->setBit(vertex);
    if (curr) {
      currBitset->setBit(vertex);
    }
    scheduledTasks = true;
  }

  inline void schedule(DenseBitset *other) {
    nextBitset->set(other);
    // currBitset->set(other);
  }

  void resize(IdType maxsize) {
    currBitset->resize(maxsize);
    nextBitset->resize(maxsize);
  }

  inline bool isScheduled(IdType vertex) { return currBitset->get(vertex); }

  void removeTasks(IdType fromvertex, IdType tovertex) {
    nextBitset->clearBits(fromvertex, tovertex);
  }

  void unschedule(IdType vertex, bool curr = true) {
    nextBitset->clearBit(vertex);
    if (curr) {
      currBitset->clearBit(vertex);
    }
  }

  void scheduleAll(bool curr = false) {
    nextBitset->setAll();
    if (curr) {
      // If possible, add to schedule already this iteration
      currBitset->setAll();
    }
    scheduledTasks = true;
  }

  IdType numTasks() { return currBitset->countSetBits(); }

  IdType numFutureTasks() { return nextBitset->countSetBits(); }

  bool anyScheduledTasks() { return scheduledTasks; }

  DenseBitset *getCurrentBitset() const { return currBitset; }

  DenseBitset *getNextBitset() const { return nextBitset; }
};

#endif
