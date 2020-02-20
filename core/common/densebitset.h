// NOTE: This bitset is taken from GraphChi's scheduler.
// Link:
// https://github.com/GraphChi/graphchi-cpp/blob/master/src/util/DenseBitset.hpp
// Be careful before using it in production and ensure that it adheres to the
// requried license.

#ifndef __DENSE_BITSET_HPP__
#define __DENSE_BITSET_HPP__
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

typedef uint32_t IdType;

class DenseBitset {
public:
  DenseBitset() : array(NULL), len(0), arrlen(0) {}

  DenseBitset(IdType size) : array(NULL), len(size) {
    resize(size);
    clear();
  }

  virtual ~DenseBitset() { free(array); }

  void resize(IdType n) {
    len = n;
    arrlen = n / (8 * sizeof(IdType)) + 1;
    array = (IdType *)realloc(array, sizeof(IdType) * arrlen);
  }

  void clear() { memset(array, 0, arrlen * sizeof(IdType)); }

  void setAll() { memset(array, 0xff, arrlen * sizeof(IdType)); }

  inline bool get(IdType b) const {
    IdType arrpos, bitpos;
    bitToPos(b, arrpos, bitpos);
    return array[arrpos] & (IdType(1) << IdType(bitpos));
  }

  // Set the bit returning the old value
  inline bool setBit(IdType b) {
    // use CAS to set the bit
    IdType arrpos, bitpos;
    bitToPos(b, arrpos, bitpos);
    const IdType mask(IdType(1) << IdType(bitpos));
    return __sync_fetch_and_or(array + arrpos, mask) & mask;
  }

  // Set the state of the bit returning the old value
  inline bool set(IdType b, bool value) {
    if (value)
      return setBit(b);
    else
      return clearBit(b);
  }

  inline void set(DenseBitset *other) {
    assert(len == other->size());
    memcpy(array, other->getArray(), arrlen * sizeof(IdType));
  }

  // Clear the bit returning the old value
  inline bool clearBit(IdType b) {
    // use CAS to set the bit
    IdType arrpos, bitpos;
    bitToPos(b, arrpos, bitpos);
    const IdType test_mask(IdType(1) << IdType(bitpos));
    const IdType clear_mask(~test_mask);
    return __sync_fetch_and_and(array + arrpos, clear_mask) & test_mask;
  }

  inline void clearBits(IdType fromb, IdType tob) { // tob is inclusive
    // Careful with alignment
    const IdType bitsperword = sizeof(IdType) * 8;
    while ((fromb % bitsperword != 0)) {
      clearBit(fromb);
      if (fromb >= tob)
        return;
      fromb++;
    }

    while ((tob % bitsperword != 0)) {
      clearBit(tob);
      if (tob <= fromb)
        return;
      tob--;
    }
    clearBit(tob);

    IdType from_arrpos = fromb / (8 * (int)sizeof(IdType));
    IdType to_arrpos = tob / (8 * (int)sizeof(IdType));
    memset(&array[from_arrpos], 0,
           (to_arrpos - from_arrpos) * (int)sizeof(IdType));
  }

  inline IdType size() const { return len; }

  inline const IdType *getArray() const { return array; }

  inline IdType countSetBits() const {
    IdType ret = 0;

    for (IdType i = 0; i < arrlen; ++i) {
      IdType val = 0x0000000000000001;
      for (unsigned j = 0; j < sizeof(IdType) * 8; ++j) {
        if (array[i] & val)
          ++ret;

        val = val << 1;
      }
    }

    return std::min(len, ret);
  }

private:
  inline static void bitToPos(IdType b, IdType &arrpos, IdType &bitpos) {
    arrpos = b / (8 * (int)sizeof(IdType));
    bitpos = b & (8 * (int)sizeof(IdType) - 1);
  }

  IdType *array;
  IdType len;
  IdType arrlen;
};

#endif
