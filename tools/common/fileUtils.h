// Modifications Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
// 
// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
// 
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "../../core/common/utils.h"
#include <fstream>
#include <stdlib.h>

// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  unsigned long n; // total number of characters
  char *Chars;     // array storing all strings
  unsigned long m; // number of substrings
  char **Strings;  // pointers to strings (all should be null terminated)
  words() {}
  words(char *C, unsigned long nn, char **S, unsigned long mm)
      : Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {
    free(Chars);
    free(Strings);
  }
};

inline bool isSpace(char c) {
  switch (c) {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ':
    return true;
  default:
    return false;
  }
}

// parallel code for converting a string to words
words stringToWords(char *Str, unsigned long n) {
  {
    parallel_for(unsigned long i = 0; i < n; i++) if (isSpace(Str[i])) Str[i] =
        0;
  }

  // mark start of words
  bool *FL = newA(bool, n);
  FL[0] = Str[0];
  {
    parallel_for(unsigned long i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1];
  }

  // offset for each start of word
  _seq<unsigned long> Off = sequence::packIndex<unsigned long>(FL, n);
  unsigned long m = Off.n;
  unsigned long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char *, m);
  { parallel_for(unsigned long j = 0; j < m; j++) SA[j] = Str + offsets[j]; }

  free(offsets);
  free(FL);
  return words(Str, n, SA, m);
}

_seq<char> readStringFromFile(char *fileName) {
  ifstream file(fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  unsigned long end = file.tellg();
  file.seekg(0, ios::beg);
  unsigned long n = end - file.tellg();
  char *bytes = newA(char, n + 1);
  file.read(bytes, n);
  file.close();
  return _seq<char>(bytes, n);
}
