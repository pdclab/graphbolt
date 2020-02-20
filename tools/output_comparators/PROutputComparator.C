// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
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

#include "../common/fileUtils.h"
#include "../../core/common/parallel.h"
#include "../../core/common/parseCommandLine.h"
#include "../../core/common/utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

struct CFData {
  long vertex;
  long inDegree;
  long outDegree;
  double pr;

  CFData() = default;

  ~CFData() {}

  inline void setVertex(string str) { vertex = stol(str); }
  inline void setInDegree(string str) { inDegree = stol(str); }
  inline void setOutDegree(string str) { outDegree = stol(str); }
  inline void setPR(string str) { pr = stod(str); }
  inline void setIndex(char *str, int index) {
    switch (index) {
    case 0:
      setVertex(string(str));
      break;
    case 1:
      setInDegree(string(str));
      break;
    case 2:
      setOutDegree(string(str));
      break;
    case 3:
      setPR(string(str));
      break;
    default:
      cout << "ERROR: Incorrect index (" << index << ")\n";
      break;
    }
  }
  string printString() {
    ostringstream streamObj;
    streamObj << setprecision(VAL_PRECISION2);
    streamObj << vertex << " " << inDegree << " " << outDegree << " " << pr;
    return streamObj.str();
  }

  bool isEqualWithThreshold(const CFData &rhs, double threshold) {
    bool ret = true;
    if (fabs(pr - rhs.pr) > threshold) {
      ret = false;
    }
    return ret;
  }
};

int parallel_main(int argc, char *argv[]) {
  commandLine P(argc, argv, "");
  char *file1_path = P.getOptionValue("-file1");
  char *file2_path = P.getOptionValue("-file2");
  char *output_file_path = P.getOptionValue("-ofile");
  bool print_output = (output_file_path != NULL);
  long thresholds = P.getOptionLongValue("-thresholds", 6);

  cout << "file1 : " << file1_path << "\n";
  cout << "file2 : " << file2_path << "\n";
  cout << "output_file : " << output_file_path << "\n";
  cout << "print_output : " << print_output << "\n";

  words W1, W2;
  cout << "Reading file 1\n";
  _seq<char> S1 = readStringFromFile(file1_path);
  cout << "String to words\n";
  W1 = stringToWords(S1.A, S1.n);

  cout << "Reading file 2\n";
  _seq<char> S2 = readStringFromFile(file2_path);
  cout << "String to words\n";
  W2 = stringToWords(S2.A, S2.n);
  cout << "Finished reading\n";

  // number of words
  long nWords1 = W1.m;
  long nWords2 = W2.m;

  // number of characters
  long nChars1 = W1.n;
  long nChars2 = W2.n;
  cout << fixed;
  //   cout << setprecision(VAL_PRECISION2);

  if (nChars1 != nChars2) {
    "WARNING: Number of characters in both files not equal\n";
  }

  if (nWords1 != nWords2) {
    "ERROR: Number of words in both files not equal\n";
    W1.del();
    W2.del();
    exit(1);
  }

  long words_per_line = 4;
  cout << "words_per_line : " << words_per_line << "\n";

  if (nWords1 % words_per_line != 0) {
    cout << "ERROR: Line format incorrect. Expected " << words_per_line
         << " words per line\n";
    W1.del();
    W2.del();
    exit(1);
  }

  long lines = nWords1 / words_per_line;

  CFData *file1_data = new CFData[lines];
  CFData *file2_data = new CFData[lines];

  parallel_for(long i = 0; i < nWords1; i++) {
    file1_data[i / words_per_line].setIndex(W1.Strings[i], i % words_per_line);
    file2_data[i / words_per_line].setIndex(W2.Strings[i], i % words_per_line);
  }

  ofstream *outputFiles = new ofstream[thresholds];
  double *threshold_array = new double[thresholds];
  long *incorrect_count = new long[thresholds];
  double temp2 = 1.0;
  bool **different_flag = new bool *[thresholds];
  for (int i = 0; i < thresholds; i++) {
    threshold_array[i] = temp2;
    temp2 = temp2 / 10.0;
    different_flag[i] = new bool[lines];
    incorrect_count[i] = 0;
    if (print_output) {
      std::ostringstream path_stream;
      path_stream.precision(thresholds);
      path_stream << std::fixed;
      path_stream << output_file_path << "_" << threshold_array[i];
      outputFiles[i].open(path_stream.str(), ios::out);
      outputFiles[i] << fixed;
      outputFiles[i] << setprecision(VAL_PRECISION2);
    }
  }

  parallel_for(long currLine = 0; currLine < lines; currLine++) {
    for (int i = 0; i < thresholds; i++) {
      if (file1_data[currLine].isEqualWithThreshold(file2_data[currLine],
                                                    threshold_array[i])) {
        different_flag[i][currLine] = false;
      } else {
        different_flag[i][currLine] = true;
      }
    }
  }

  for (long currLine = 0; currLine < lines; currLine++) {
    for (int i = 0; i < thresholds; i++) {
      if (different_flag[i][currLine] == true) {
        incorrect_count[i]++;
        if (print_output) {
          outputFiles[i] << file1_data[currLine].printString() << " \t "
                         << file2_data[currLine].printString() << "\n";
        }
      }
    }
  }
  cout << setprecision(thresholds);
  for (int i = 0; i < thresholds; i++) {
    cout << threshold_array[i] << " : " << incorrect_count[i] << "\n";
  }

  for (int i = 0; i < thresholds; i++) {
    if (print_output) {
      outputFiles[i].close();
    }
    delete[] different_flag[i];
  }
  W1.del();
  W2.del();
  delete[] different_flag;
  delete[] outputFiles;

  return 0;
}