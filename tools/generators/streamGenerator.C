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

#include "../../core/common/parseCommandLine.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#define newA(__E, __n) (__E *)malloc((__n) * sizeof(__E))
#define renewA(__E, __T, __n) (__E *)realloc(__T, (__n) * sizeof(__E))

// A simple stream generator to mock a stream of edges arriving in a pipe
int main(int argc, char **argv) {
  commandLine config(argc, argv, "");
  string output_pipe = config.getOptionValue("-outputPipe", "/tmp/");
  string edge_operations_file =
      config.getOptionValue("-edgeOperationsFile", "/tmp/");
  if (output_pipe.compare("/tmp/") == 0) {
    cout << "Incorrect arguments. Missing value for \"-outputPipe\"\n";
  }
  if (edge_operations_file.compare("/tmp/") == 0) {
    cout << "Incorrect arguments. Missing value for \"-edgeOperationsFile\"\n";
  }

  ofstream named_pipe;
  named_pipe.open(output_pipe, ios::binary);

  ifstream input_file;
  input_file.open(edge_operations_file);

  int num_lines;
  string line;
  string temp;
  bool bad_input = false;

  cout << "Enter a value equal to or less than zero to end program" << endl;
  cout << "Enter number of lines to send: " << endl;
  cin >> num_lines;
  if (num_lines <= 0) {
    named_pipe.close();
    input_file.close();
    exit(1);
  }

  // Create a pipe and stream in the 'k' edge operations from the
  // "edge_operations_file" to "output_pipe" Once 'k' edges have been streamed
  // to the pipe. Get the size of the next batch 'k'. If there are less than 'k'
  // edges in the edge_operations file, the streamGenerator will detect EOF and
  // exit.
  do {
    for (int i = 0; i < num_lines; ++i) {
      std::getline(input_file, line);
      if ((line.length() == 1) || (line[0] == '%') || (line[0] == '#')) {
        continue;
      }
      if (!input_file.good()) {
        cout << "End of Input File reached" << endl;
        bad_input = true;
        break;
      }
      named_pipe << line << endl;
    }
    if (bad_input == false) {
      cout << "Enter number of lines to send: " << endl;
      cin >> num_lines;
    }
  } while (bad_input == false && num_lines > 0);

  named_pipe.close();
  input_file.close();
  return 0;
}
