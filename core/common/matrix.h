#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "utils.h"
#include <iomanip>

#define PRINT_WIDTH_MATRIX 18
// #define PRINT_WIDTH_MATRIX 10

template <class T>
void multiplyMatrices(const T *A, int A_rows, int A_cols, const T *B,
                      int B_rows, int B_cols, T *C, int C_rows, int C_cols) {
  // perform C = A*B
  // TODO : Add additional checks here
  if (A_cols != B_rows) {
    cout << "ERROR : MATRIX MULTIPLICATION NOT POSSIBLE" << endl;
    return;
  }

  if (A_rows != C_rows) {
    cout << "ERROR : MATRIX MULTIPLICATION NOT POSSIBLE" << endl;
    return;
  }

  if (B_cols != C_cols) {
    cout << "ERROR : MATRIX MULTIPLICATION NOT POSSIBLE" << endl;
    return;
  }

  for (int i = 0; i < C_rows; i++) {
    for (int j = 0; j < C_cols; j++) {
      C[i * C_cols + j] = 0;
      for (int k = 0; k < A_cols; k++) {
        C[i * C_cols + j] += A[i * A_cols + k] * B[k * B_cols + j];
      }
    }
  }
}

template <class T>
void addMatrices(T *A, int A_rows, int A_cols, T *B, int B_rows, int B_cols,
                 T *C, int C_rows, int C_cols) {
  if ((A_rows != B_rows) || (A_rows != C_rows)) {
    cout << "ERROR : MATRIX ADDITION NOT POSSIBLE" << endl;
    return;
  }
  if ((A_cols != B_cols) || (A_cols != C_cols)) {
    cout << "ERROR : MATRIX ADDITION NOT POSSIBLE" << endl;
    return;
  }
  for (int i = 0; i < C_rows; i++) {
    for (int j = 0; j < C_cols; j++) {
      C[i * C_cols + j] = A[i * A_cols + j] + B[i * B_cols + j];
    }
  }
}

template <class T>
void subtractMatrices(T *A, int A_rows, int A_cols, T *B, int B_rows,
                      int B_cols, T *C, int C_rows, int C_cols) {
  if ((A_rows != B_rows) || (A_rows != C_rows)) {
    cout << "ERROR : MATRIX ADDITION NOT POSSIBLE" << endl;
    return;
  }
  if ((A_cols != B_cols) || (A_cols != C_cols)) {
    cout << "ERROR : MATRIX ADDITION NOT POSSIBLE" << endl;
    return;
  }
  for (int i = 0; i < C_rows; i++) {
    for (int j = 0; j < C_cols; j++) {
      C[i * C_cols + j] = A[i * A_cols + j] - B[i * B_cols + j];
    }
  }
}

template <class T>
void scalarMultiplyToMatrix(T *A, int A_rows, int A_cols, T *B, int B_rows,
                            int B_cols, T scalarVal) {

  if ((A_rows != B_rows) || (A_cols != B_cols)) {
    cout << "ERROR : SCALAR MATRIX MULTIPLICATION NOT POSSIBLE" << endl;
    return;
  }

  for (int i = 0; i < B_rows; i++) {
    for (int j = 0; j < B_cols; j++) {
      B[i * B_cols + j] = A[i * A_cols + j] * scalarVal;
    }
  }
}

template <class T>
void getCoFactor(T *A, int A_rows, int A_cols, T *B, int B_rows, int B_cols,
                 int m, int n) {
  // Find co-factor for element A[m][n] amd store the coFactor matrix in B
  // Check if A and B are square matrices
  if ((A_rows != A_cols) || (B_rows != B_cols)) {
    cout << "ERROR : COFACTOR NOT POSSIBLE : Not a square matrix" << endl;
    return;
  }

  int B_i = 0, B_j = 0;

  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      //  If A[i][j] is not in the mth row and the nth column, add it to the
      //  B[B_i][B_j]
      if (i != m && j != n) {
        B[B_i * B_cols + B_j] = A[i * A_cols + j];
        B_j++;
        if (B_j == B_cols) {
          B_j = 0;
          B_i++;
        }
      }
    }
  }
}

template <class T> T determinantOfMatrix(T *A, int A_rows, int A_cols) {
  double result = 0.0;
  if (A_rows != A_cols) {
    cout << "ERROR : DETERMINANT NOT POSSIBLE : Not a square matrix" << endl;
    return 0;
  }

  if (A_rows == 1)
    return A[0];

  T coFactorMatrix[(A_rows - 1) * (A_cols - 1)];

  int currentSign = 1;

  // Iterate for each element of first row
  for (int i = 0; i < A_cols; i++) {
    getCoFactor(A, A_rows, A_cols, coFactorMatrix, A_rows - 1, A_cols - 1, 0,
                i);
    result += (currentSign * A[0 * A_cols + i] *
               determinantOfMatrix(coFactorMatrix, A_rows - 1, A_cols - 1));
    currentSign = -currentSign;
  }

  return result;
}

template <class T>
void adjointOfMatrix(T *A, int A_rows, int A_cols, T *B, int B_rows,
                     int B_cols) {
  if ((A_rows != B_rows) || (A_cols != B_cols)) {
    cout << "ERROR : ADJOINT NOT POSSIBLE : dimension error" << endl;
    return;
  }
  int currentSign = 1;
  T coFactorMatrix[(A_rows - 1) * (A_cols - 1)];

  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      getCoFactor(A, A_rows, A_cols, coFactorMatrix, A_rows - 1, A_cols - 1, i,
                  j);

      currentSign = ((i + j) % 2 == 0) ? 1 : -1;

      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      B[j * A_rows + i] =
          (currentSign) *
          (determinantOfMatrix(coFactorMatrix, A_rows - 1, A_cols - 1));
    }
  }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
template <class T>
void inverseOfMatrix(T *A, int A_rows, int A_cols, T *B, int B_rows,
                     int B_cols) {
  if ((A_rows != B_rows) || (A_cols != B_cols)) {
    cout << "ERROR : INVERSE NOT POSSIBLE : dimension error" << endl;
    return;
  }

  T determinant = determinantOfMatrix(A, A_rows, A_cols);
  if (determinant == 0) {
    cout << "ERROR : INVERSE NOT POSSIBLE : Singular matrix" << endl;
    exit(1);
    return;
  }

  T adjointMatrix[A_rows * A_cols];

  adjointOfMatrix(A, A_rows, A_cols, adjointMatrix, A_rows, A_cols);

  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      B[i * B_cols + j] = adjointMatrix[i * A_cols + j] / determinant;
    }
  }
  return;
}

template <class T>
T *getTransposeProduct(const T *A, int A_rows, int A_cols, T *B, int B_rows,
                       int B_cols) {
  multiplyMatrices(A, A_rows, A_cols, A, A_cols, A_rows, B, B_rows, B_cols);
}

template <class T> void printMatrix(T *A, int A_rows, int A_cols) {
  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      cout << setw(PRINT_WIDTH_MATRIX) << A[i * A_cols + j];
    }
    cout << "\n";
  }
}

template <class T> void printMatrixComma(T *A, int A_rows, int A_cols) {
  for (int i = 0; i < A_rows; i++) {
    for (int j = 0; j < A_cols; j++) {
      cout << A[i * A_cols + j] << ",";
      // return;
    }
  }
}

void additionTest() {
  double A[12] = {1, 2, 3, 2, 4, 5, 6, 2, 7, 8, 9, 4};
  double B[12] = {4, 5, 5, 2, 1, 2, 3, 4, 4, 5, 2, 1};
  double C[12];
  addMatrices(A, 3, 4, B, 3, 4, C, 3, 4);
  cout << "A ";
  printMatrix<double>(A, 3, 4);
  cout << "B ";
  printMatrix<double>(B, 3, 4);
  cout << "A+B ";
  printMatrix<double>(C, 3, 4);

  addMatrices(A, 4, 3, B, 4, 3, C, 4, 3);
  cout << "A ";
  printMatrix<double>(A, 4, 3);
  cout << "B ";
  printMatrix<double>(B, 4, 3);
  cout << "A+B ";
  printMatrix<double>(C, 4, 3);
}

void multiplicationTest() {
  double A[12] = {1, 2, 3, 2, 4, 5, 6, 2, 7, 8, 9, 4};
  double B[12] = {4, 5, 5, 2, 1, 2, 3, 4, 4, 5, 2, 1};
  double C[16];
  double D[9];
  multiplyMatrices(A, 4, 3, B, 3, 4, C, 4, 4);
  cout << "A ";
  printMatrix<double>(A, 4, 3);
  cout << "B ";
  printMatrix<double>(B, 3, 4);
  cout << "AxB ";
  printMatrix<double>(C, 4, 4);

  multiplyMatrices(A, 3, 4, B, 4, 3, D, 3, 3);
  cout << "A ";
  printMatrix<double>(A, 3, 4);
  cout << "B ";
  printMatrix<double>(B, 4, 3);
  cout << "AxB ";
  printMatrix<double>(D, 3, 3);
}

void adjointTest() {
  double A[16] = {5, -2, 2, 7, 1, 0, 0, 3, -3, 1, 5, 0, 3, -1, -9, 4};
  double B[16];
  printMatrix<double>(A, 4, 4);
  adjointOfMatrix(A, 4, 4, B, 4, 4);
  cout << "Adjoint ";
  printMatrix<double>(B, 4, 4);
}

void determinantTest() {
  double A[16] = {1, 0, 2, -1, 3, 0, 0, 5, 2, 1, 4, -3, 1, 0, 5, 0};
  printMatrix<double>(A, 4, 4);
  cout << "Determinant : " << determinantOfMatrix<double>(A, 4, 4) << endl;
}

void inverseTest() {
  double A[9] = {0.000149344612024,  -0.000052379192225, -0.000166212557819,
                 -0.000052379192225, 0.000115858437771,  0.000136789232738,
                 -0.000166212557819, 0.000136789232738,  0.000248242440191};
  double C[9];
  printMatrix<double>(A, 3, 3);
  inverseOfMatrix(A, 3, 3, C, 3, 3);
  cout << "Inverse \n";
  printMatrix<double>(C, 3, 3);
}

void matrixTests() {
  cout << "Matrix test \n\n";
  inverseTest();
  // multiplicationTest();
  // additionTest();
  // adjointTest();
  // determinantTest();
}

#endif