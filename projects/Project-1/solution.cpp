// Name: Robby Lawrence
// StudentID: 000691931
// NetID: rlawren9
// Description: program uses A = LU decomposition to find the inverse and determinant of the 6 x 6 Hilbert matrix

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

vector<vector<long double> > generate_hilbert(int n) {
  vector<vector<long double> > H(n, vector<long double>(n));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            H[i][j] = 1.0 / (i + j + 1);
        }
    }
    return H;
}

long double lu_decomp(vector<vector<long double> >& A, int n) {
  // we calculate the determinant and return it, modifying A as we gos
  long double det = 1.0;
  for(int i = 0; i < n; i++) {
    // handle the upper triangular part first
    for(int k = i; k < n; k++) {
      long double sum = 0;
      for(int j = 0; j < i; j++) {
        sum += (i > j ? A[i][j] : 1.0) * A[j][k];
      }
      A[i][k] = A[i][k] - sum;
      if(i == k) det *= A[i][i];
    }

    // Lower triangular part (excluding diagonal)
    for(int k = i + 1; k < n; k++) {
      long double sum = 0;
      for(int j = 0; j < i; j++) {
        sum += A[k][j] * A[j][i];
      }
      A[k][i] = (A[k][i] - sum) / A[i][i];
    }
  }
  return det;
}

// function for forward substitution
vector<long double> forward_sub(const vector<vector<long double> >& A, const vector<long double>& b, int n) {
    vector<long double> y(n);
    for(int i = 0; i < n; i++) {
        long double sum = 0;
        for(int j = 0; j < i; j++) {
            sum += A[i][j] * y[j];
        }
        y[i] = b[i] - sum;
    }
    return y;
}

// function for back substitution
vector<long double> back_sub(const vector<vector<long double> >& A, const vector<long double>& y, int n) {
    vector<long double> x(n);
    for(int i = n - 1; i >= 0; i--) {
        long double sum = 0;
        for(int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / A[i][i];
    }
    return x;
}

vector<vector<long double> > multiply_matrices(const vector<vector<long double> >& A, const vector<vector<long double> >& B, int n) {
    vector<vector<long double> > result(n, vector<long double>(n, 0));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

int main() {
    const int n = 6;
    vector<vector<long double> > A = generate_hilbert(n);
    vector<vector<long double> > original_A = A; // save original A for later

    cout << fixed << setprecision(16);
    cout << "A, 6 x 6 Hilbert matrix:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';
    // compute and print LU factorization stored in A
    long double det = lu_decomp(A, n);
    cout << "Modified A matrix (L stored in strictly lower triangular part, U in upper triangular part):\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';
    // print out the previously calculated determinant
    cout << scientific << "Determinant of A: " << det << '\n' << '\n';
    cout << fixed;
    // use functions to find Ainv
    vector<vector<long double> > Ainv(n, vector<long double>(n));
    for(int j = 0; j < n; j++) {
        vector<long double> b(n, 0);
        b[j] = 1;  // matrix column
        vector<long double> y = forward_sub(A, b, n);
        vector<long double> x = back_sub(A, y, n);
        // store each solution as a column of Ainv
        for(int i = 0; i < n; i++) {
            Ainv[i][j] = x[i];
        }
    }

    cout << "A inverse:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << Ainv[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';

    // output B = A * Ainv
    vector<vector<long double> > B = multiply_matrices(original_A, Ainv, n);
    cout << "B = A * Ainv:\n";
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << B[i][j] << " ";
        }
        cout << '\n';
    }

    return 0;
}
