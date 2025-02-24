#include<bits/stdc++.h>
using namespace std;

long double Q1(long double n, long double n_1) {
    return 3.0 * n_1 - 2.0 * n;
}

long double f(long double x) {
    return sin(5.0 * x);
}

long double X(long double x, long double n) {
    return (3.0 * x) / n;
}

vector<vector<long double>> create_matrix_A(vector<vector<long double>> &L, vector<vector<long double>> &U) {
    int n = L.size();
    vector<vector<long double>> A(n, vector<long double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            for (int k = 0; k < i; k++) {
                A[i][j] += L[i][k] * U[k][j];
            }
        }
    }
    return A;
}

vector<vector<long double>> create_matrix(long double n) {
    int m = (int)n;
    vector<vector<long double>> A(m, vector<long double>(m, 0));
    A[0][0] = 1.0;
    A[0][1] = 2.0;
    A[m - 1][m - 1] = 1.0;
    A[m - 1][m - 2] = 2.0;
    for (int i = 1; i < m - 1; i++) {
        A[i][i - 1] = 1.0;
        A[i][i] = 4.0;
        A[i][i + 1] = 1.0;
    }
    return A;
}

vector<long double> create_b(long double n) {
    int m = (int)n;
    vector<long double> b(m, 0.0);
    b[0] = (n / 3.0) * ((-5.0 / 2.0) * f(X(0.0, n)) + 2.0 * f(X(1.0,n)) + (1.0 / 2.0) * f(X(2.0, n)));
    b[m - 1] = (n / 3.0) * ((5.0 / 2.0) * f(X(n, n)) - 2.0 * f(X(1.0,n)) - (1.0 / 2.0) * f(X(n - 2.0, n)));
    for (int t = 1; t < m - 1; t++) {
        long double x = X(t, n);
        b[t] = n * (f(X(t + 1.0, n)) - f(X(t - 1.0, n)));
    }
    return b;
}

void doLittle(vector<vector<long double>> &A, vector<vector<long double>> &L, vector<vector<long double>> &U, int ind) {
    int n = A.size();
    L[ind][ind] = 1.0;
    for (int i = ind; i < n; i++) {
        long double ans = 0.0;
        for (int j = 0; j < ind; j++) {
            ans += L[ind][j] * U[j][i];
        }
        U[ind][i] = (A[ind][i] - ans) / L[ind][ind];
    }
    for (int i = ind + 1; i < n; i++) {
        long double ans = 0.0;
        for (int j = 0; j < ind; j++) {
            ans += L[i][j] * U[j][ind];
        }
        L[i][ind] = (A[i][ind] - ans) / U[ind][ind];
    }
}

vector<long double> forward_substitution(vector<vector<long double>> &A, vector<long double> &b) {
    int n = A.size();
    vector<long double> x(n, 0.0);
    for (int i = 0; i < n; i++) {
        long double ans = 0.0;
        for (int j = 0; j < i; j++) {
            ans += A[i][j] * x[j];
        }
        x[i] = (b[i] - ans) / A[i][i];
    }
    return x;
}

vector<long double> backward_substitution(vector<vector<long double>> &A, vector<long double> &b) {
    int n = A.size();
    vector<long double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        long double ans = 0.0;
        for (int j = i + 1; j < n; j++) {
            ans += A[i][j] * x[j];
        }
        x[i] = (b[i] - ans) / A[i][i];
    }
    return x;
}

vector<long double> LUDecomposition(vector<vector<long double>> &A, vector<long double> &b) {
    int m = A.size();
    
    vector<vector<long double>> L(m, vector<long double>(m, 0.0)), U(m, vector<long double>(m, 0.0));
    
    for (int i = 0; i < m; i++) {
        doLittle(A, L, U, i);
    }

    cout << "L Matrix:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cout << L[i][j] << " ";
        }
        cout << endl;
    }

    cout << "U Matrix:\n";
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cout << U[i][j] << " ";
        }
        cout << endl;
    }

    vector<long double> x = forward_substitution(L, b);
    
    vector<long double> y = backward_substitution(U, x);

    return y;
}

vector<long double> Thomas_algorithm(vector<vector<long double>> &A, vector<long double> &b) {
    int n = A.size();
    vector<long double> x(n, 0.0);

    // Forward elimination
    for (int i = 1; i < n; i++) {
        long double factor = A[i][i - 1] / A[i - 1][i - 1];
        A[i][i] -= factor * A[i - 1][i];
        b[i] -= factor * b[i - 1];
        A[i][i - 1] = 0.0; // Set lower diagonal to zero
    }

    // Backward substitution
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (b[i] - A[i][i + 1] * x[i + 1]) / A[i][i];
    }

    return x;
}

int main() {
    long double n = 65.0, constant = 2.9689;

    // Example initialization for Q1 and c
    vector<long double> a(65, constant);
    for (int i = 2; i < a.size(); i++) {
        a[i] = Q1(a[i - 2], a[i - 1]);
    }

    // Print last value of 'a'
    cout << "a[n-1]: " << a[a.size() - 1] << endl;

    // Re-initialize with different constant
    constant = 2.96875;
    vector<int> c(65, constant);
    for (int i = 2; i < c.size(); i++) {
        c[i] = Q1(c[i - 2], c[i - 1]);
    }

    // Print last value of 'c'
    cout << "c[n-1]: " << c[c.size() - 1] << endl;

    // Matrix and vector creation
    n = 4.0;
    vector<vector<long double>> A = create_matrix(n);
    vector<long double> b = create_b(n);

    // Print matrix A
    cout << "Matrix A:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    // Print vector b
    cout << "Vector b:\n";
    for (int i = 0; i < n; i++) {
        cout << b[i] << " ";
    }
    cout << endl;

    // Solve using LU Decomposition
    vector<long double> x_lu = LUDecomposition(A, b);
    
    cout << "Solution using LU Decomposition:\n";
    for (int i = 0; i < n; i++) {
        cout << x_lu[i] << " ";
    }
    cout << endl;

    // Solve using Thomas Algorithm
    vector<long double> x_thomas = Thomas_algorithm(A, b);
    
    cout << "Solution using Thomas Algorithm:\n";
    for (int i = 0; i < n; i++) {
        cout << x_thomas[i] << " ";
    }
    cout << endl;

    return 0;
}

