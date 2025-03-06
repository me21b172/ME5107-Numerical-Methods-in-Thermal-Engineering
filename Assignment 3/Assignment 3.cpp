#include<bits/stdc++.h>
using namespace std;
vector<vector<double>> matrix_inverse(vector<vector<double>> A){
    // "There are some issues. Go through the code again"
    int n = A.size();
    vector<vector<double>> B(n,vector<double>(n,0.0));
    for(int i=0;i<n;i++)
        B[i][i] = 1.0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i!=j) {
                double ratio = A[j][i]/A[i][i];
                for(int k=0;k<n;k++){
                    B[j][k] = B[j][k] - B[i][k]*(ratio);
                    A[j][k] = A[j][k] - A[i][k]*(ratio);
                }
            }
            else{
                double ratio = A[i][i];
                for(int k=0;k<n;k++){
                    A[i][k] /= ratio;
                    B[i][k] /= ratio;
                }
            }

        }
    }
    return B;
}
vector<vector<double>> matrix_multiplication(vector<vector<double>>& A, vector<vector<double>>& B){
    vector<vector<double>> C(A.size(),vector<double>(B[0].size(),0.0));
    for(int i=0;i<A.size();i++){
        for(int j=0;j<B[0].size();j++){
            for(int k=0;k<B.size();k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}
vector<vector<double>> matrix_subtraction(vector<vector<double>>& A, vector<vector<double>>& B){
    vector<vector<double>> C(A.size(),vector<double>(A[0].size(),0.0));
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A[0].size();j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

void print_matrix(vector<vector<double>>& A){
    for(int i=0;i<A.size();i++){
        for(int j=0;j<A[0].size();j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}
vector<vector<double>> matrix_D(int N){
    vector<vector<double>> D(N,vector<double>(N,0.0));
    for(int i=0;i<N-1;i++){
        D[i][i] = -4.0;
        D[i][i+1] = 1.0;
        D[i+1][i] = 1.0;
    }
    D[N-1][N-1] = -4.0;
    return D;
}
vector<vector<double>> matrix_A(int N){
    vector<vector<double>> D(N,vector<double>(N,0.0));
    for(int i=0;i<N;i++){
        D[i][i] = 1.0;
    }
    return D;
}
vector<vector<double>> matrix_B(int N){
    vector<vector<double>> D(N,vector<double>(N,0.0));
    for(int i=0;i<N;i++){
        D[i][i] = 1.0;
    }
    return D;
}
vector<double> vector_f(int N){
    double Nblk = (double)N;
    vector<double> f(N,1/Nblk);
    f[0] = 1.0;f[N-1] = 2.0;
    return f;
}
vector<vector<double>> thomas_algorithm(vector<vector<vector<vector<double>>>>& A,vector<vector<double>>& f){
    // Forward Elimination
    int n = A.size();
    for(int i=1;i<n;i++){
        vector<vector<double>> inv_D = matrix_inverse(A[i-1][i-1]);
        vector<vector<double>> multipler = matrix_multiplication(A[i][i-1],inv_D),f_i_minus_1 = {f[i-1]};
        vector<vector<double>> sub = matrix_multiplication(multipler,A[i-1][i]),temp = matrix_multiplication(multipler, f_i_minus_1);
        A[i][i] = matrix_subtraction(A[i][i],sub);
        f[i] = temp[0];
    }
    // Backward Substitution
    vector<vector<double>> x(n,vector<double>(f[0].size(),0.0)),inv_D = matrix_inverse(A[n-1][n-1]),f_i_minus_1 = {f[n-1]};
    vector<vector<double>> temp = matrix_multiplication(inv_D,f_i_minus_1);
    x[n-1] = temp[0];
    for(int i=n-2;i>=0;i--){
        vector<vector<double>> inv_D = matrix_inverse(A[i][i]),x_i_plus_1 = {x[i+1]};
        vector<vector<double>> sub = matrix_multiplication(A[i][i+1],x_i_plus_1) ,f_i_minus_1 = {f[i]};
        sub = matrix_subtraction(f_i_minus_1,sub);
        temp = matrix_multiplication(inv_D,sub);
        x[i] = temp[0];
    }
    return x;
}
vector<vector<double>> tridaigonal_matrix_solver(int N,int Nblk){
    vector<vector<vector<vector<double>>>> A(N,vector<vector<vector<double>>>(N,vector<vector<double>>(Nblk,vector<double>(Nblk,0.0))));
    vector<vector<double>> f(N,vector<double>(Nblk));
    for(int i=0;i<N-1;i++){
        A[i][i] = matrix_D(Nblk);
        A[i+1][i] = matrix_B(Nblk);
        A[i][i+1] = matrix_A(Nblk);
        f[i] = vector_f(Nblk);
    }
    A[N-1][N-1] = matrix_D(Nblk);
    f[N-1] = vector_f(Nblk);
    vector<vector<double>> x = thomas_algorithm(A,f);
    return x;
}
int main(){
    // freopen("output.txt","w",stdout);
    int Nblk = 5;
    vector<int> N = {10, 20, 30};
    
    // Loop over each value in N
    for (int k = 0; k < N.size(); k++){
        vector<vector<double>> x = tridaigonal_matrix_solver(N[k], Nblk);
        // cout << "Solution matrix dimensions: " << x.size() << " x " << x[0].size() << endl;
        
        // Example: print the middle block of the solution
        int midBlock = N[k] / 2; // or another valid index within x's range
        if(midBlock < x.size()){
            for(int j = 0; j < Nblk; j++){
                cout << x[midBlock][j] << " ";
            }
            cout << endl;
        } else {
            cout << "Invalid block index" << endl;
        }
    }
    return 0;
}