#include<bits/stdc++.h>
using namespace std;
vector<vector<double>> matrix_inverse(vector<vector<double>>& A){
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
vector<vector<double>> tridaigonal_matrix_solution(double N){

}
int main(){
    double N;
    // cin>>N;
    vector<vector<double>> x = tridaigonal_matrix_solution(N);
    // vector<vector<double>> b ={ {1,2,3},{4,5,6},{7,8,9}};
    vector<vector<double>> b = {{1,2,3},{5,6,1},{2,1,9}};
    vector<vector<double>> A = matrix_inverse(b);
    vector<vector<double>> C = matrix_multiplication(A,b);
    // for(int i=0;i<C.size();i++){
    //     for(int j=0;j<C[0].size();j++){
    //         cout<<C[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    for(int i=0;i<C.size();i++){
        for(int j=0;j<C[0].size();j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
    
}