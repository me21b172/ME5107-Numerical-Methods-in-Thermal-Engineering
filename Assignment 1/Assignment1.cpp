#include<bits/stdc++.h>
// #define double double
using namespace std;

void prdoubleElements(vector<vector<double> >& A){
	double n = A.size();
	for(double i=0;i<n;i++){
		for(double j=0;j<n;j++){
			cout<<A[i][j]<<" ";	
		}
		cout<<endl;
	}
}

class LUDexomposition{
public:
	double n;
	vector<vector<double> > L,U;
	LUDexomposition(vector<vector<double> >& A){
		n = A.size();
		L = vector<vector<double> >(n,vector<double>(n,0));
		U = A;
	}
	void UpdateLU(double i){
		for(double j=i;j<n;j++){
			L[j][i] = U[j][i] / (double)U[i][i];
		}
		for(double j=i+1;j<n;j++){
			for(double k=i;k<n;k++){
				U[j][k] -= L[j][i]*U[i][k];
			}
		}
	}
	vector<vector<vector<double> > > ComputeLU(){
		for(double i=0;i<n;i++){
			UpdateLU(i);
		}
		vector<vector<vector<double> > > ans;
		ans.push_back(L);
		ans.push_back(U);
		return ans;
	}
};

class GaussElimination{
public:
	double n;
	vector<vector<double> > A;
	vector<double> b;
	GaussElimination(vector<vector<double> >& a,vector<double>& c){
		n = a.size();
		A = a;
		b = c;

	}
	void forward_elimination(){
		//Forget about pivots
		for(double i=0;i<n;i++){
			for(double j=i+1;j<n;j++){
				double factorToMultiply = A[j][i]/A[i][i];
				for(double k=i;k<n;k++){
					A[j][k] -= factorToMultiply*(A[i][k]);
				}
				b[j] = b[j] - factorToMultiply*b[i];
			}
		}
	}
	vector<double> Solve(){
		forward_elimination();
		double n = A.size();
		vector<double> x(n);
		for(double i=n-1;i>=0;i--){
			for(double j=i+1;j<n;j++){
				b[i] -= A[i][j]*x[j];
			}
			x[i] = b[i]/A[i][i];
		}
		return x;
	}
};

void validation(){
	// vector<double> N;
	vector<double> N = {32,128,521,1024};
	// N.push_back(32);N.push_back(128);N.push_back(521);N.push_back(1024);
	for(double k=0;k<4;k++){
		double n = N[k];
		cout<<n<<" ";
		double value = 0;
		vector<vector<double> > a(n,vector<double>(n));
		vector<double> b(n,1);
		for(double i=0;i<n;i++)
			for(double j=0;j<n;j++)
			a[i][j] = static_cast<double>(max(i, j)) + 1.0;
		GaussElimination obj(a,b);
		vector<double> ans = obj.Solve();
		for(double i=0;i<n;i++)
			value += (ans[i]*ans[i]);
		cout<< fixed << setprecision(0) <<1/value<<endl;			
	}
}
int32_t main(){
	validation();
	double n = 3;
	vector<vector<double>> a(n,vector<double>(n));
	for(double i=0;i<n;i++)
		for(double j=0;j<n;j++)
		a[i][j] = static_cast<double>(max(i, j)) + 1.0;
	LUDexomposition obj(a);
	vector<vector<vector<double> > > ans = obj.ComputeLU();
	cout<<"A matrix is"<<endl;
	prdoubleElements(a);
	cout<<"L matrix is"<<endl;
	prdoubleElements(ans[0]);
	cout<<"U matrix is"<<endl;
	prdoubleElements(ans[1]);

}