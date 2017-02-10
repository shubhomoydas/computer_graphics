#include "glMatrix.h"
#include <iostream>
#include <iomanip>

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

namespace smd {
	
	// unit test
	void testGaussSeidel() {
		Matrix A(4,4);
		double **a = A.data;
		a[0][0] = 10.0;  a[0][1] = -1.0;  a[0][2] =  2.0; a[0][3] =  0.0;
		a[1][0] = -1.0;  a[1][1] = 11.0;  a[1][2] = -1.0; a[1][3] =  3.0;
		a[2][0] =  2.0;  a[2][1] = -1.0;  a[2][2] = 10.0; a[2][3] = -1.0;
		a[3][0] =  0.0;  a[3][1] =  3.0;  a[3][2] = -1.0; a[3][3] =  8.0;

		Matrix b(1,4);
		double *d_b = b.data[0];
		d_b[0] = 6.0; d_b[1] = 25; d_b[2] = -11; d_b[3] = 15;

		Matrix x0(1,4);
		x0.setZeros();

		Matrix x(1,4);
		x.setZeros();

		gauss_seidel(A, b, x0, x, 10);

		std::cout << "Final result: " << std::endl << x << std::endl;
	}

	double row_X_vec(Matrix& A, double *x, int r) {
		double s = 0.0;
		double *d_x = x;
		double *d_A = A.data[r];
		for (int i = 0; i < A.cols; i++, d_x++, d_A++) {
			s = s + (*d_A) * (*d_x);
		}
		return s;
	}

	// A is a square matrix, b and x0 are row vectors.
	Matrix& gauss_seidel(Matrix& A, Matrix& b, Matrix& x0, Matrix& x, int iters) {
		int n = A.rows;
		//Matrix x(1,x0.cols);
		x.copyFrom(x0);
		double *d_x = x.data[0];
		double *d_b = b.data[0];
		double **d_A = A.data;
		for (int k = 0; k < iters; k++) {
			for (int i = 0; i < n; i++) {
				d_x[i] = (1/d_A[i][i])*(d_b[i] - row_X_vec(A, d_x, i) + d_A[i][i]*d_x[i]);
			}
			//std::cout << "Solution at iter " << k << std::endl << x << std::endl;
		}
		return x;
	}

	int gaussj(double **a, int n, double **b, int m);
	
	// invert a square matrix
	int Matrix::inverse(Matrix& b) {
		if (rows != cols || b.rows != b.cols || rows != rows) {
			std::cout << "invert: Matrix dimensions do not match" << std::endl;
			exit(-1);
		}
		Matrix *m = new Matrix(rows);
		b.copyFrom(*this);
		int ret = gaussj(b.data, b.rows, m->data, m->rows);
		if (ret == -1) {
			//std::cout << "Singular Matrix:" << std::endl << (*this) << std::endl;
			return -1;
		}
		//std::cout << "m=" << std::endl << (*m) << std::endl;
		delete m;
		return 0;
	}

	/*
	Borrowed from William H Press Numerical Recipies in C
	Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[0..n-1][0..n-1]
	is the input matrix. b[0..n-1][0..m-1] is input containing the m right-hand side vectors. On
	output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
	vectors.
	*/
	int gaussj(double **a, int n, double **b, int m) {
		int *indxc,*indxr,*ipiv;
		int i,icol,irow,j,k,l,ll;
		double big,dum,pivinv,temp;

		indxc = new int[n]; //ivector(1,n);
		indxr = new int[n]; //ivector(1,n);
		ipiv =  new int[n]; //ivector(1,n);
		for (j=0;j<n;j++) ipiv[j]=-1;
		for (i=0;i<n;i++) {
			big=0.0;
			for (j=0;j<n;j++)
				if (ipiv[j] != 0)
					for (k=0;k<n;k++) {
						if (ipiv[k] == -1) {
							if (fabs(a[j][k]) >= big) {
								big=fabs(a[j][k]);
								irow=j;
								icol=k;
							}
						}
					}
					++(ipiv[icol]);

					if (irow != icol) {
						for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
						for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
					}
					indxr[i]=irow;
					indxc[i]=icol;
					if (a[icol][icol] == 0.0) {
						//nrerror("gaussj: Singular Matrix");
						//std::cout << "gaussj: Singular Matrix" << std::endl;
						delete[] ipiv; //free_ivector(ipiv,1,n);
						delete[] indxr; //free_ivector(indxr,1,n);
						delete[] indxc; //free_ivector(indxc,1,n);
						return -1;
					}
					pivinv=1.0/a[icol][icol];
					a[icol][icol]=1.0;
					for (l=0;l<n;l++) a[icol][l] *= pivinv;
					for (l=0;l<m;l++) b[icol][l] *= pivinv;

					for (ll=0;ll<n;ll++)
						if (ll != icol) {
							dum=a[ll][icol];
							a[ll][icol]=0.0;
							for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
							for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
						}
		}
		for (l=n-1;l>=0;l--) {
			if (indxr[l] != indxc[l])
				for (k=0;k<n;k++)
					SWAP(a[k][indxr[l]],a[k][indxc[l]]);
		}
		delete[] ipiv; //free_ivector(ipiv,1,n);
		delete[] indxr; //free_ivector(indxr,1,n);
		delete[] indxc; //free_ivector(indxc,1,n);
		return 0;
	}

	std::ostream& Matrix::print(std::ostream& out) const {
		std::ios state(nullptr);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				std::cout << std::setprecision(4) << std::setw(10) << data[i][j];
			std::cout << std::endl;
		}
		std::cout.copyfmt(state);
		return out;
	}

	std::ostream& operator<< (std::ostream& out, const Matrix& matrix) {
		return matrix.print(out);
	}

} // namespace smd
