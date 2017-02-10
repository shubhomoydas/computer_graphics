#ifndef OSUGLMATRIX
#define OSUGLMATRIX

#include <iostream>
#include <stdio.h>
#include <array>

namespace smd {

	class Matrix {
	public:
		double **data;
		int rows, cols;
		Matrix() : data(0), rows(0), cols(0) {}
		Matrix(int dim) {
			rows = dim;
			cols = dim;
			data = createData(rows,cols);
			//setZeros(); // commented for performance
		}

		Matrix(int _rows, int _cols, double **_a):
			rows(_rows), cols(_cols), data(_a) {}
	
		Matrix(int _rows, int _cols):
			rows(_rows), cols(_cols) {
			data = createData(rows, cols);
		}
	
		Matrix* clone() {
			Matrix *cl = new Matrix(rows, cols);
			cl->copyData(data);
			return cl;
		}

		Matrix& setIdentity() {
			setZeros();
			setDiag(1);
			return *this;
		}

		Matrix& copyFrom(const Matrix& src) {
			if (!(rows == src.rows && cols == src.cols)) {
				std::cout << "Illegal copy. Dimensions do not match." << std::endl;
				exit(-1);
			}
			copyData(src.data);
			return *this;
		}

		/**
		 * Only works for square matrices
		 *
		 * this[i,i] = val
		 */
		Matrix& setDiag(double val) {
			if (rows != cols) {
				std::cout << "Need rows==cols for setDiag()" << std::endl;
				exit(-1);
			}
			for (int i = 0; i < rows; i++) {
				data[i][i] = val;
			}
			return *this;
		}
		/** this[i,i] = vals[i] */
		Matrix& setDiag(double* vals, int length) {
			for (int i = 0; i < rows && i < length; i++) {
				data[i][i] = vals[i];
			}
			return *this;
		}

		/** this[row,col] = val */
		Matrix& set(int row, int col, double val) {
			data[row][col] = val;
			return *this;
		}

		/** this[.,.] = this[.,.] * d */
		Matrix& mul(double d) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++)
					data[i][j] *= d;
			}
			return *this;
		}
	
		/** this[.,.] = this[.,.] / d */
		Matrix& div(double d) {
			return mul(1 / d);
		}
	
		/** this[.,.] = this[.,.] + d */
		Matrix& add(double d) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++)
					data[i][j] += d;
			}
			return *this;
		}
	
		/** this = this + B*scale */
		Matrix& add(Matrix& b, double scale) {
			if (rows == b.rows && cols == b.cols) {
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++)
						data[i][j] += scale*b.data[i][j];
				}
			} else {
				printf("Matrix dimensions do not match in add()!\n");
				exit(-1);
			}
			return *this;
		}
	
		/** this = this + B */
		Matrix& add(Matrix& b) {
			return add(b, 1);
		}

		/** this = this - B */
		Matrix& sub(Matrix& b) {
			return add(b, -1);
		}

		/** dest = transpose(this) */
		Matrix& transpose(Matrix& dest) {
			if (rows != dest.cols || cols != dest.cols) {
				printf("Matrix dimensions do not match in transpose()!\n");
				exit(-1);
			}
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					dest.data[j][i] = data[i][j];
				}
			}
			return dest;
		}
	
		/** this = transpose(this) */
		Matrix& transpose() {
			double **newData = createData(cols, rows);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					newData[j][i] = data[i][j];
				}
			}
			free();
			data = newData;
			int temp = rows;
			rows = cols; cols = temp;
			return *this;
		}
	
		/** dest = this * B */
		Matrix& mul(Matrix& b, Matrix& dest) {
			if (!(cols == b.rows && rows == dest.rows && b.cols == dest.cols)) {
				printf("Matrix dimensions do not match in mul()!\n");
				exit(-1);
			}
			int p = b.cols;
			double **bData = b.data;
			double **destData = dest.data;
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < p; j++) {
					destData[i][j] = 0;
					for (int k = 0; k < cols; k++) {
						destData[i][j] += (data[i][k] * bData[k][j]);
					}
				}
			}
			return dest;
		}

		/** this = this * B */
		Matrix& mul(Matrix& b) {
			Matrix dest(rows, cols);
			mul(b, dest);
			copyData(dest.data);
			return *this;
		}
	
		Matrix& setVal(const double val) {
			int l = rows*cols;
			if (l==0) return *this;
			double *ptr = data[0];
			for (int i = 0; i < l; i++, ptr++) {
				*ptr = val;
			}
			return *this;
		}

		Matrix& setZeros() {
			int l = rows*cols;
			if (l==0) return *this;
			double *ptr = data[0];
			for (int i = 0; i < l; i++, ptr++) {
				*ptr = 0;
			}
			return *this;
		}

		int inverse(Matrix& b);

		~Matrix() {
			free();
		}
		std::ostream& print(std::ostream& out) const;
		friend std::ostream& operator<< (std::ostream& stream, const Matrix& matrix);
	private:
		double **copyData(double **source) {
			double **src = source;
			double **dst = data;
			for (int i = 0; i < rows; i++, src++, dst++) {
				double *s = *src;
				double *d = *dst;
				for (int j = 0; j < cols; j++, s++, d++) *d=*s;
			}
			return data;
		}
		double **createData(int rows, int cols) {
			double **data = new double*[rows];
			double *tdata = new double[rows*cols];
			for (int i = 0; i < rows; i++) {
				data[i] = tdata;
				tdata = tdata+cols;
			}
			return data;
		}
		void free() {
			if (data) {
				delete[] data[0];
				delete[] data;
			}
			data = 0;
		}

	};

	Matrix& gauss_seidel(Matrix& A, Matrix& b, Matrix& x0, Matrix& x, int iters);
	void testGaussSeidel();

} // namespace smd

#endif OSUGLMATRIX
