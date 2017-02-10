#include <iostream>
#include "glDraw.h"

int quaternionTests(int argc , char **argv);
int matrixTests(int argc , char **argv);
int lookAtTests(int argc , char **argv);
int inverseTests(int argc , char **argv);

int __main(int argc , char **argv) {
	//return matrixTests(argc, argv);
	//return quaternionTests(argc, argv);
	//return lookAtTests(argc, argv);
	return inverseTests(argc, argv);
}

int inverseTests(int argc , char **argv) {
	smd::Matrix m1(4), m2(4);
	m1.setIdentity();
	m1.inverse(m2);
	std::cout << "Matrix: " << std::endl << m1 << std::endl;
	std::cout << "Inverse: " << std::endl << m2 << std::endl;

	smd::Matrix m3(2), m4(2);
	m3.set(0,0,2); m3.set(0,1,3);
	m3.set(1,0,4); m3.set(1,1,5);
	m3.inverse(m4);
	std::cout << "Matrix: " << std::endl << m3 << std::endl;
	std::cout << "Inverse: " << std::endl << m4 << std::endl;

	smd::Matrix m5(3), m6(3);
	double **d = m5.data;
	d[0][0] = 1; d[0][1] = 0; d[0][2] = 0;
	d[1][0] = 2; d[1][1] = 1; d[1][2] = 0;
	d[2][0] = 0; d[2][1] = 2; d[2][2] = 3;
	std::cout << "Matrix: " << std::endl << m5 << std::endl;
	m5.inverse(m6);
	std::cout << "Inverse: " << std::endl << m6 << std::endl;

	return 0;
}

int lookAtTests(int argc , char **argv) {
	double from[3] = {0,0,0};
	double at[3]   = {1,1,1};
	double up[3]   = {0,1,0};
	smd::GraphicsManager gfxManager;
	//gfxManager.getGraphicsContext().setDebug(true, true);
	gfxManager.lookAt(from, at, up);
	return 0;
}

int quaternionTests(int argc , char **argv) {
	smd::Quaternion q;
	smd::Matrix M(4);

	q = smd::getRotationQuaternion(120, 1, 1, 1, q);
	std::cout << "Quaternion: " << q << std::endl;

	M = smd::getQuaternionToRotationMatrix(q, M);
	std::cout << "Matrix: \n" << M << std::endl;

	return 0;
}

int matrixTests(int argc , char **argv)
{
	smd::Matrix m1(4), m2(4), m3(4);
	m1.setIdentity(); m2.setIdentity(); m3.setIdentity();
	m2.setDiag(3);
	std::cout << "M2 is: \n" << m2 << std::endl;
	m2.set(1,2,5);
	std::cout << "M2 with (1,2)=5 is: \n" << m2 << std::endl;
	m2.transpose();
	std::cout << "M2 transposed is: \n" << m2 << std::endl;
	std::cout << "M1 is: \n" << m1 << std::endl;
	m2.add(m1);
	std::cout << "M2 = M1+M2 is: \n" << m2 << std::endl;
	m2.mul(m1,m3);
	std::cout << "M2*M1 is: \n" << m3 << std::endl;
	m1.mul(m2,m3);
	std::cout << "M1*M2 is: \n" << m3 << std::endl;
	double vals[3] = {1.2,3.4,4.3};
	m2.setDiag(vals,3);
	std::cout << "M2 with arbitrary diagonal values is: \n" << m2 << std::endl;
	m2.mul(m1);
	std::cout << "M2 = M1*M2 is: \n" << m2 << std::endl;
	m2.set(2,1,0);
	std::cout << "M2 is: \n" << m2 << std::endl;
	m2.mul(m2);
	std::cout << "M2 = M2*M2 is: \n" << m2 << std::endl;
	smd::Vector v1, v2;
	v1.set(1,1,1);
	std::cout << "V = " << v1 << std::endl;
	smd::matrix_X_vector(m2,v1,v2);
	std::cout << "M2 x V = " << v2 << std::endl;
	return 0;
}
