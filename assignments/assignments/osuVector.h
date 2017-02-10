#ifndef OSU_VECTOR_DEFINED
#define OSU_VECTOR_DEFINED


class osuVector4 {
public:
	osuVector4();
	osuVector4(double d0, double d1, double d2, double d3);
	 ~osuVector4();

	//---------------------------------------------------
	//index operators
	double &operator[](unsigned int i);
	double operator[](unsigned int i) const;


	//----------------------------------------------------
	// in place arithmetic
	osuVector4 &operator/=(double d);
	osuVector4 &operator-(osuVector4 &a);
	osuVector4 &operator/(double d);

public:
	double data[4];
};

void normalize(osuVector4 &a);
osuVector4 cross(const osuVector4 &a, const osuVector4 &b);
double dot (const osuVector4 &a, const osuVector4 &b);



///////////////////////////////////////////////////////////
class osuVector3 {
public:
	osuVector3();
	osuVector3(double d0, double d1, double d2);
	 ~osuVector3();

	//---------------------------------------------------
	//index operators
	double &operator[](unsigned int i);
	double operator[](unsigned int i) const;


	//----------------------------------------------------
	// in place arithmetic
	osuVector3 &operator/=(double d);
	osuVector3 &operator-(osuVector3 &a);
	osuVector3 &operator/(double d);

public:
	double data[3];
};

#endif