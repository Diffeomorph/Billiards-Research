// ************ version 19/1/12 ***************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <complex>

#define PI 3.1415926535897932385
#define r32 0.86602540378443864676

#define MAXCELLS 10 // max # of cells in billiard
#define MAXSEGS 10000  // max # of segs in a cell
#define MAXGEOMS 10000 // max # of different seg geoms
#define MAXQUADS 10000 // max # of  quadratic functions required
#define DIM 2  // dimension of the space
#define MAXCON (DIM*2)  // maximum number of domain constraints
#define EPS 1.0e-17 // Recommend 10^-17 times maximum free flight
#define double long double // Warning: platform dependent!
#define PLOTPOINTS 1000 // points for plotting conic
// Space taken by cells is MAXCELLS*MAXSEGS roughly
// Space taken by geoms is MAXGEOMS*MAXCON*DIM*DIM roughly

class Vec;
class Mat;
class CMat; // complex matrices; only defined for D=2, limited functions
class Part;
class Quad;
class Geom;

#if DIM>2
class Vec{  // vector of size DIM
friend class Mat;
public:
void set(double *);
double get(int) const;
Vec& operator=(const Vec &);
const double operator*(const Vec &) const; // dot product
const Vec operator*(const double &) const; // mult by scalar
const Vec operator+(const Vec &) const;
const Vec operator-(const Vec &) const;
void print(ofstream &);
void print(ostream &);
private:
double x[DIM];};

class Mat{   // square matrix of size DIM
public:
void set(double *);                // set from 1D array
const Vec operator*(const Vec &) const;  // computes A.v
const Mat operator*(const Mat &) const;
const bool operator==(const Mat &) const; 
// diel has: const bool Mat::operator==(const Mat &) const;
private:
double x[DIM][DIM];};

#else

class Vec{  // 2D vector
friend class Mat;
public:
	void set(double *);
	void set(double,double);  // overloaded set fn for D=2 only
	double get(int) const;          // gets ith component
	Vec& operator=(const Vec &);
	const double operator*(const Vec &) const;
	const Vec operator*(const double &) const;
	const Vec operator+(const Vec &) const;
	const Vec operator-(const Vec &) const;
	void print(ofstream &);
	void print(ostream &);
private:
	double x1,x2;};

class Mat{   // square matrix of size 2
public:
void set(double *);  
void set(double,double,double,double);  // overloaded set fn for D=2 only
const Vec operator*(const Vec &) const; 
const Mat operator*(const Mat &) const;
const bool operator==(const Mat &) const;
const Mat inverse() const;              // presently avail only for D=2
private:
double x11,x12,x21,x22;};

class CMat{   // square matrix of size 2
public:
void set(complex<double> *);
void set(complex<double>,complex<double>,complex<double>,complex<double>);  
complex<double> get(int,int) const;
const CMat operator*(const CMat &) const;
private:
complex<double> x11,x12,x21,x22;};

#endif

class Part{
public:
void set(Vec,Vec,Vec,Vec,int,double);  
void sett(double);
int flight();                   // free flight, returns collision segment
int flight(double);             // free flight with maximum time condition
int collision(int);            // collision operation, may change cell
void print(ofstream &);
void print(ostream &);
Vec getqa() const;
Vec getpa() const;
Vec getqr() const;
Vec getpr() const;
double getca(const int) const; // added to get the collision angle 
double findArclength(const double) const; // added to find the arclength
int getcell() const;
double gett() const;
double rcurv(const int) const;
// private:                     // optional privacy
Vec pa;  // a=absolute, r=relative to cell, pa not in use
Vec qa;
Vec pr;
Vec qr;
int cell;   // current cell
double t;	// time
double ca;  // collision angle
};  

class Quad{  // quadratic functions of Vec, v.A.v+B.v+C
friend int Part::flight();
friend int Part::flight(double);
public:
void set(Mat,Vec,double);
Mat getA() const;
Vec getB() const;
double getC() const;
double eval(const Vec &) const;   // v.A.v+B.v+C
Vec evald(const Vec &) const;     // 2A.v+B
Vec evalq(const Vec &) const;     // A.v+B
Vec evalp(const Vec &) const;     // A.v
void evali(Mat &) const;          // M=inverse(A).M
double rcurv(const Vec &) const;  // rad of curv 
// private:     // optional privacy
Mat A;
Vec B;
double C;
private:
int t;};    // t=0 iff A=0, t=1 iff A=I, else t=2.

class Geom{
friend int Part::flight();
friend int Part::flight(double);
friend int Part::collision(int);
friend double Part::rcurv(const int) const;
friend bool indom(Vec,int,int);
friend void bplot(ofstream &,int);
friend double Part::getca(const int) const; // to gain access to the *bound action
public:
void set(int, int, int, int, int *);
private:
Quad *bound;
int collrule;
Quad *colq;
int numcon;  // number of constraints 0..MAXCON
Quad *cons[MAXCON];};  // Domain cons[i](qr)>=0 for all i=0..numcon-1

bool indom(Vec,int,int);  // Determines whether in domain
void bread(ifstream &);     // Reads billiard data
void bplot(ofstream &,int);  // plots cell boundary
void globalinit();        // Initialise global variables
