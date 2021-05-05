// ********** Core version 19/1/12 **************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <complex>    

using namespace std;

#include "mb.h"


// ******************** Global variables *************************

Quad quadbank[MAXQUADS];   // quadratic functions used in geometry
Geom geoms[MAXGEOMS];      // types of segment geometry
int lens[MAXCELLS];        // number of segments in each cell
int segs[MAXCELLS][MAXSEGS];  // which geom for each segment
int collcells[MAXCELLS][MAXSEGS];  // new cell after collision at segment
Part part;  // the billiard particle
Part part1,part2; // extra copies
int cro=-1; // collrule override - control coll rule in code not in geometry

bool verbose, pmap, pflow, pabs;       // output options, set in globalinit
double tflow;
ofstream fmap, fflow, fplot;

Mat ZeroMat, IdMat, Rot, Sym;                 // matrix constants
Vec ZeroVec;

// ************** Functions for Vec and Mat ***********************

#if DIM>2

inline void Vec::set(double *vi){
	int i;
	for(i=0;i<DIM;i++){
		x[i]=vi[i];}}

inline double Vec::get(int i) const{
	return x[i];}

inline Vec& Vec::operator=(const Vec &that){
	int i;
	if(this != &that){
		for(i=0;i<DIM;i++){
			x[i]=that.x[i];}}
	return *this;}	

inline const double Vec::operator*(const Vec &that) const {
	int i;
	double ret=0.0;
	for(i=0;i<DIM;i++){
		ret+=x[i]*that.x[i];}
	return ret;}

inline const Vec Vec::operator*(const double &that) const {
	int i;
	Vec ret;
	for(i=0;i<DIM;i++){
		ret.x[i]=x[i]*that;}
	return ret;}
	
inline const Vec Vec::operator+(const Vec &that) const {
        int i;
        Vec ret;
        for(i=0;i<DIM;i++){
                ret.x[i]=x[i]+that.x[i];}
        return ret;}

inline const Vec Vec::operator-(const Vec &that) const {
        int i;
        Vec ret;
        for(i=0;i<DIM;i++){
                ret.x[i]=x[i]-that.x[i];}
        return ret;}

inline void Vec::print(ofstream & str){
	int i;
	for(i=0;i<DIM;i++){
		str << x[i] << " ";}}

inline void Vec::print(ostream & str){
        int i;
	for(i=0;i<DIM;i++){
	         str << x[i] << " ";}}

inline void Mat::set(double *mi){
	int i,j;
	for(i=0;i<DIM;i++){
		for(j=0;j<DIM;j++){
			x[i][j]=mi[j+DIM*i];}}}

inline const Vec Mat::operator*(const Vec &that) const {
	int i,j;
	Vec ret;
	for(i=0;i<DIM;i++){
		ret.x[i]=0.0;
		for(j=0;j<DIM;j++){
			ret.x[i]+=x[i][j]*that.x[j];}}
	return ret;}

inline const Mat Mat::operator*(const Mat &that) const {
        int i,j,k;
        Mat ret;
        for(i=0;i<DIM;i++){
		for(j=0;j<DIM;j++){
	                ret.x[i][j]=0.0;
	                for(k=0;k<DIM;k++){
       	                	ret.x[i][j]+=x[i][k]*that.x[k][j];}}}
        return ret;}

inline const bool Mat::operator==(const Mat &that) const {
	int i,j;
	bool b=true;
	for(i=0;i<DIM;i++){
		for(j=0;j<DIM;j++){
			if(x[i][j]!=that.x[i][j]){
				b=false;}}}
	return b;}

#else

inline void Vec::set(double *xi){
	x1=xi[0];
	x2=xi[1];}

inline void Vec::set(double x1i, double x2i){
	x1=x1i;
	x2=x2i;}

inline double Vec::get(int i) const{
	return i==0?x1:x2;}

inline Vec& Vec::operator=(const Vec &that){
        if(this != &that){
		x1=that.x1;
		x2=that.x2;}
        return *this;}

inline const double Vec::operator*(const Vec &that) const {
	return x1 * that.x1 + x2 * that.x2;}

inline const Vec Vec::operator*(const double &that) const {
	Vec ret;
	ret.x1=x1*that;
	ret.x2=x2*that;
	return ret;}

inline const Vec Vec::operator+(const Vec &that) const {
	Vec ret;
	ret.x1=x1+that.x1;
	ret.x2=x2+that.x2;
        return ret;}

inline const Vec Vec::operator-(const Vec &that) const {
        Vec ret;
        ret.x1=x1-that.x1;
        ret.x2=x2-that.x2;
        return ret;}

inline void Vec::print(ofstream &str){
	str << x1 << " " << x2 << " ";}

inline void Vec::print(ostream &str){
        str << x1 << " " << x2 << " ";}

inline void Mat::set(double *mi){
	x11=mi[0];
	x12=mi[1];
	x21=mi[2];
	x22=mi[3];}

inline void CMat::set(complex<double> *mi){
        x11=mi[0];
        x12=mi[1];
        x21=mi[2];
        x22=mi[3];}

inline void Mat::set(double x11i,double x12i,double x21i,double x22i){
	x11=x11i;
	x12=x12i;
	x21=x21i;
	x22=x22i;}

inline void CMat::set(complex<double> x11i,complex<double> x12i,complex<double> x21i,complex<double> x22i){
        x11=x11i;
        x12=x12i;
        x21=x21i;
        x22=x22i;}

inline complex<double> CMat::get(int i,int j) const {
	return i==0?(j==0?x11:x12):(j==0?x21:x22);}

inline const Mat Mat::inverse() const {
	double det;
	Mat ret;
	det=x11*x22-x12*x21;
	ret.x11=x22/det;
	ret.x12=-x12/det;
	ret.x21=-x21/det;
	ret.x22=x11/det;
	return ret;}

inline const Vec Mat::operator*(const Vec &that) const {
        Vec ret;
	ret.x1=x11*that.x1+x12*that.x2;
	ret.x2=x21*that.x1+x22*that.x2;
        return ret;}

inline const Mat Mat::operator*(const Mat &that) const {
        Mat ret;
        ret.x11=x11*that.x11+x12*that.x21;
	ret.x12=x11*that.x12+x12*that.x22;
	ret.x21=x21*that.x11+x22*that.x21;
        ret.x22=x21*that.x12+x22*that.x22;
        return ret;}

inline const CMat CMat::operator*(const CMat &that) const {
        CMat ret;
        ret.x11=x11*that.x11+x12*that.x21;
        ret.x12=x11*that.x12+x12*that.x22;
        ret.x21=x21*that.x11+x22*that.x21;
        ret.x22=x21*that.x12+x22*that.x22;
        return ret;}

inline const bool Mat::operator==(const Mat &that) const {
	return (x11==that.x11 && x12==that.x12 && x21==that.x21 && x22==that.x22);}

#endif

// ****************** Functions for Part *********************************

inline void Part::set(Vec pai,Vec qai,Vec pri,Vec qri,int celli,double ti){ 
        pa=pai;
        qa=qai;
        pr=pri;
        qr=qri;
        cell=celli;
		t=ti;
}

inline void Part::sett(double ti){
		t=ti;}

int Part::flight(){
	int seg;
	seg=flight(1.0/EPS);
	if(seg==-1){
                cerr << "Part::flight: Escaped to infinity!" << endl;
                exit(1);}
	return seg;}

int Part::flight(double tmax){
        int i,seg,walls;
        double a,b,c,dt,dtmin,t1,t2,t3,t4;
        Vec tmp,qn;
        Geom g;
        Quad q;
        seg=-1;                    // returns this if stop before collision
        dtmin=tmax-t;
        walls=0;
        if(dtmin<0.0){
                cerr << "Part::flight: Negative time" << endl;
                exit(1);}
        for(i=0;i<lens[cell];i++){
                g=geoms[segs[cell][i]];
                q=*g.bound;
                if(q.t==0){                   // flat
                        b=q.B*pr;
                        c=q.B*qr+q.C;
                        dt=-c/b;
                        if(fabs(c)<EPS){
                                walls++;
                                continue;}
                        if(dt<EPS){  // ignore negative time!
                                continue;}
                        qn=qr+pr*dt;
                        if(indom(qn,cell,i) && dt<dtmin){
                                dtmin=dt;
                                seg=i;}
                        continue;}
                tmp=q.A*qr;
                a=q.A*pr*pr;
                b=tmp*pr+q.B*pr*0.5;
                c=tmp*qr+q.B*qr+q.C;
                if(fabs(a)<EPS){      // not flat, but eq linear
                        if(fabs(c)<EPS){
                                walls++;
                                continue;}
                        dt=-0.5*c/b;
                        if(dt<EPS){
                                continue;}
                        qn=qr+pr*dt;
                        if(indom(qn,cell,i) && dt<dtmin){
                                dtmin=dt;
                                seg=i;}
                        continue;}
                t1=b*b-a*c;
                if(t1<0.0){           // eq quadratic, no solns
                        continue;}
                t1=sqrt(t1);
		t2=b>0?-b-t1:-b+t1;
                if(a*b<0){
                        dt=t2/a;
                        if(dt>EPS){
                                qn=qr+pr*dt;
                                if(indom(qn,cell,i) && dt<dtmin){
                                        dtmin=dt;
                                        seg=i;}}}
                if(fabs(c)<EPS){
                        walls++;
                        continue;}
                if(b*c<0){
                        dt=c/t2;
                        if(dt>EPS){
                                qn=qr+pr*dt;
                                if(indom(qn,cell,i) && dt<dtmin){
                                        dtmin=dt;
                                        seg=i;}}}
		if(b==0.0){           // case fails above
		        dt=a>0?t1/a:-t1/a;
		        if(dt>EPS){
		                qn=qr+pr*dt;
		                if(indom(qn,cell,i) && dt<dtmin){
	                                dtmin=dt;
	                                seg=i;}}}}
        dt=dtmin;
        if(pflow){
                t1=tflow*(floor(t/tflow)+1.0);
                while(t1<=t+dt){
                        if(pabs){
                                fflow << t1 << " ";
                                tmp=qa+pa*(t1-t);
                                tmp.print(fflow);
                                fflow << endl;}
                        else{
                                fflow << t1 << " ";
                                tmp=qr+pr*(t1-t);
                                tmp.print(fflow);
                                fflow << endl;}
                        t1+=tflow;}}
        qr=qr+pr*dt;
        qa=qa+pa*dt;
        t+=dt;
        if(walls>1){
                cerr << "Part::flight warning: ignoring corner" << endl;
				cerr << "qr=";
				qr.print(cerr);
				cerr << endl << "pr=";
				pr.print(cerr);
				cerr << endl;}
        if(verbose){
                cout << "Part::flight exiting: qr=";
                qr.print(cout);
                cout << " qa=";
                qa.print(cout);
                cout << " pr=";
                pr.print(cout);
                cout << " pa=";
                pa.print(cout);
                cout << " dt=" << dt << " t=" << t << endl;}
        return seg;}

int Part::collision(int seg){
	int cr;
	double n,cp,f2,p2,tmp;
	if(seg==-1){
		cerr << "Part:collision: Ignoring mid-air collision" << endl;
		return -1;}
	Geom g=geoms[segs[cell][seg]];
	Vec gradf;
	cr=cro>-1?cro:g.collrule;  // collrule override if not -1
	if(cr<4){
		if(cr%2){   // normal collision for cr = 1 or 3
			gradf=g.bound->evald(qr);
			pr=pr-gradf*((pr*gradf*2.0)/(gradf*gradf));
			pa=Rot*pr;}
		if(cr/2){   // translation/rotation for cr = 2 or 3
			qr=g.colq->evalq(qr);
			pr=g.colq->evalp(pr);
#if DIM==2
			g.colq->evali(Rot);
#endif
		}}
	else{
		if(cr==4){  // refraction
			gradf=g.bound->evald(qr);
			n=g.colq->C;     // refractive index
			if(n==0.0){
				cerr << "Part::collision: Zero refractive index!" << endl;
				exit(1);}
			cp=pr*gradf;
			if(cp<0){       // gradf points away from n material
				n=1.0/n;}
			f2=gradf*gradf;
			p2=pr*pr;
			tmp=n*n*(cp*cp-f2*p2)+f2*p2;
			if(tmp<0){      // total internal reflection
				pr=pr-gradf*((cp*2.0)/f2);
				pa=Rot*pr;
				cr=1;}
			else{           // note magnitude of pr varies
				pr=pr*n*n+gradf*(((cp>0?1.0:-1.0)*sqrt(tmp)-n*cp)*n/f2);
				pa=Rot*pr;}}
		else{}}              //exotica, cr=5,6,...
	if(pmap){
		fmap << t << " ";
		if(pabs){
			qa.print(fmap);}
		else{
			qr.print(fmap);}
		fmap << endl;}
	if(cr!=1){      // don't change cells for reflection
		cell=collcells[cell][seg];}
	return cr;}

void Part::print(ofstream &str){
	if(verbose){
		str << "qr= "; 
		qr.print(str);
		str << " qa= ";
		qa.print(str);
		str << " pr= ";
		pr.print(str);
		str << " pa= ";
		pa.print(str);
		str << " cell= " << cell;
		str << " t= " << t << endl;}
	else{
		qr.print(str);
		qa.print(str);
		pr.print(str);
		pa.print(str);
		str << cell << " " << t << endl;}}

void Part::print(ostream &str){
	if(verbose){
	        str << "qr= ";
		qr.print(str);
       		str << " qa= ";
		qa.print(str);
       		str << " pr= ";
		pr.print(str);
       		str << " pa= ";
		pa.print(str);
       		str << " cell= " << cell;
       		str << " t= " << t << endl;}
	else{
        	qr.print(str);
                qa.print(str);
                pr.print(str);
                pa.print(str);
                str << cell << " " << t << endl;}}

inline Vec Part::getqa() const {
	return qa;}

inline Vec Part::getpa() const {
	return pa;}

inline Vec Part::getqr() const {
	return qr;}

inline Vec Part::getpr() const {
	return pr;}

inline int Part::getcell() const {
	return cell;}

inline double Part::gett() const {
	return t;}

double Part::rcurv(const int seg) const {
	Geom g=geoms[segs[cell][seg]];
	return g.bound->rcurv(qr);}

// added function to get the collision angle ***MESS****
inline double Part::getca(int seg) const { 
		Geom g=geoms[segs[cell][seg]]; // define gradf and g so that we can use them later in the function.
		Vec gradf;
		gradf=g.bound->evald(qr);
		double cl=(pr*gradf)/(sqrt(pr*pr)*sqrt(gradf*gradf)); // this is the cosine of the modulus of the correct angle
			if (pr.get(0)*gradf.get(1) - pr.get(1)*gradf.get(0) > 0){ // use cross product to determine whether angle is positive or negative
				return sin(acos(cl)); // want to sine of the angle to pass to the output
			} else {
				return sin(acos(cl)-PI); 
		}
	}

// function to get the arclength
inline double Part::findArclength(double sradius) const { // takes argument of the small circle radius
		Vec temp = getqa(); // holds the vector of the current absolute position
		if (temp.get(0) == 0){  // checks whether x=0, then runs through cases
			if (temp.get(1) == 1){
				return PI;
			} else if (temp.get(1) == -1){
				return 2*PI;  
			} else {
				return PI*(2 + sradius); // r is the radius of the small circle ***needs to be chosen***	r=0.4 i have chosen	
			}
		} else {
			if (fabs((temp.get(0))*(temp.get(0))+(temp.get(1))*(temp.get(1))-1.0) < EPS){ // check whether it hits outer or inner circle, if it hits the outer circle..... 
				if (temp.get(0) < 0){
					return atan(temp.get(1)/temp.get(0))+3*PI/2;
				} else {
					return atan(temp.get(1)/temp.get(0))+PI/2;
				}
			} else {
				if (temp.get(0) < 0){
					return sradius*atan((-(temp.get(1)+1-sradius))/temp.get(0))+sradius*PI/2+2*PI; // ***r needs to chosen*** I have chosen r=0.4
				} else {
					return sradius*atan((-(temp.get(1)+1-sradius))/temp.get(0))+sradius*3*PI/2+2*PI;
				} 
			}
		}
	}		

// ******************** Functions for Quad *****************************

void Quad::set(Mat Ai,Vec Bi,double Ci){
	A=Ai;
	B=Bi;
	C=Ci;
	if(A==ZeroMat){
		t=0;}
	else if(A==IdMat){
		t=1;}
	else{
		t=2;}}

Mat Quad::getA() const {
	return A;}

Vec Quad::getB() const {
	return B;}

double Quad::getC() const {
	return C;}

double Quad::eval(const Vec &v) const {   //  Computes v.A.v+B.v+C
	if(t==0){
		return (B*v+C);}
	else if(t==1){
		return ((v+B)*v+C);}
	else{
		return (A*v*v+B*v+C);}}

Vec Quad::evald(const Vec &v) const {    // Computes 2A.v+B
	if(t==0){
		return B;}
	else if(t==1){
		return (v*2.0+B);}
	else{
		return (A*v*2.0+B);}}
	
Vec Quad::evalq(const Vec &v) const {    // Computes A.v+B
        if(t==0){
                return B;}
        else if(t==1){
                return (v+B);}
        else{
                return (A*v+B);}}

Vec Quad::evalp(const Vec &v) const {    // Computes A.v
        if(t==0){
                return ZeroVec;}
        else if(t==1){
                return v;}
        else{
                return A*v;}}

#if DIM==2
void Quad::evali(Mat &M) const {          // M=inverse(A).M  DIM=2 only
	if(t==1){
		return;}
	else if(t==0){
		cerr << "Quad::evali Division by zero matrix" << endl;
		exit(1);}
	M=A.inverse()*M;}
#endif

double Quad::rcurv(const Vec &x) const {   // Computes rad of curv at x
	Vec v;
	double den,v2;
	if(t==0){
		return 1.0/EPS;}           // Infinite radius of curvature
	v=evald(x);
	den=-2.0*(Sym*A*Sym*v*v);
	if(den==0.0){
		return 1.0/EPS;}
	v2=v*v;
	return sqrt(v2*v2*v2)/den;} 
	

// ********************* Miscellaneous functions ************************

inline void Geom::set(int vp,int cr,int cp,int nc,int *cnp){
	int i;
	bound=quadbank+vp;
	collrule=cr;
	colq=quadbank+cp;
	numcon=nc;
	if(numcon>MAXCON){
	        cerr << "Geom::set: numcon too big:" << numcon << endl;
	        exit(1);}
	for(i=0;i<nc;i++){
		cons[i]=quadbank+cnp[i];}}

inline bool indom(Vec q,int cell,int seg){
	Geom g=geoms[segs[cell][seg]];
	bool b=true;
	int i;
	for(i=0;i<g.numcon;i++){
		if(g.cons[i]->eval(q)<0.0){
			b=false;
			break;}}
	return b;}

void bread(ifstream &file){
	bool isflat;
	int d,bi,cr,ci,nc,nquads,ngeoms,cells,s,i,j,k,tc[MAXCON];
	double x[DIM*DIM],  C;
	Mat A;
	Vec B,tv,cons[MAXCON];
	string line;
	if(!file){
		cerr << "bread: Cannot open file" << endl;
		exit(1);}
	do{
		getline(file,line);}
	while(line[0]=='%');      // Ignore initial lines starting with %
	file >> d;                // First, check the dimension
	if (d!=DIM){
		cerr << "bread: Wrong dimension:" << d << endl;
		exit(1);}
	file >> nquads;           // Number of Quads
        if(nquads > MAXQUADS){
                cerr << "bread: Too many quadratic functions" << endl;
                exit(1);}
        if(verbose){
                cout << "Dim=" << d << " Quads=" << nquads << endl;}
	for(i=0;i<nquads;i++){      // For each quad...
                if(verbose){
                        cout << "A[" << i << "]= {";}
                for(j=0;j<DIM*DIM;j++){
                        file >> x[j];        //  DIM^2 doubles for A
                        if(verbose){
                                cout << x[j] << " ";}}
                A.set(x);
                if(verbose){
                        cout << "}, " << "B[" << i << "]= {";}
                for(j=0;j<DIM;j++){
                        file >> x[j];
                        if(verbose){
                                cout << x[j] << " ";}}   //  DIM doubles for B
                B.set(x);
                file >> C;                   // double for C
		quadbank[i].set(A,B,C);
		if(verbose){
			cout << "}, C[" << i << "]= " << C << endl;}}
	file >> ngeoms;                       // int ngeoms
	if(ngeoms > MAXGEOMS){
		cerr << "bread: Too many geometries" << endl;
		exit(1);}
	if(verbose){
		cout << " Geoms=" << ngeoms << endl;}
	for(i=0;i<ngeoms;i++){    // For each geometry...
		file >> bi >> cr >> ci >> nc;   // 4 ints
		if(verbose){
			cout << "Geom[" << i << "]: Bound= " << bi;
			cout << ", Collrule= " << cr << ", Collq= " << ci;
			cout << ", Numcon= " << nc << endl << "Con = {";}
		for(j=0;j<nc;j++){    // For each constraint...
			file >> tc[j];
			if(verbose){
				cout << tc[j] << " ";}}
		geoms[i].set(bi,cr,ci,nc,tc);
		if(verbose){
			cout << "}" << endl;}}
	file >> cells;                      // int: number of cells
	if(cells>MAXCELLS){
		cerr << "bread: Too many cells" << endl;
		exit(1);}
	if(verbose){
		cout << "ncells=" << cells;}
	for(i=0;i<cells;i++){               // for each cell...
		file >> s;               // int: number of segs
		if(s>MAXSEGS){
			cerr << "bread: Too many segments" << endl;
			exit(1);}
		lens[i]=s;
		if(verbose){
			cout << ", cell[" << i << "] nsegs=" << s << endl;}
		for(j=0;j<s;j++){
			file >> segs[i][j];        // 2 ints for each seg
			file >> collcells[i][j];
			if(verbose){
				cout << "seg " << j << ": " << segs[i][j] << " " << collcells[i][j] << endl;}}}}

#if DIM==2
void bplot(ofstream & str, int cell){   // only for D=2
	int i,j,k,f;
	double x,y,th,ap,bp,Cp,d;
	Vec q0,q1,q2,qh,Bp;
	Mat Ap;
	Geom g;
	Quad q;
	for(i=0;i<lens[cell];i++){
		g=geoms[segs[cell][i]];	
		q=*g.bound;
		q0.set(0.0,0.0);
		while((Cp=q.eval(q0))==0.0){
			x=q0.get(0);
			if(x==0.0){
				q0.set(1.0,0.0);}
			else{
				y=q0.get(1);
				q0.set(x*0.5-y*r32,y*0.5+x*r32);}}
		Ap=q.getA();
		Bp=q.evald(q0);
		for(j=0;j<PLOTPOINTS;j++){
			th=j*PI/PLOTPOINTS;
			qh.set(cos(th),sin(th));
			ap=Ap*qh*qh;
			bp=Bp*qh*0.5;
			d=bp*bp-ap*Cp;
			if(d<0.0){
				continue;}
			d=sqrt(d);
			d=bp>0?-bp-d:-bp+d;
			if(ap!=0.0){
				q1=q0+qh*(d/ap);
				f=1;
				for(k=0;k<g.numcon;k++){
					if(!indom(q1,cell,i)){
						f=0;}}
				if(f==1){
					q1.print(str);
					str << endl;}}
			if(d!=0.0){
				q1=q0+qh*(Cp/d);
				f=1;
                                for(k=0;k<g.numcon;k++){
	                                if(!indom(q1,cell,i)){
        	                                f=0;}}
                                if(f==1){
                                        q1.print(str);
                                        str << endl;}}}}}
#endif

#if DIM==2
void globalinit(){
	ZeroMat.set(0,0,0,0);
	IdMat.set(1,0,0,1);
	ZeroVec.set(0,0);
	Rot=IdMat;
	Sym.set(0,1,-1,0);
	pabs=false;
        pflow=false;
	verbose=false;
	pmap=false;
//        tflow=1000;
        }
#else
void globalinit(){
	int i,j;
	double a[DIM*DIM];
	for(i=0;i<DIM;i++){
		for(j=0;j<DIM;j++){
			a[i*DIM+j]=0;}}
	ZeroMat.set(a);
	ZeroVec.set(a);
	for(i=0;i<DIM;i++){
		a[i*DIM+i]=1;}
	IdMat.set(a);
	Rot=IdMat;
        pabs=false;
        pflow=false;
        verbose=false;
        pmap=false;}
#endif


main(int argc, char *argv[]){
	char filename[50];
	int i,j,ncoll;
	double x,y,th,r; // added r to be the small radius
	Vec q,p;
	ifstream infile;
	globalinit();
	pmap=false;           
	if(argc!=9){
		cerr << "Usage: mb infile outfile celli xi yi thi collisions smallradius" << endl;
	exit(1);}
	infile.open(argv[1]);
	bread(infile);       // reads in billiard geometry
	infile.close();
	strcpy(filename,argv[1]);  // plots a 2D billiard
	strcat(filename,".plot");
	fplot.open(filename);     
	bplot(fplot,0);
	fplot.close();
	fmap.open(argv[2]);  // replace fmap by fflow if using flow option above
	x=atof(argv[4]);     // initial x value
	y=atof(argv[5]);     // initial y value
	q.set(x,y);
	th=atof(argv[6]);    // initial direction angle
	p.set(cos(th),sin(th));
	part.set(p,q,p,q,atoi(argv[3]),0.0);
	ncoll=atoi(argv[7]);
	r=atof(argv[8]); // small radius
//	fmap << "0 " << x << " " << y << endl;
for(i=0;i<ncoll;i++){
		part.collision(j=part.flight());
		part.print(fmap);}  // output to fmap or fflow occurs automatically
	fmap.close();}


