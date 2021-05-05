#include <iostream>
#include <cstdlib>


using namespace std;


main(int argc, char *argv[]){
        double r;
		if (r > 0.5) {
			cerr << " the radius must be less than or equal to 0.5" << endl;
			exit(1);
		}
        if(argc!=2){
                cerr << "Usage: cwcmaker r" << endl;
                exit(1);}
        r=atof(argv[1]);
        
        cout << "% cwc radius= " << r << endl << endl << "2 2" << endl;
        cout << "-1 0 0 -1 0 0 1" << endl;
        cout << "1 0 0 1 0" << 2*r-2 << " " << 1-2*r << endl;
        cout << "2" << endl;
        cout << "0 1 0 1 1" << endl;
        cout << "1 1 0 1 0" << endl;
        cout << "1 2" << endl;
        cout << "0 0 1 0" << endl;
}



