#include <iostream>
#include <boost/random.hpp>
#include <time.h>
#include <math.h>

using namespace std;
//using namespace boost;
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

int main () {
    using namespace boost::numeric::ublas;

    vector<double> v (3);

v(1)=1;
v(2)=2;
v(3)=3;

cout <<v*v<<endl;

}
