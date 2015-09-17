#include <iostream>
//#include <boost/random.hpp>
#include <time.h>
#include <math.h>

using namespace std;
//using namespace boost;
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

int main () {
    
using namespace boost::numeric::ublas;

    boost::numeric::ublas::vector<double> v (3);

v(1)=1;
v(2)=2;
v(0)=3;

for(int i=0;i<1000000;i++)
{
//double b=v(1)*v(1)+v(2)*v(2)+v(0)*v(0);
double b=inner_prod (v , v);
}
/*
double v1=1;
double v2=2;
double v3=3;
for(int i=0;i<1000000;i++)
{
double b=v1*v1+v2*v2+v3*v3;
}
*/
//std::cout << 2.0 * v << std::endl;
//cout << inner_prod (v , v) << v(2);

//std::cout <<(v * 2.0)<<std::end//<<v(2)<<v(2)*v(0);
}
