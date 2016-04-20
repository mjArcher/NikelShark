#include <iostream>
using namespace std;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

int main ()
{
  cout << sgn(-2) << endl;

  cout << min(10, 5) << endl;

  double a=1, b=2, c=3;

  if (a < b < c)
    cout<< "this works " <<endl;

  double b1 = 0.0;
  if(b1 == 0.0)
    cout << "This is zero" << endl;

  //another way of writing the same thing
  double c1 = 0.0;
  if( (b1 == 0.0) && (c1 == 0))
    cout << "This also works" << endl;


  cout << 1/2 * 2. << endl;
  cout << 13/12 << endl;

  const double d[3] = {3./10., 3./5., 1./10.};
  cout << d[0] << endl;

  /* cout << 0/0 << endl; */

  /* cout << b1||c1 << endl; */
  /* int i; */
  /* cout << "Please enter an integer value: "; */
  /* cin >> i; */
  /* cout << "The value you entered is " << i; */
  /* cout << " and its double is " << i*2 << ".\n"; */
  /* return 0; */
}
