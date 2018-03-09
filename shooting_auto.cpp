/*
  Author:         Junjie Liao
  File name:      shooting_auto.cpp
  Created date:   Feb. 22 2018
  Last modified:  Mar. 08 2018

*/

//*****************************************************************************
//  This programs implements the shooting method.
//  Note that this program hasn't considered if it implements the excitation
//  energy yet
//
//  All distance quantities are in unit of a0=0.529e-10[m]
//  All energy quantities are in unit of Eh=27.2[eV]
//
//  Don't consider the excitation energies yet, just focused on the first that
//  is found.
//*****************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;


// to do list
// - make V={+V0, r<=R; 1/r^3, r>R}, plot E(R)
// - find out under what situatoins, E(R) is a horizontal line
//   with 1/R^3 for example

const double PRECISION = 0.01;
const int    MAX_ITER  = 1000;
const double  dr = 0.0001;  // step size
// theoretically, this value is larger the better, but for computational
// purpose, it's better adjusted to just the "right" size
const int nstep = int(50/dr);
const int pieces = int(1/0.01);

// functions declarations
double V_eff(double r, double R);
bool root_found(double curr, double last);
double schrodinger(double r, double u, double E, double R);
void   update(double *euler_list, double *r_list, double *du_list, double *E_list, \
              double dr, int nstep, double E, double R);
void   my_euler(double *euler_list, double r, double du, double u, double dr, \
                double E, double R);


// main function
int main(void)
{
  // variables used for writing data to file
  ofstream f_out ("data.txt");
  // check if file opens successfully
  if (f_out.is_open())
  {
    f_out << "R\tE" <<endl;
  }
  else
  {
    cout << "Failed to open the file, terminates!" << endl;
    return 1;
  }

  // variables used during calculation
  double R;
  double E_max, E_min, E_mid, E_curr;
  double last, curr;
  bool found;
  int curr_iter = 0;

  // allocate memory for lists
  double *r_list      = new double[nstep+1];    // list for r, radial position
  double *du_list     = new double[nstep+1];    // list for u', 1st derivative
  double *u_list      = new double[nstep+1];    // list for u, 2nd derivative
  double *euler_list  = new double[3];          // list for r_next, u_next, u'_next
  double *R_list      = new double[pieces];   // list for R

  // break intervals from dr to 1 into pieces
  R_list[0] = dr;
  for (int i=1; i<pieces; i++)
  {
    R_list[i] = R_list[i-1] + dr;
  }  // end for-loop

  // initial conditons
  r_list[0]  = 0;
  du_list[0] = 1;
  u_list[0]  = 0;

  // loop over R_list
  for (int i=0; i<=pieces; i++)
  {
    R = R_list[i];

    // start-energy
    // note that it would be more preferrable to set it to the largest negative
    // double precision number in C++, but for pratical purpose, it's enough to
    // set a number as followed
    E_curr = -100000.0;

    // update r_list, du_list, u_list
    update(euler_list, r_list, du_list, u_list, dr, nstep, E_curr, R);

    cout << "searching for bound state energy with R=" << R << endl;
    // "binary search" for root
    last = u_list[nstep];
    found = false;
    while (!found)
    {
      // update energy
      E_max = E_curr;

      // if (curr_iter == 0)
      // {
      //   E_curr /= 2;
      // }
      E_curr /= 2;
      // update corresponding arrays
      update(euler_list, r_list, du_list, u_list, dr, nstep, E_curr, R);

      // check if root is found
      curr = u_list[nstep];
      if (root_found(curr, last))
      {
        E_min = E_curr;
        // approach in order to meet the precision limit
        while (!(fabs(E_max-E_min) <= PRECISION))
        {
          E_mid = (E_max+E_min)/2;
          update(euler_list, r_list, du_list, u_list, dr, nstep, E_mid, R);
          if (u_list[nstep] > 0)
          {
            E_max = E_mid;
          }
          else
          {
            E_min = E_mid;
          }
        }  // end while
        found = true;
        cout << "the bound energy is: " << (E_max+E_min)/2 << endl;
      }

      // update sign
      last = u_list[nstep];

      // check if running time is too high
      if (curr_iter <= MAX_ITER)
      {
        curr_iter += 1;
      }
      else
      {
        found = true;
        cout << "Program aborts, exceeds maximum iteration!!" << endl;
      }
    }  // end while


  }  // end for-loop


  // release memory
  delete[] r_list;
  delete[] du_list;
  delete[] u_list;
  delete[] R_list;

  return 0;
}  // end main()


//*****************************************************************************
//  This function determines if a root exists within a certain interval in a
//  more sophisticated way
//*****************************************************************************
bool root_found(double curr, double last)
{
  // most obvious case: sign changed
  if ((curr * last) < 0)
  {
    return true;
  }

  return false;
}

//*****************************************************************************
//  This function calculates the effective potential at a given
//  position r; if r <= R, returns a positive constant
//*****************************************************************************
double V_eff(double r, double R)
{
  if (r > R) {
    return -1/pow(r, 3.0);
  }

  return -1/pow(R, 3.0);
}

//*****************************************************************************
//  This function returns the 2nd derivative of the schrodinger equation:
//    1/2 * u'' + (E-V_eff)*u = 0
//*****************************************************************************
double schrodinger(double r, double u, double E, double R)
{
  return -2*(E-V_eff(r, R))*u;
}

//*****************************************************************************
//  This function implements the euler method for solving the schrodinger
//  equation
//*****************************************************************************
void my_euler(double *euler_list, double r, double du, double u, double dr, \
              double E, double R)
{

  // calculates r_next, u_next and du_next in ordered
  euler_list[0] = r + dr;
  euler_list[1] = u + du * dr;
  euler_list[2] = du + dr * schrodinger(r, u, E, R);

}

//*****************************************************************************
//  This function updates r_list, du_list, and u_list, with a new energy E.
//*****************************************************************************
void update(double *euler_list, double *r_list, double *du_list, double *u_list, \
            double dr, int nstep, double E, double R)
{
  for (int i = 1; i <= nstep; i++)
  {
    my_euler(euler_list, r_list[i-1], du_list[i-1], u_list[i-1], dr, E, R);
    r_list[i] = euler_list[0];
    u_list[i] = euler_list[1];
    du_list[i] = euler_list[2];
  } // for loop

}
