/*
  Author:         Junjie Liao
  File name:      shooting_auto.cpp
  Created date:   Feb. 22 2018
  Last modified:  Mar. 13 2018

Comment:
  - compile:
    make
  - execute:
    ./shooting

*/

//*****************************************************************************
//  This programs implements the shooting method.
//  Note that this program hasn't considered if it implements the excitation
//  energy yet
//
//  All distance quantities are in unit of a0=0.529e-10[m]
//  All energy quantities are in unit of Eh=27.2[eV]
//
//*****************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;


const double PRECISION = 0.01;
const double REL_PRECI = 0.01;
const int    MAX_ITER  = 100;
const double INCRE     = 0.001;
const double dr        = 0.0001;
const int    nstep     = int(50/dr);
const int    pieces    = 1;

// functions declarations
bool   root_found(double curr, double last);
double V_eff(double r, double R, double C);
double schrodinger(double r, double u, double E, double R);
void   update(double *euler_list, double *r_list, double *du_list, double *E_list, \
              double E, double R, double V);
void   my_euler(double *&euler_list, double r, double du, double u, double dr, \
                double E, double R, double V);
double binary_search_E_bound(double *euler_list, double *r_list, double*du_list, \
                             double *u_list, double R, double C);


// main function
int main(int argc, char* argv[])
{
  // variables used for writing data to file
  ofstream f_out ("data.dat");
  // check if file opens successfully
  if (!f_out.is_open())
  {
    cout << "Failed to open the file, terminates!" << endl;
    return 1;
  }

  // variables used during calculation
  double R, C;
  double E_b, E_tmp;
  double last_E, last_C;

  // allocate memory for lists
  double *r_list      = new double[nstep+1];    // list for r, radial position
  double *du_list     = new double[nstep+1];    // list for u', 1st derivative
  double *u_list      = new double[nstep+1];    // list for u, 2nd derivative
  double *euler_list  = new double[3];          // list for r_next, u_next, u'_next
  double *R_list      = new double[pieces];     // list for R

  // break intervals from 0.001 to 0.507 with step 0.002
  R_list[0] = 0.507;
  for (int i=1; i<pieces; i++)
  {
    R_list[i] = (2*i+1)*INCRE;
  }  // end for-loop

  // initial conditons
  r_list[0]  = 0;
  du_list[0] = 1;
  u_list[0]  = 0;

  //****************************************************************************
  //  Main body - loops to find C given R
  //****************************************************************************
  // loop over R_list
  for (int i=pieces-1; i>=0; i--)
  {
    R = R_list[i];

    cout << "searching for bound state energy with R=" << R << endl;

    // start-C
    C = 0.053;
    last_E = 1;

    // loop to find the C for which bound state energy is approximately -0.5
    while (true)
    {
      // "binary search energy"
      E_b = binary_search_E_bound(euler_list, r_list, du_list, u_list, R, C);
      cout << "found: " << E_b << " with C: " << C << endl;
      // check if it's close to 0.5
      if (fabs((E_b+0.5)/0.5) <= REL_PRECI)
      {
        // approach C to meet the precision
        C = (C+last_C)/2;
        while (!(fabs(C-last_C) <= PRECISION))
        {
          E_tmp = binary_search_E_bound(euler_list, r_list, du_list, u_list, R, C);
          if (E_tmp > E_b)
          {
            last_C = C;
          }
        }

        cout << "The bound energy is: " << E_b << " with C=: " << C << endl;
        f_out << R << "\t" << C << "\t" << E_b << endl;
        break;
      }

      // check if C is out of range
      if (C/2 <= PRECISION)
      {
        cout << "Couldn't find such C within certin precision." << endl;
        break;
      }

      // update last if it isn't the first round
      if (last_E < 0)
      {
        if (!((E_b-last_E)<=0))
        {
          cout << "Wrong searching direction for C in descending order" << endl;
        }
        else
        {
          C /= 2;
        }
      }
      // update last to the newest E_b found
      last_E = E_b;
      last_C = C;

    } // end C-loop

  } // end R-loop

    // close file stream
    f_out.close();

    // release memory
    delete[] r_list;
    delete[] du_list;
    delete[] u_list;
    delete[] euler_list;
    delete[] R_list;

  return 0;
}  // end main()


//*****************************************************************************
//  This function determines if a root exists within a certain interval
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
//  This function calculates the effective potential at a given position r; if
//  r <= R, returns a positive energy, otherwise, returns -1/r^3.
//*****************************************************************************
double V_eff(double r, double R, double C)
{
  if (r > R) {
    return -1/pow(r, 3.0);
  }
  return C;
}

//*****************************************************************************
//  This function returns the 2nd derivative of the schrodinger equation:
//    1/2 * u'' + (E-V_eff)*u = 0
//*****************************************************************************
double schrodinger(double r, double u, double E, double R, double V)
{
  return -2*(E-V_eff(r, R, V))*u;
}

//*****************************************************************************
//  This function implements the euler method for solving the schrodinger
//  equation
//*****************************************************************************
void my_euler(double *&euler_list, double r, double du, double u, double dr, \
              double E, double R, double V)
{
  // calculates r_next, u_next and du_next in ordered
  euler_list[0] = r + dr;
  euler_list[1] = u + du * dr;
  euler_list[2] = du + dr * schrodinger(r, u, E, R, V);
}

//*****************************************************************************
//  This function updates r_list, du_list, and u_list, with a new energy E.
//*****************************************************************************
void update(double *euler_list, double *r_list, double *du_list, double *u_list, \
            double E, double R, double V)
{
  for (int i = 1; i <= nstep; i++)
  {
    my_euler(euler_list, r_list[i-1], du_list[i-1], u_list[i-1], dr, E, R, V);
    r_list[i] = euler_list[0];
    u_list[i] = euler_list[1];
    du_list[i] = euler_list[2];
  }
}

//*****************************************************************************
//  This function runs a "binary search" to find the bound state energy within
//  a certain precision; since the bound state energy must be negative,
//*****************************************************************************
double binary_search_E_bound(double *euler_list, double *r_list, double*du_list, \
                             double *u_list, double R, double C)
{
  double last, curr;
  double E_max, E_min, E_mid, E_curr = -5;
  int curr_iter = 0;

  update(euler_list, r_list, du_list, u_list, E_curr, R, C);

  last = u_list[nstep];

  while (true)
  {
    // update energy
    E_max = E_curr;

    E_curr /= 2;
    // update corresponding arrays
    update(euler_list, r_list, du_list, u_list, E_curr, R, C);

    // check if root is found
    curr = u_list[nstep];
    if (root_found(curr, last))
    {
      E_min = E_curr;
      // approach E to meet the precision
      while (!(fabs(E_max-E_min) <= PRECISION))
      {
        E_mid = (E_max+E_min)/2;
        update(euler_list, r_list, du_list, u_list, E_mid, R, C);
        if (u_list[nstep] > 0)
        {
          E_max = E_mid;
        }
        else
        {
          E_min = E_mid;
        }
      }  // end precision_while-loop
      return (E_max+E_min)/2;
    }

    // update sign
    last = u_list[nstep];

    // check if running time is too high
    if (curr_iter <= MAX_ITER)
    {
      curr_iter += 1;
      continue;
    }

    return 1.0;
  }
}
