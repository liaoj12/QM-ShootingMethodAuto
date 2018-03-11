/*
  Author:         Junjie Liao
  File name:      shooting_auto.cpp
  Created date:   Feb. 22 2018
  Last modified:  Mar. 09 2018

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
//  Haven't consider the excitation energies yet, just focused on the first that
//  is found.
//*****************************************************************************

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;


const double PRECISION = 0.01;
const int    MAX_ITER  = 100;
const double INCRE     = 0.5;
const double dr        = 0.0001;
const int    nstep     = int(50/dr);
const int    pieces    = int(1/1);
const int    cpie      = 1;

// functions declarations
double V_eff(double r, double R, double C);
bool root_found(double curr, double last);
double schrodinger(double r, double u, double E, double R);
void   update(double *euler_list, double *r_list, double *du_list, double *E_list, \
              double dr, int nstep, double E, double R, double V);
void   my_euler(double *euler_list, double r, double du, double u, double dr, \
                double E, double R, double V);


// main function
int main(void)
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
  double E_max, E_min, E_mid, E_curr, E_b;
  double last, curr;
  bool found;
  int curr_iter = 0;

  // allocate memory for lists
  double *r_list      = new double[nstep+1];    // list for r, radial position
  double *du_list     = new double[nstep+1];    // list for u', 1st derivative
  double *u_list      = new double[nstep+1];    // list for u, 2nd derivative
  double *euler_list  = new double[3];          // list for r_next, u_next, u'_next
  double *R_list      = new double[pieces];     // list for R
  double *C_list      = new double[cpie];       // list for C

  // break intervals from 0? to 1 into pieces for R
  R_list[0] = 0.507;
  for (int i=1; i<pieces; i++)
  {
    R_list[i] = R_list[i-1] + INCRE;
  }  // end for-loop

  // break intervals from -1000 to 1000, with stepsize 1
  C_list[0] = 0;
  for (int i=1; i<cpie; i++)
  {
    C_list[i] = C_list[i-1] + 1;
  }

  // initial conditons
  r_list[0]  = 0;
  du_list[0] = 1;
  u_list[0]  = 0;

  //****************************************************************************
  //  Main body - loops to find everything
  //****************************************************************************
  // loop over R_list
  for (int i=0; i<pieces; i++)
  {
    R = R_list[i];

    cout << "searching for bound state energy with R=" << R << endl;
    for (int j=cpie-1; j>=0; j--)
    {
      cout << "trying searching with C=" << C << endl;
      C = C_list[j];
      found = false;

      // start-energy
      E_curr = -10;

      // update r_list, du_list, u_list
      update(euler_list, r_list, du_list, u_list, dr, nstep, E_curr, R, C);

      last = u_list[nstep];
      // "binary search"
      while (!found)
      {
        // update energy
        E_max = E_curr;

        E_curr /= 2;
        // update corresponding arrays
        update(euler_list, r_list, du_list, u_list, dr, nstep, E_curr, R, C);

        // check if root is found
        curr = u_list[nstep];
        if (root_found(curr, last))
        {
          E_min = E_curr;
          // approach in order to meet the precision limit
          while (!(fabs(E_max-E_min) <= PRECISION))
          {
            E_mid = (E_max+E_min)/2;
            update(euler_list, r_list, du_list, u_list, dr, nstep, E_mid, R, C);
            if (u_list[nstep] > 0)
            {
              E_max = E_mid;
            }
            else
            {
              E_min = E_mid;
            }
          }  // end precision_while-loop

          //********************************************************************
          // check if the root found is -0.5; if so, found;
          // if not:
          //    if E_b < -0.5, keep descending
          //    else, just skip with that value of C
          //********************************************************************
          E_b = (E_max+E_min)/2;
          cout << E_b << endl;
          if (fabs(E_b+0.5) <= PRECISION)
          {
            cout << "the bound energy is: " << E_b << endl;
            f_out << R << "\t" << C << endl;
            curr_iter = 0;
            break;
          }
          else
          {
            if (E_b < -0.5)
            {
              continue;
            }
            break;
          }
        }

        // update sign
        last = u_list[nstep];

        // check if running time is too high, if so, skips the current setup
        if (curr_iter <= MAX_ITER)
        {
          curr_iter += 1;
        }
        else
        {
          cout << "Can't find E for R=" << R << " and C=" << C << endl;
          break;
        }
      }  // end find-C_while-loop
    }  // end C-loop
  }  // end R-loop

    // close file stream
    f_out.close();

    // release memory
    delete[] r_list;
    delete[] du_list;
    delete[] u_list;
    delete[] euler_list;
    delete[] R_list;
    delete[] C_list;

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
void my_euler(double *euler_list, double r, double du, double u, double dr, \
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
            double dr, int nstep, double E, double R, double V)
{
  for (int i = 1; i <= nstep; i++)
  {
    my_euler(euler_list, r_list[i-1], du_list[i-1], u_list[i-1], dr, E, R, V);
    r_list[i] = euler_list[0];
    u_list[i] = euler_list[1];
    du_list[i] = euler_list[2];
  } // for loop

}
