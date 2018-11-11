 /*
 
Original Copyright (C) 2007 Lucas K. Wagner

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
*/

#ifndef NCOL_SLAT_WF_DATA_H_INCLUDED
#define NCOL_SLAT_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "ncol_clark_updates.h"
//#include "Orbital_rotation.h"
template <class T> class ncol_Slat_wf;

class ncol_Slat_wf_data : public Wavefunction_data
{
public:

  ncol_Slat_wf_data():molecorb(NULL) {}

  ~ncol_Slat_wf_data()
  {
    if(molecorb != NULL ) delete molecorb;
  }

  int optimize_mo; //!< whether to optimize MO basis
  virtual int valSize() {
    if(optimize_mo) return 0;
    if(optimize_det) return 2*detwt.GetDim(0);
    else return nfunc*2;
  }
  virtual void lockInParms(Array1<doublevar> & parms){
    if(optimize_mo)
     error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
     //orbrot->lockIn(parms);
  }
  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);
  virtual void linearParms(Array1 <bool> & is_linear);

  virtual void renormalize();

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);

  int showinfo(ostream & os);

  int writeinput(string &, ostream &);

  int nparms()
  {
    if(optimize_mo) {
      error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");}
      //return orbrot->nparms();
    else if(optimize_det){ return ncsf-1; } 
    else return 0;
  }

  ncol_Slat_wf_data * clone() const{
    return new ncol_Slat_wf_data(*this);
  }

private:
  void init_mo();
  friend class ncol_Slat_wf<doublevar>;    
  friend class ncol_Slat_wf<dcomplex>;
  friend class Cslat_wf;

  Array2 < Array1 <int> > occupation_alpha;
  Array2 < Array1 <int> > occupation_beta;
  Array2 < Array1 <int> > occupation_orig_alpha;  //(function,det)
  Array2 < Array1 <int> > occupation_orig_beta;  
  Array1 <log_value<doublevar> > detwt;
  Array1 < Array1 <int> > totoccupation; //occupations of spinors
  Array1 <int> totoccupation_alpha;
  Array1 <int> totoccupation_beta;

  int max_occupation_changes;

  int nelectrons; // maybe don't need this. We really only need the total number of electrons
  Array1 <doublevar> spin_dir_phi;
  Array1 <doublevar> spin_dir_theta;  // stores the spin directions for each electron
 
  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc;      //!<Number of separate functions this represents

  Array1 <int> orbitals_for_optimize_mo; //!< which orbitals should be optimized 
  string mo_place; //!< where to place the new mo's
  int optimize_det; //!< whether to optimize determinant coefficients
  Array1 <Array1 <doublevar> > CSF;
  int ncsf;
  int sort; //whether to sort det. weights by size


  General_MO_matrix * genmolecorb;
  MO_matrix * molecorb;
  int use_complexmo;
  bool use_clark_updates; //!<Use Bryan Clark's updates.
  ncol_Excitation_list excitations;
  Complex_MO_matrix * cmolecorb;
  //Orbital_rotation * orbrot;   //The optimize_mo option is not yet supported for noncollinear-spin wf

};

#endif //NCOL_SLAT_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
