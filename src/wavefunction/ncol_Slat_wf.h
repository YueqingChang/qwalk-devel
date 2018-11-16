/*
 
Copyright (C) 2007 Lucas K. Wagner

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

#ifndef NCOL_SLAT_WF_H_INCLUDED

#define NCOL_SLAT_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
#include "MatrixAlgebra.h"
#include "MO_matrix.h"
#include "ncol_clark_updates.h"
class Wavefunction_data;
class ncol_Slat_wf_data;
class System;


//----------------------------------------------------------------------

template <class T> class ncol_Slat_wf_storage : public Wavefunction_storage //TODO: check this!
{
public:
  virtual ~ncol_Slat_wf_storage()
  {}
private:
  friend class ncol_Slat_wf<T>;

  //dimensions are [value gradient lap, MO]
  Array2 <T>  moVal_alpha_temp;
  Array2 <T>  moVal_beta_temp;

  // Added by Matous
  Array2 <T>  moVal_alpha_temp_2;
  Array2 <T>  moVal_beta_temp_2;
 
  Array2 < Array2 <dcomplex> > inverse_temp;
  Array2 <log_value<dcomplex> > detVal_temp;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
template <class T> class ncol_Slat_wf : public  Wavefunction
{

public:

  ncol_Slat_wf()
  {}


  virtual int nfunc() {
    return nfunc_;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);


  virtual void getLap(Wavefunction_data *, int, Wf_return &);
  virtual void evalTestPos(Array1 <doublevar> & pos, Sample_point *, Array1 <Wf_return> & wf);


  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  // Added by Matous
  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);


  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *,Templated_MO_matrix<T> * molecorb);


  void getDetVal(Wavefunction_data * wfdata,
                     Array2 < log_value <dcomplex> > detval_return)
  {
    assert(detval_return.GetDim(0)>=nfunc_);
    assert(detval_return.GetDim(1)>=ndet);

    ncol_Slat_wf_data * dataptr;
    recast(wfdata, dataptr);

    for(int f=0; f< nfunc_; f++) {
      for(int det=0; det < ndet; det++){
        detval_return(f,det) = dataptr->detwt(det)*detVal(f,det);
      }   
    }

  }

  void getDetSpinDir(Wavefunction_data * wfdata, Array3 <doublevar> spin_p, Array3 <doublevar> spin_t)
  {
    assert(spin_p.GetDim(0)>=nfunc_);
    assert(spin_p.GetDim(1)>=ndet);
    assert(spin_p.GetDim(2)>=nelectrons);
    assert(spin_t.GetDim(0)>=nfunc_);
    assert(spin_t.GetDim(1)>=ndet);
    assert(spin_t.GetDim(2)>=nelectrons);

    for(int f=0;f<nfunc_;f++){
      for(int det=0;det<ndet;det++){
        for(int e=0;e<nelectrons;e++){
          spin_p(f,det,e) =  spindir_phi(f,det,e);
          spin_t(f,det,e) =  spindir_theta(f,det,e);

        }  
      }
    } 

  }

  //--
private:


  void updateInverse(ncol_Slat_wf_data *, int e);
  int updateValNoInverse(ncol_Slat_wf_data *, int e); 
  //!< update the value, but not the inverse.  Returns 0 if the determinant is zero and updates aren't possible
  
  void calcVal(ncol_Slat_wf_data *, Sample_point *);
  void updateVal(ncol_Slat_wf_data *, Sample_point *, int);
  void calcLap(ncol_Slat_wf_data *, Sample_point *);
  void updateLap(ncol_Slat_wf_data *, Sample_point *, int);
  void getDetLap(int e, Array3<log_value <dcomplex> > & vals );
  

  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  ncol_Slat_wf_data * parent;
  Templated_MO_matrix<T> * molecorb;
  //lazy updates of the determinant(which saves a lot of time in pseudopotentials, etc)
  int inverseStale;
  int lastValUpdate;
  Array2<log_value<dcomplex> > lastDetVal;
  
  //Saved variables for electron updates
  Array3 <T>  moVal_alpha;  //5,e,i
  Array3 <T>  moVal_beta;

  Array2 <T> updatedMoVal_alpha;  //nmo,5
  Array2 <T> updatedMoVal_beta;

  Array2 < Array2 <dcomplex> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array2 <log_value<dcomplex> > detVal; //function #, determinant #


  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc_;      //!<Number of functions this class represents.

  int nelectrons;
  //don't need spin
  Array3 <doublevar> spindir_phi; // (wf,det,electron)
  Array3 <doublevar> spindir_theta;

/*
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down       //TODO check these
  Array1 <int> spin;       //!< lookup table for the spin of a given electron
  Array1 <int> rede;       //!< number of the electron within its spin channel
  Array1 <int> opspin;
*/

  Array2 <T> work1, work2; //!< Work matrices

};


/*
 
Copyright (C) 2007 Lucas K. Wagner

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

#include "Qmc_std.h"
#include "MatrixAlgebra.h"
#include "Sample_point.h"
#include "ncol_Slat_wf_data.h"

//----------------------------------------------------------------------


template <class T> 
inline void ncol_Slat_wf<T>::generateStorage(Wavefunction_storage * & wfstore)
{
  wfstore=new ncol_Slat_wf_storage<T>;
  ncol_Slat_wf_storage<T> * store;
  recast(wfstore, store);


  store->moVal_alpha_temp.Resize (5,   nmo);
  store->moVal_beta_temp.Resize (5,   nmo);
  
  // Added by Matous
  store->moVal_alpha_temp_2.Resize (5,   nmo);
  store->moVal_beta_temp_2.Resize (5,   nmo);

  store->detVal_temp.Resize(nfunc_, ndet);  //no need to keep the 2 spin channels
  store->inverse_temp.Resize(nfunc_, ndet);
  for(int i=0; i< nfunc_; i++)
  {
    for(int det=0; det < ndet; det++)
    {
      //for(int s=0; s<2; s++) //TODO: check this! no need to loop over the spin dof
      //{
        store->inverse_temp(i,det).Resize(nelectrons,nelectrons);
        store->inverse_temp(i,det)=0;
        //store->detVal_temp(i,det,s)=1;
      //}
    }
  }
}


//----------------------------------------------------------------------

template <class T> inline void ncol_Slat_wf<T>::init(Wavefunction_data * wfdata,
    Templated_MO_matrix<T> * orb)
{

  molecorb=orb;
  ncol_Slat_wf_data * dataptr;
  recast(wfdata, dataptr);
  recast(wfdata, parent);

  nfunc_=dataptr->nfunc;
  nmo=dataptr->nmo;
  //cout << "nmo=" << nmo << endl; 
  ndet=dataptr->ndet;

  //nelectrons.Resize(2);
  //nelectrons.Resize(1);
  nelectrons=dataptr->nelectrons;  // now we are not reading in the tote num from sys
  //cout << "nelectrons=" << nelectrons << endl;

  //int tote=nelectrons(0)+nelectrons(1);
  int tote=nelectrons;
  ndim=3;

  //spin.Resize(tote);
  //don't know if we need the following 2 arrays or not
  
  spindir_phi.Resize(nfunc_,ndet,tote);
  spindir_theta.Resize(nfunc_,ndet,tote);

  for(int f=0; f < nfunc_; f++){
    for(int d=0; d < ndet; d++) {
      for(int e=0; e < nelectrons; e++) {
        spindir_phi(f,d,e) = dataptr->spin_dir_phi(f,d,e);
        spindir_theta(f,d,e) = dataptr->spin_dir_theta(f,d,e);
      }
    }
  }

  //cout << "in ncol_Slat_wf:273, the spin directions are: " <<  spindir_phi(0) << ", " << spindir_theta(0)
  //<< ",... " << endl; 
 
  /*
  for(int e=nelectrons(0); e< nelectrons(0)+nelectrons(1); e++) {
    rede(e)=e-nelectrons(0);
    spin(e)=1;
    opspin(e)=0;
  }
  */
  //Properties and intermediate calculation storage.
  moVal_alpha.Resize(5,   tote, nmo); 
  moVal_beta.Resize(5,   tote, nmo); 
  //cout << "nmo=" << nmo;
  updatedMoVal_alpha.Resize(nmo,5); //these are the intermediate storage for the moVals for each e
  updatedMoVal_beta.Resize(nmo,5);

  /*
  detVal.Resize (nfunc_, ndet, 2);  // nfunc,ndet,spin 
  inverse.Resize(nfunc_, ndet, 2);
  */

  detVal.Resize (nfunc_,ndet);   // no need to keep the spin dof, just one big det.
  inverse.Resize(nfunc_,ndet);


  for(int i=0; i< nfunc_; i++) {
    for(int det=0; det < ndet; det++) {
      inverse(i,det).Resize(nelectrons, nelectrons);
      inverse(i,det)=0;
      for(int e=0; e< nelectrons; e++) {  //initialize the inverse matrice(s) to be identity
        inverse(i,det)(e,e)=1;
        inverse(i,det)(e,e)=1;
      }

      detVal(i,det)=dcomplex(1.0); // initialize detVal to be 1.0
    }
  }


  electronIsStaleVal.Resize(tote);
  electronIsStaleLap.Resize(tote);

  electronIsStaleVal=0;
  electronIsStaleLap=0;
  updateEverythingVal=1;
  updateEverythingLap=1;
  sampleAttached=0;
  dataAttached=0;
  
  inverseStale=0;
  lastValUpdate=0;
}

//----------------------------------------------------------------------

/*!
 */
template<class T> inline void ncol_Slat_wf<T>::notify(change_type change, int num)
{
  switch(change) {
    case electron_move:
      electronIsStaleVal(num)=1;
      electronIsStaleLap(num)=1;
      break;
    case all_electrons_move:
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    case wf_parm_change:
    case all_wf_parms_change:
      if(parent->optimize_mo  ) {
        updateEverythingVal=1;
        updateEverythingLap=1;
      }
      break;
    case sample_attach:
      sampleAttached=1;
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    case data_attach:
      dataAttached=1;
      updateEverythingVal=1;
      updateEverythingLap=1;
      break;
    default:
      updateEverythingVal=1;
      updateEverythingLap=1;
  }
}


//----------------------------------------------------------------------

template<class T>inline void ncol_Slat_wf<T>::saveUpdate(Sample_point * sample, int e,
                                                    Wavefunction_storage * wfstore) {
  
  ncol_Slat_wf_storage<T> * store;
  recast(wfstore, store);
  
  //presumably, if we care enough to save the update, we care enough
  //to have the inverse up to date
  if(inverseStale) {
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }

  // the spin direction of electron e  
  //doublevar phi = spindir_phi(e);
  //doublevar theta = spindir_theta(e);
  
  int ndet_save=ndet;
  if(parent->use_clark_updates) ndet_save=1;
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det<ndet_save; det++) {
      //store->inverse_temp(f,det,s)=inverse(f,det,s);
      store->inverse_temp(f,det)=inverse(f,det);
    }
    for(int det=0; det < ndet; det++) {
      //store->detVal_temp(f,det,s)=detVal(f,det,s);
      store->detVal_temp(f,det)=detVal(f,det);
    }
  }
  
  
  int norb_a=moVal_alpha.GetDim(2);
  for(int d=0; d< 5; d++) {
    for(int i=0; i< norb_a; i++) {
      store->moVal_alpha_temp(d,i)=moVal_alpha(d,e,i);
    }
  }
  int norb_b=moVal_beta.GetDim(2);
  for(int d=0; d< 5; d++) {
    for(int i=0; i< norb_b; i++) {
      store->moVal_beta_temp(d,i)=moVal_beta(d,e,i);
    }
  }
  


}

//----------------------------------------------------------------------

template<class T>inline void ncol_Slat_wf<T>::restoreUpdate(Sample_point * sample, int e,
                            Wavefunction_storage * wfstore) {

  ncol_Slat_wf_storage<T> * store;
  recast(wfstore, store);

  //doublevar phi=spindir_phi(e);
  //doublevar theta=spindir_theta(e);

  inverseStale=0;
  
  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal_alpha.GetDim(2); i++) {  //norb
      moVal_alpha(j,e,i)=store->moVal_alpha_temp(j,i);  // compute the mo value
    }
  }

  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal_beta.GetDim(2); i++) {  //norb
      moVal_beta(j,e,i)=store->moVal_beta_temp(j,i);  // compute the mo value
    }
  }


  int ndet_save=ndet;
  if(parent->use_clark_updates) ndet_save=1; //clark updates is not yet supported
  

  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet_save; det++) {
      inverse(f,det)=store->inverse_temp(f,det);
    }
    for(int det=0; det < ndet; det++) {
      detVal(f,det)=store->detVal_temp(f,det); 
    }
  }
  //It seems to be faster to update the inverse than to save it and
  //recover it.  However, it complicates the implementation too much.
  //For now, we'll disable it.
  //updateInverse(parent,e);
  
  electronIsStaleVal(e)=0;
  electronIsStaleLap(e)=0;

}

//----------------------------------------------------------------------

// Added by Matous
template <class T>inline void ncol_Slat_wf<T>::saveUpdate(Sample_point * sample, int e1, int e2,
                         Wavefunction_storage * wfstore) {

  ncol_Slat_wf_storage<T> * store;
  recast(wfstore, store);
  
  //presumably, if we care enough to save the update, we care enough
  //to have the inverse up to date
  if(inverseStale) {
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }
  
  //doublevar phi1=spindir_phi(e1), phi2=spindir_phi(e2);
  //doublevar theta1=spindir_theta(e1), theta2=spindir_theta(e2);
 
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det<ndet; det++) {
     // now the inverse matrix is just one large matrix with dim NxN
     // TODO: check this, not sure
      store->inverse_temp(f,det)=inverse(f,det);
      store->detVal_temp(f,det)=detVal(f,det);
     /*
      if ( s1 == s2 ) {  // if the two electrons have the same spins
        store->inverse_temp(f,det,s1)=inverse(f,det,s1);
        store->detVal_temp(f,det,s1)=detVal(f,det,s1);
      }
      else {  // if the two electrons have opposite spins
        store->inverse_temp(f,det,s1)=inverse(f,det,s1);
        store->inverse_temp(f,det,s2)=inverse(f,det,s2);
        store->detVal_temp(f,det,s1)=detVal(f,det,s1);
        store->detVal_temp(f,det,s2)=detVal(f,det,s2);
      }
     */

    }
  }
  
  
  for(int d=0; d< 5; d++) {
    for(int i=0; i< moVal_alpha.GetDim(2); i++) {
      store->moVal_alpha_temp(d,i)=moVal_alpha(d,e1,i);
      store->moVal_alpha_temp_2(d,i)=moVal_alpha(d,e2,i);
    }
  }
  for(int d=0; d< 5; d++) {
    for(int i=0; i< moVal_beta.GetDim(2); i++) {
      store->moVal_beta_temp(d,i)=moVal_beta(d,e1,i);
      store->moVal_beta_temp_2(d,i)=moVal_beta(d,e2,i);
    }
  }
  

}

//----------------------------------------------------------------------

// Added by Matous
template<class T> inline void ncol_Slat_wf<T>::restoreUpdate(Sample_point * sample, int e1, int e2,
                            Wavefunction_storage * wfstore)
{

  ncol_Slat_wf_storage<T> * store;
  recast(wfstore, store);
  //cout << "ncol_Slat_wf:: restoreUpdate" << endl;

  
  //doublevar phi1=spindir_phi(e1), phi2=spindir_phi(e2);
  //doublevar theta1=spindir_theta(e1), theta2=spindir_theta(e2);

  inverseStale=0;
  
  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal_alpha.GetDim(2); i++) {
      moVal_alpha(j,e1,i)=store->moVal_alpha_temp(j,i);
      moVal_alpha(j,e2,i)=store->moVal_alpha_temp_2(j,i);
    }
  }
  for(int j=0; j<5; j++) {
    for(int i=0; i<moVal_beta.GetDim(2); i++) {
      moVal_beta(j,e1,i)=store->moVal_beta_temp(j,i);
      moVal_beta(j,e2,i)=store->moVal_beta_temp_2(j,i);
    }
  }
  for(int f=0; f< nfunc_; f++) {
    for(int det=0; det < ndet; det++) {
      inverse(f,det)=store->inverse_temp(f,det);
      detVal(f,det)=store->detVal_temp(f,det);
    }
  }
  
  electronIsStaleVal(e1)=0;
  electronIsStaleLap(e1)=0;
  electronIsStaleVal(e2)=0;
  electronIsStaleLap(e2)=0;

}

//----------------------------------------------------------------------

template <class T> inline void ncol_Slat_wf<T>::updateVal(Wavefunction_data * wfdata,
                        Sample_point * sample)
{

  //updateEverythingVal=1; //TODO: for debug purpose

  assert(sampleAttached);
  assert(dataAttached);

  if(updateEverythingVal==1) {
    calcVal(parent, sample);
    updateEverythingVal=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons; e++) {
      if(electronIsStaleVal(e)) {
        updateVal(parent, sample, e);
        electronIsStaleVal(e)=0;
      }
    }
  }

}

//----------------------------------------------------------------------

template <class T> inline void ncol_Slat_wf<T>::updateLap( Wavefunction_data * wfdata,
                        Sample_point * sample)
{
  assert(sampleAttached);
  assert(dataAttached);

  
  if(updateEverythingLap==1) {
    calcLap(parent, sample);
    updateEverythingVal=0;
    updateEverythingLap=0;
    electronIsStaleLap=0;
    electronIsStaleVal=0;
  }
  else {
    for(int e=0; e< nelectrons; e++) {
      if(electronIsStaleLap(e)) {
        updateLap(parent, sample, e);   //updateLap for each electron 
        electronIsStaleLap(e)=0;
        electronIsStaleVal(e)=0;
      }
    }
    
  }

}




//-----------------------------------------------------------------------


template <> inline int ncol_Slat_wf<dcomplex>::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
  error("parmderiv not supported for complex orbitals and ncol wave functions yet");
  return 0;
}

template <> inline int ncol_Slat_wf<doublevar>::getParmDeriv(Wavefunction_data *  wfdata, 
			  Sample_point * sample ,
			  Parm_deriv_return & derivatives){
 error("parmderiv not supported for ncol wave functions yet!"); 

  /*  
  if(inverseStale) { 
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
    inverseStale=0;
  }
  
  int nparms=parent->nparms();
  int tote=nelectrons(0)+nelectrons(1);
  
  derivatives.gradient.Resize(nparms);
  derivatives.hessian.Resize(nparms, nparms);
  derivatives.gradderiv.Resize(nparms,tote,4);
  derivatives.val_gradient.Resize(tote,4);


  Wf_return lap(1,5);
  for(int e=0; e< tote; e++) {
    getLap(wfdata,e,lap);
    for(int d=1; d< 5; d++) {
      derivatives.val_gradient(e,d-1)=lap.amp(0,d);
    }
  }
  
  if(parent->optimize_mo) {
    parent->orbrot->getParmDeriv<doublevar>(parent->detwt,moVal,inverse, detVal, derivatives);
    return 1;
  }
  else if(parent->optimize_det) {
    log_value<doublevar> detsum=0;
    Array1 <log_value<doublevar> > detvals(ndet);
    Array3 <log_value <doublevar> > detgrads(ndet,tote,5);
    Array2 <log_value <doublevar> > totgrads(tote,5);
    for(int det=0; det < ndet; det++) {
      log_value<doublevar> thisdet=detVal(0,det,0)*detVal(0,det,1);
      detvals(det)=parent->detwt(det)*thisdet;
    }

    Array3 <log_value<doublevar> > tmp_detgrads;
    for(int e=0; e< tote; e++) { 
      getDetLap(e,tmp_detgrads);
      for(int det=0; det < ndet; det++) { 
        for(int d=1; d< 5; d++) {
          detgrads(det,e,d)=tmp_detgrads(0,det,d);
        }
      }
    }

    detsum=sum(detvals);
    detsum.logval*=-1;
    derivatives.gradient=0.0;
    derivatives.gradderiv=0.0;
    //---------------  set up temporary variables
    for(int e=0; e< tote; e++) {
      for(int d=1; d< 5; d++) { 

        Array1 <log_value<doublevar> > tmpgrad(ndet);
        for(int det=0; det < ndet; det++) 
          tmpgrad(det)=parent->detwt(det)*detgrads(det,e,d);
        totgrads(e,d)=sum(tmpgrad);
        totgrads(e,d)*=detsum;
      }
    }

    //---------------
    int det=parent->CSF(0).GetDim(0)-1;
    derivatives.gradderiv=0.0;
    for(int csf=1; csf < parent->ncsf; csf++) { 
      for(int j=1;j<parent->CSF(csf).GetDim(0);j++){
        doublevar coeff=parent->CSF(csf)(j);
        int index=csf-1;
        log_value<doublevar> thisdet=detVal(0,det,0)*detVal(0,det,1);
        derivatives.gradient(index)+=coeff*thisdet.val();
        for(int e=0; e< tote; e++) {
          for(int d=1; d< 5; d++) {
            derivatives.gradderiv(index,e,d-1)+=coeff
            *(detgrads(det,e,d).val()-totgrads(e,d).val()*thisdet.val())*detsum.val();
            //cout << "coeff " << coeff << " detgrad " << detgrads(det,e,d).val()
            //  << " totgrad " << totgrads(e,d).val() << " thisdet " << thisdet.val()
            //  << " detsum " << detsum.val() << endl;
            //cout << "deriv " << derivatives.gradderiv(index,e,d-1) << endl;
          }
        }
        det++;
      }
    }
    for(int csf=0; csf< nparms; csf++) {
      derivatives.gradient(csf)*=detsum.val(); 
    }
    derivatives.hessian=0;
    //for(int csf=0; csf < parent->ncsf; csf++) { 
    //  for(int e=0; e< tote; e++) { 
    //    for(int d=0; d< 4; d++) { 
    //      cout << "deriv " << derivatives.gradderiv(csf,e,d) << endl;
    //    }
    //  }
    //}
    return 1;
  }
  else { 
    derivatives.gradient=0;
    derivatives.hessian=0;
    return 1;
  }
  
  return 0;
  */

}


//------------------------------------------------------------------------


template <class T> inline void ncol_Slat_wf<T>::calcVal(ncol_Slat_wf_data * dataptr, Sample_point * sample)
{
  //Hmm, I don't completely understand why, but something is not 
  //completely clean, so we can't just cycle through and updateVal
  //This is actually probably the best way to do it anyway, since it 
  //should in theory be faster, and it gives us a clean start
  calcLap(dataptr, sample);

}

//------------------------------------------------------------------------
inline doublevar real_qw(doublevar & a) { return a; } 
inline doublevar real_qw(dcomplex & a) { return real(a); } 

template <class T>inline void ncol_Slat_wf<T>::updateInverse(ncol_Slat_wf_data * dataptr, int e) { 

  int maxmatsize=nelectrons;
  Array1 <dcomplex> modet(maxmatsize);

  //doublevar phi = spindir_phi(e);
  //doublevar theta = spindir_theta(e);

  int ndet_update=ndet;
  if(parent->use_clark_updates) ndet_update=1;
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet_update; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(real_qw(detVal(f,det).logval) < -1e200) {  
        //cout << "updateInverse::largevalue! " << endl;
        Array2 <dcomplex> allmos(nelectrons, nelectrons);
        for(int e=0; e< nelectrons; e++) {
          int curre=e;
          doublevar phi_e = spindir_phi(f,det,e);
          doublevar theta_e = spindir_theta(f,det,e);
          if(dataptr->optimize_mo){
            error("orbital rotation is not yet supported for noncollinear wave functions");
          }
          else{
            for(int i=0; i< nelectrons; i++) {
              allmos(e,i)=cos(phi_e/2.)*moVal_alpha(0,curre, dataptr->occupation_alpha(f,det)(i))
                          +sin(phi_e/2.)*dcomplex(cos(theta_e),sin(theta_e))*moVal_beta(0,curre, dataptr->occupation_beta(f,det)(i));
            }
          }
        }

#ifdef SUPERDEBUG
        cout << "Slat_wf::updateInverse: near-zero determinant " 
          << " f " << f << " det " << det << " old det " << detVal(f,det).logval << " "; 
#endif


        detVal(f,det)=
          TransposeInverseMatrix(allmos,inverse(f,det), nelectrons);
      

#ifdef SUPERDEBUG
        cout << "Slat_wf::updateInverse: near-zero determinant " 
          << " f " << f << " det " << det << " new det " << detVal(f,det).logval << " ";
#endif

        

      }
      else { 
        if(dataptr->optimize_mo){
          error("orbital rotation is not yet supported for noncollinear wave functions");
        }
        else{
          doublevar phi = spindir_phi(f,det,e); //TODO: check this
          doublevar theta = spindir_theta(f,det,e);
          for(int i = 0; i < nelectrons; i++) {
            modet(i)=cos(phi/2.)*moVal_alpha(0,e,dataptr->occupation_alpha(f,det)(i))
                     +sin(phi/2.)*dcomplex(cos(theta),sin(theta))*moVal_beta(0,e, dataptr->occupation_beta(f,det)(i));
          }
        }
        dcomplex ratio=1./InverseUpdateColumn(inverse(f,det),
            modet, e, 
            nelectrons);

        detVal(f,det)=ratio*detVal(f,det);
      }
    }
    if(parent->use_clark_updates) { 
      error("clark_updates is not yet supported for noncollinear wave functions");
      /*Array2 <T> M(nelectrons(s),updatedMoVal.GetDim(0));
      for(int i=0; i< nelectrons(s); i++){ 
        int elec=i+s*nelectrons(0);
        for(int j=0; j< updatedMoVal.GetDim(0); j++) { 
          M(i,j)=moVal(0,elec,j);
        }
      }
      Array1 <T> ratios;
      parent->excitations.clark_updates(inverse(0,0,s),M,s,ratios);
      for(int d=0; d< ndet; d++) { 
        detVal(0,d,s)=ratios(d)*detVal(0,0,s); 
      }*/

    }

  }
  
}

//------------------------------------------------------------------------

template <class T> inline int ncol_Slat_wf<T>::updateValNoInverse(ncol_Slat_wf_data * dataptr, int e) { 
  int maxmatsize=nelectrons;
  Array1 <dcomplex> modet(maxmatsize);


  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      //fill the molecular orbitals for this
      //determinant
      if(real_qw(detVal(f,det).logval) < -1e200) return 0;
    }
  }
  
  
  for(int f=0; f< nfunc_; f++)  {
    for(int det=0; det< ndet; det++)  {
      doublevar phi = spindir_phi(f,det,e);
      doublevar theta = spindir_theta(f,det,e); 
      //fill the molecular orbitals for this
      //determinant
      if(dataptr->optimize_mo){
        error("orbital rotation is not yet supported for noncollinear wave functions");
        /*Array1<T> orb;
        orb.Resize(dataptr->orbrot->Nact(det,s));
        for(int i=0;i<orb.GetDim(0);i++){
          orb(i)=moVal(0,e,dataptr->occupation(f,det,s)(i));
        }
        dataptr->orbrot->rotMoVals<T>(det,s,orb);
        for(int i=0;i<nelectrons(s);i++){
          modet(i)=orb(i);
        }*/
      }
      else{
        for(int i = 0; i < nelectrons; i++) {
          modet(i)=cos(phi/2.)*moVal_alpha(0,e,dataptr->occupation_alpha(f,det)(i))
                   +sin(phi/2.)*dcomplex(cos(theta),sin(theta))*moVal_beta(0,e, dataptr->occupation_beta(f,det)(i));
        }
      }
      
      dcomplex ratio=1./InverseUpdateColumn(inverse(f,det),
              modet, e, 
              nelectrons);
      
#ifdef SUPERDEBUG
      dcomplex tmpratio=InverseGetNewRatio(inverse(f,det,s),
                                            modet, rede(e),
                                            nelectrons(s));

      cout << "Slat_wf::updateValNoInverse: " << "ratio " << ratio 
        << " inv ratio " << tmpratio << " old detVal " << detVal(f,det,s).logval << endl;
#endif

      detVal(f,det)=ratio*detVal(f,det);
      
    }
  }
  return 1;
}

//------------------------------------------------------------------------
/*!

*/
template <class T> inline void ncol_Slat_wf<T>::updateVal( ncol_Slat_wf_data * dataptr, Sample_point * sample,int e) {
  //TODO: bug here!!!
  //cout << "for electron = " << e << ", start ncol_Slat_wf::updateVal 929" << endl;
  if(inverseStale && lastValUpdate!=e) {
    //cout << "updateVal923" << endl; 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(dataptr, lastValUpdate);
  }
  if(inverseStale && lastValUpdate==e) { 
    //cout << "updateVal929" << endl; 
    inverseStale=0;
    detVal=lastDetVal;
  }

  assert(dataptr != NULL);
  sample->updateEIDist();
  

  //update all the mo's that we will be using.
  //cout << "start updateVal for the molecular orbitals" << endl;

  molecorb->updateVal(sample,e,0,updatedMoVal_alpha); 
  molecorb->updateVal(sample,e,1,updatedMoVal_beta); 


  for(int i=0; i< updatedMoVal_alpha.GetDim(0); i++)
  {
    moVal_alpha(0,e,i)=updatedMoVal_alpha(i,0);
  }
  for(int i=0; i< updatedMoVal_beta.GetDim(0); i++)
  {
    moVal_beta(0,e,i)=updatedMoVal_beta(i,0);
  }

  inverseStale=1;
  lastValUpdate=e;
  lastDetVal=detVal;

  //TODO: for debug purposes
  updateInverse(dataptr,e);
  inverseStale=0;
  //TODO: for debug purposes

  if(!parent->use_clark_updates) { 
    if(!updateValNoInverse(dataptr, e)) { 
      inverseStale=0;
      updateInverse(dataptr,e);
    }
  }
  else {
    updateInverse(dataptr,e);
    inverseStale=0;
  } 
//  for(int d=0; d< ndet; d++) { 
//    cout << "orig " << detVal(0,d,s).val()  
//       << " update " << ratios(d)*detVal(0,0,s).val()
//       << " ratio " << detVal(0,d,s).val()/detVal(0,0,s).val() 
//       << " computed ratio " << ratios(d) << endl;
//  }
  //----------------------------------
}

//------------------------------------------------------------------------


template <class T>inline void ncol_Slat_wf<T>::getVal(Wavefunction_data * wfdata, int e,
                     Wf_return & val)
{

  //Array1 <doublevar> si(nfunc_, 0.0);
  Array2 <log_value<dcomplex> > vals(nfunc_,1,dcomplex(0.0));

  assert(val.amp.GetDim(0) >=nfunc_);
  assert(val.amp.GetDim(1) >= 1);
 
  
  ncol_Slat_wf_data * dataptr;
  recast(wfdata, dataptr);
  
  
  for(int f=0; f< nfunc_; f++) {
    Array1 <log_value<dcomplex> > detvals(ndet);
    for(int det=0; det < ndet; det++) {
      detvals(det) = dataptr->detwt(det)*detVal(f,det);  //compute the product of det's
    }
    log_value<dcomplex> totval=sum(detvals);
    //vals(f,0)=totval.logval;
    //si(f)=totval.sign;
    vals(f,0)=totval;
  }

  val.setVals(vals);

}
//----------------------------------------------------------------------

//----------------------------------------------------------------------

template <class T> inline void ncol_Slat_wf<T>::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){
  val.phase(0, 0)=0;
  val.amp(0, 0)=0;
  val.cvals(0,0)=0;
} 

//----------------------------------------------------------------------



template <class T> inline void ncol_Slat_wf<T>::calcLap(ncol_Slat_wf_data * dataptr, Sample_point * sample)
{
  inverseStale=0;
  for(int e=0; e< nelectrons; e++)  {
    //int s=spin(e);
    sample->updateEIDist();

    //update all the mo's that we will be using, using the lists made in
    //Slat_wf_data(one for each spin).
    //cout << "mo_updatelap " << endl;
   
    molecorb->updateLap(sample, e, 0, updatedMoVal_alpha); //listnum=0 for alpha orbitals
    molecorb->updateLap(sample, e, 1, updatedMoVal_beta); //listnum=1 for beta orbitals


       
    //cout << "done " << endl;
    for(int d=0; d< 5; d++)  {
      for(int i=0; i< updatedMoVal_alpha.GetDim(0); i++) {  //nmo
        moVal_alpha(d,e,i) = updatedMoVal_alpha(i,d);
        moVal_beta(d,e,i) = updatedMoVal_beta(i,d); 
      }
    }
  }

  int maxmatsize=nelectrons;
  Array2 <dcomplex> modet(maxmatsize, maxmatsize);
  //ofstream matout("matrix_out", ios::app);
  //matout.precision(15);
  //matout << "initial_matrix " << nelectrons(0) << " rows are electrons, columns are orbital values " << endl;
  for(int f=0; f< nfunc_; f++)   {
    for(int det=0; det < ndet; det++ ) {
      //for(int s=0; s< 2; s++ ) {

        for(int e=0; e< nelectrons; e++) {
          doublevar phi = spindir_phi(f,det,e);
          doublevar theta = spindir_theta(f,det,e);

          //int curre=s*nelectrons(0)+e;
          
          //----NEW----
          if(dataptr->optimize_mo){
            error("optimize_mo is not yet supported for noncollinear wave functions.");
            
          //----------   
          }
          else{
            for(int i=0; i< nelectrons; i++) {
              int curre=e; //TODO: check curre
              //cout << "occupation_alpha(" << f << "," << det << ")(" << i << ")="
              //     << dataptr->occupation_alpha(f,det)(i) << endl;
              //cout << "occupation_beta(" << f << "," << det << ")(" << i << ")=" 
              //     << dataptr->occupation_beta(f,det)(i) << endl;
 
              modet(e,i)= cos(phi/2.)*moVal_alpha(0,curre, dataptr->occupation_alpha(f,det)(i))
               +sin(phi/2.)*dcomplex(cos(theta),sin(theta))*moVal_beta(0,curre,dataptr->occupation_beta(f,det)(i));
            }
          }
        }
        
        if(nelectrons > 0) { 
          detVal(f,det)=
          TransposeInverseMatrix(modet,inverse(f,det), nelectrons);
        }
        else detVal(f,det)=dcomplex(1.0,0.);
#ifdef SUPERDEBUG
        cout << "ncol_Slat_wf::calcLap: f " << f<< " det " << det 
          << " detVal " << detVal(f,det).logval << endl;
#endif
        //if(f==0 && det==0 && s==0) matout << "determinant " << detVal(f,det,s) 
         //   << " should be " << modet(0,0)*modet(1,1)-modet(0,1)*modet(1,0) << endl;
      //} //spin loop
    }
  }
  //cout << "done " << endl;
}

//------------------------------------------------------------------------


template <class T> void ncol_Slat_wf<T>::getDetLap(int e, Array3<log_value <dcomplex> > &  vals ) { 
  vals.Resize(nfunc_,ndet,5);
  

  //int opp=opspin(e);

  //Prepare the matrices we need for the inverse.
  //There is likely a way to do this via updates, but we'll
  //leave it for now since it doesn't seem to cost too much.
  Array2 <T> & lapvec=work2;
  Array1 <T> tmplapvec(nelectrons);
  int n=moVal_alpha.GetDim(2);  // n=nelectrons = nmo
  //cout << "ncol_Slat_wf:1146 n= " << n << endl;
  if(parent->use_clark_updates) { 
    error("clark updates is not yet supported for ncol wfs");
  /*lapvec.Resize(nelectrons(s),n);
    for(int e1=0; e1< nelectrons(s); e1++) {
      int shift=e1+s*nelectrons(0);
      for(int j=0; j< n; j++) {  
        lapvec(e1,j)=moVal(0,shift,j);
      }
    }*/
  }


  for(int f=0; f< nfunc_; f++) {
    Array1 <log_value <dcomplex> > detvals(ndet);
    for(int det=0; det < ndet; det++) {
      detvals(det) = detVal(f,det);

      vals(f,det,0)=detvals(det);
    }
    
    Array1 <log_value <dcomplex> > detgrads(ndet);
    for(int i=1; i< 5; i++) {
      if(!parent->use_clark_updates) {   //Sherman-Morrison updates
        for(int det=0; det < ndet; det++) {
          doublevar phi = spindir_phi(f,det,e);
          doublevar theta = spindir_theta(f,det,e);
          dcomplex temp=0.;
          if(parent->optimize_mo){  
            error("orbital rotation is not yet supported for ncol wfs");
          }
          else{
            for(int j=0; j<nelectrons; j++) {
  
            //cout << "partial_" << i << " chi= " 
            //       << cos(phi/2.)*moVal_alpha(i , e, parent->occupation_alpha(f,det)(j) )+sin(phi/2.)*dcomplex(cos(theta),sin(theta))*moVal_beta(i,e,parent->occupation_beta(f,det)(j)) << endl;
            //cout << moVal_beta(i,e,parent->occupation_beta(f,det)(j)) << endl;
            //cout << "phi=" << phi << " theta=" << theta << " " << sin(phi/2.) << " " << dcomplex(cos(theta),sin(theta)) << endl;
              temp+=(cos(phi/2.)*moVal_alpha(i , e, parent->occupation_alpha(f,det)(j) )
                    +sin(phi/2.)*dcomplex(cos(theta),sin(theta))
                    *moVal_beta(i,e,parent->occupation_beta(f,det)(j)) )
                    *inverse(f,det)(e, j);
              //cout << "inverse(0,0)(" << e << "," << j << ")=" 
              //     << inverse(f,det)(e, j) << endl;  //TODO: the value of inverse has a bug
            }
          }

          detgrads(det)=temp; 

          //cout << "detgrads=" << temp  << endl;
          detgrads(det)*=detVal(f,det);
        }
      } //-------
      else {  //clark updates
        error("clark updates is not yet supported for ncol wfs");
      /*Array2 <T> &  tmpinverse=work1;
        tmpinverse=inverse(f,0,s);

        for(int j=0; j< n; j++) 
          lapvec(rede(e),j)=moVal(i,e,j);

        for(int j=0; j< nelectrons(s); j++) 
          tmplapvec(j)=lapvec(rede(e),parent->occupation(f,0,s)(j));

        T baseratio=1.0/InverseUpdateColumn(tmpinverse,tmplapvec,
            rede(e),nelectrons(s));
        Array1 <T> ratios;
        parent->excitations.clark_updates(tmpinverse,lapvec,s,ratios);
        detgrads(0)=baseratio*detVal(f,0,s);
        for(int d=1; d< ndet; d++) { 
          //detgrads(d)=baseratio*ratios(d)*detVal(f,0,s);
          detgrads(d)=ratios(d)*detgrads(0);
        }*/
      } //------Done clark updates
  
      //for(int d=0; d< ndet; d++) {
      //  detgrads(d)*=detVal(f,d,opp);
      //}  //TODO: check if the above is needed

      //--------------------------------
      for(int d=0; d< ndet; d++) {
        vals(f,d,i)=detgrads(d);
      }
    }

  }
}

/*!
*/

template <class T> void ncol_Slat_wf<T>::getLap(Wavefunction_data * wfdata,
                     int e, Wf_return & lap)
{
  //sum over all the determinants
  
  //Array1 <doublevar> si(nfunc_, 0.0);
  //Array2 <doublevar> vals(nfunc_,5,0.0);
  Array2 <log_value <dcomplex> > vals(nfunc_,5);
  
  Array3 <log_value<dcomplex> > detvals;
  getDetLap(e,detvals);
  Array1 <log_value<dcomplex> > tempsum(ndet);
  for(int f=0; f< nfunc_; f++) {
    for(int i=0; i< 5; i++) {
      for(int d=0;d < ndet; d++) {
        tempsum(d)=parent->detwt(d)*detvals(f,d,i);
      }
      vals(f,i)=sum(tempsum);
    }
    log_value<dcomplex> inv=vals(f,0);
    inv.logval*=-1;
    for(int i=1; i< 5; i++) vals(f,i)*=inv;
  }

  lap.setVals(vals);
  
}

//-------------------------------------------------------------------------

/*!
*/
template <class T> inline void ncol_Slat_wf<T>::updateLap(ncol_Slat_wf_data * dataptr,
                        Sample_point * sample,
                        int e ) {
  
  if(inverseStale && lastValUpdate!=e) { 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(dataptr, lastValUpdate);
  }
  if(inverseStale && lastValUpdate==e) { 
    detVal=lastDetVal;
    inverseStale=0;
  }
  assert(dataptr != NULL);

  //int s=spin(e); //no need

  sample->updateEIDist();


  //update all the mo's that we will be using.
  molecorb->updateLap(sample,e,0,updatedMoVal_alpha);
  molecorb->updateLap(sample,e,1,updatedMoVal_beta);

  for(int d=0; d< 5; d++)
    for(int i=0; i< updatedMoVal_alpha.GetDim(0); i++){
      moVal_alpha(d,e,i)=updatedMoVal_alpha(i,d);
      moVal_beta(d,e,i)=updatedMoVal_beta(i,d);
    }
  updateInverse(dataptr,e);
}

//-------------------------------------------------------------------------

template <class T> inline void ncol_Slat_wf<T>::evalTestPos(Array1 <doublevar> & pos, 
    Sample_point * sample, Array1 <Wf_return> & wf) {
  
  error("evalTestPos is not yet supported for ncol wfs");
/*
  if(inverseStale) { 
    inverseStale=0;
    detVal=lastDetVal;
    updateInverse(parent, lastValUpdate);
  }

  //int nspin=2;
  Array1<Array2 <T> > movals(nspin);
  Array1 <doublevar> oldpos(ndim);
  Array1 <T> modet(nmo);
  
  sample->getElectronPos(0,oldpos);
  for(int s=0; s< nspin; s++) {
    movals(s).Resize(nmo,1);
    sample->setElectronPosNoNotify(0,pos);
    sample->updateEIDist();
    molecorb->updateVal(sample,0,s,movals(s));
  }
  sample->setElectronPosNoNotify(0,oldpos);

  int tote=sample->electronSize();
  wf.Resize(tote);

  for(int e=0; e< tote; e++) { 
    wf(e).Resize(nfunc_,1);
    Array2 <log_value<T> > vals(nfunc_,1,T(0.0));
     
    int s=spin(e);
    int opps= s==0?1:0;
    Array1 <log_value <T> > new_detVals(ndet);
    int f=0;
    if(!parent->use_clark_updates) {  //Sherman-morrison updates
      for(int det=0; det< ndet; det++)  {
        if(parent->optimize_mo){
          Array1<T> orb;
          orb.Resize(parent->orbrot->Nact(det,s));
          for(int i=0;i<orb.GetDim(0);i++){
            orb(i)=movals(s)(parent->occupation(f,det,s)(i),0);
          }
          parent->orbrot->rotMoVals<T>(det,s,orb);
          for(int i=0;i<nelectrons(s);i++){
            modet(i)=orb(i);
          }
        }else{
          for(int i = 0; i < nelectrons(s); i++) {
            modet(i)=movals(s)(parent->occupation(f,det,s)(i),0);
          }
        }
        T ratio=1./InverseGetNewRatio(inverse(f,det,s),
            modet, rede(e),
            nelectrons(s));
        new_detVals(det)=parent->detwt(det)*detVal(f,det,s);
        new_detVals(det)*=ratio;
        new_detVals(det)*=detVal(f,det,opps);
      }
    }
    else { //Clark updates 
      Array2 <T> & motmp=work1;
      Array2 <T> & invtmp=work2;
      invtmp=inverse(f,0,s);
      int n=moVal.GetDim(2);
      motmp.Resize(nelectrons(s),n);
      for(int e1=0; e1< nelectrons(s); e1++) {
        int shift=e1+s*nelectrons(0);
        for(int j=0; j< n; j++) {  
          motmp(e1,j)=moVal(0,shift,j);
        }
      }
      for(int j=0; j< n; j++) motmp(rede(e),j)=movals(s)(j,0);
      Array1 <T> tmpvec(nelectrons(s));
      for(int j=0; j< nelectrons(s); j++) 
        tmpvec(j)=movals(s)(parent->occupation(f,0,s)(j),0);
      T baseratio=1.0/InverseUpdateColumn(invtmp,tmpvec,rede(e),nelectrons(s));
      Array1 <T> ratios;
      parent->excitations.clark_updates(invtmp,motmp,s,ratios);
      for(int d=0; d< ndet; d++) {
        new_detVals(d)=parent->detwt(d)*detVal(f,0,s);
        new_detVals(d)*=baseratio*ratios(d);
        new_detVals(d)*=detVal(f,d,opps);
      }
      
    }
    log_value<T> totval=sum(new_detVals);
    vals(f,0)=totval;
    wf(e).setVals(vals);

  }
*/
  

}


//----------------------------------------------------------------------
#endif //NCOL_SLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
