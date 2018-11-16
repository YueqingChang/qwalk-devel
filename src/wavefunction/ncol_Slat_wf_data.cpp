#include "Qmc_std.h"
#include "qmc_io.h"
#include "ncol_Slat_wf_data.h"
#include "Wavefunction_data.h"
#include "ncol_Slat_wf.h"
#include <algorithm>
#include "MatrixAlgebra.h"
#include <map>
#include <utility>  // for make_pair

void ncol_Slat_wf_data::read(vector <string> & words, unsigned int & pos,
    System * sys)
{


  vector <string> strdetwt;
  vector <string> strstates_alpha;
  vector <string> strstates_beta;
  vector <vector <string> > statevec_alpha;
  vector <vector <string> > statevec_beta;
  vector <string> strspindir_phi;
  vector <string> strspindir_theta;
  vector <vector <string> > spindir_phi; // nfunc, ndet*nelectrons
  vector <vector <string> > spindir_theta;

  unsigned int startpos=pos;

  pos=startpos;


  while(readsection(words, pos, strstates_beta, "STATES_BETA"))
  {
    statevec_beta.push_back(strstates_beta);
  }


  pos=startpos;

  while(readsection(words, pos, strstates_alpha, "STATES_ALPHA"))
  {
    statevec_alpha.push_back(strstates_alpha);
  }


  pos=startpos;

  while(readsection(words, pos, strspindir_phi, "SPIN_DIR_PHI"))
  {
    spindir_phi.push_back(strspindir_phi);
  }

  
  pos=startpos;

  while(readsection(words, pos, strspindir_theta, "SPIN_DIR_THETA"))
  {
    spindir_theta.push_back(strspindir_theta);
  }


  //save them to arrays
  nfunc = statevec_alpha.size();
  nelectrons = sys->nelectrons(0);



  if(readvalue(words, pos=startpos, mo_place, "OPTIMIZE_MO") ) {
    error("orbital rotation is not yet supported for ncol wfs");
    //optimize_mo=1; // optimize_mo is not supported now
  }
  else optimize_mo=0;

  if(haskeyword(words, pos=startpos, "OPTIMIZE_DET")) {
    optimize_det=1;
    if(optimize_mo) error("Can't optimize both mo and det right now");
  }
  else optimize_det=0;

  sort=1;
  if(haskeyword(words, pos=startpos, "NOSORT")) {
    sort=0;
  }

  pos=startpos;
  vector <vector <string> > csfstr;
  vector <string> csfsubstr;

  while(readsection(words, pos, csfsubstr , "CSF")){
    csfstr.push_back(csfsubstr);
  }
  ncsf=csfstr.size();
  if(ncsf){
    CSF.Resize(ncsf);
    int counter=0;
    for(int csf=0;csf<ncsf;csf++){
      if(csfstr[csf].size()<2)
        error(" Wrong number of elements in the CSF number ",csf+1);
      CSF(csf).Resize(csfstr[csf].size());
      for(int j=0;j<CSF(csf).GetDim(0);j++){
        CSF(csf)(j)=atof(csfstr[csf][j].c_str());
        if(j>0)
          counter++;
      }
    }
    ndet=counter;
  

    if(readsection(words, pos, strdetwt, "DETWT")){
      error("Slater determinant: CSF is incompatible with DETWT.");
    }
    counter=0;
    detwt.Resize(ndet);
    for(int csf=0;csf<ncsf;csf++)
      for(int j=1;j<CSF(csf).GetDim(0);j++){
        detwt(counter++)=CSF(csf)(0)*CSF(csf)(j);
      }
  }
  else{
    pos=startpos;
    readsection(words, pos, strdetwt, "DETWT");
    ndet=strdetwt.size();
    detwt.Resize(ndet);
    for(int det=0; det < ndet; det++){
      // here is to convert the input detwt string to dcomplex
      string s = strdetwt[det].c_str();
      istringstream is('(' + s + ')');
      dcomplex c;
      is >> c;
      detwt(det)=c;
      cout << "detwt(" << det << ")=" << detwt(det).val() << endl;
    }
    CSF.Resize(ndet);
    ncsf=ndet;
    for(int det=0; det < ndet; det++) {
      CSF(det).Resize(2);
      CSF(det)(0)=detwt(det).val();
      CSF(det)(1)=1.0;
    }
  }
  if(fabs(CSF(0)(0).real()*CSF(0)(0).real()+CSF(0)(0).imag()*CSF(0)(0).imag()) < 1e-10)
    error("Cannot deal with the first CSF having zero weight.");


  //no sorting when ndet=1;
  if( ndet==1 && sort)
    sort=0;

  pos=startpos;
  vector <string> mowords;
  if(readsection(words,pos, mowords, "ORBITALS"))  {

    allocate(mowords, sys, molecorb);

    nmo=molecorb->getNmo();
    genmolecorb=molecorb;
    use_complexmo=0;
  }
  else if(readsection(words, pos=0, mowords, "CORBITALS")) {
    //cout << sys->nelectrons(0) << " " << sys->nelectrons(1) << endl;
    allocate(mowords, sys, cmolecorb);
    nmo=cmolecorb->getNmo();
    genmolecorb=cmolecorb;
    use_complexmo=1;
    //cout << "finish cmolecorb" << endl;
    if(optimize_mo) error("Can't optimize MO's with complex MO's");
    if(optimize_det)
      error("Don't support optimizing determinants with complex MO's yet");
  }
  else {
    error("Need ORBITALS or CORBITALS section in SLATER wave function");
  }

  //cout  << "finish reading the orb file" << endl;

  spin_dir_phi.Resize(nfunc,ndet,nelectrons);
  spin_dir_theta.Resize(nfunc,ndet,nelectrons);
  for(int f=0; f< nfunc; f++){
    for(int d=0; d< ndet; d++){
      for(int e=0; e< nelectrons; e++){
        spin_dir_phi(f,d,e) = atof(spindir_phi[f][d*nelectrons+e].c_str());
        spin_dir_theta(f,d,e) = atof(spindir_theta[f][d*nelectrons+e].c_str());
        cout << "spindir_phi(" << f << "," << d << "," << e << ")=" << spin_dir_phi(f,d,e) << " ";
        cout << "spindir_theta(" << f << "," << d << "," << e << ")=" << spin_dir_theta(f,d,e) << " ";
      }
    }
  }
  //cout << "nelectrons=" << nelectrons << endl;
  //nelectrons = sys-> nelectrons(); //Now we no longer need up or down spin flags
  pos=startpos;
  /* We don't need the NSPIN flag in sys file!
  vector <string> nspinstr;
  if(readsection(words, pos, nspinstr, "NSPIN"))
  {
    if(nspinstr.size() != 2)
      error("NSPIN must have 2 elements");
    nelectrons(0)=atoi(nspinstr[0].c_str());
    nelectrons(1)=atoi(nspinstr[1].c_str());
    if(nelectrons(0)+nelectrons(1) != sys->nelectrons(0)+sys->nelectrons(1)) {
      error("NSPIN must specify the same number of electrons as the SYSTEM "
          "in SLATER.");
    }
  }
  */

  //unsigned int canonstates=ndet*(nelectrons(0)+nelectrons(1));
  unsigned int canonstates=ndet*nelectrons*2;
  //cout << "canonstates=" << canonstates << endl;
  //cout << statevec_alpha[0].size() << endl;  TODO: somehow adding the two sizes doesn't work.
  /*
  for(int i=0; i< nfunc; i++)
  {
    if( canonstates != statevec_alpha[i].size()+statevec_beta[i].size())  //TODO: check this
    {
      error("in STATES section, expecting to find ", canonstates,
          " states(as calculated from NELEC), but found ",
          statevec_alpha[i].size()+statevec_beta[i].size(), " instead.");
    }
  }
  */

  //int tote=nelectrons(0)+nelectrons(1);
  int tote = nelectrons;
  ndim=3;
  //cout << "tote=" << tote << endl;

  //Input parameters // TODO: check if this needs to be occupation_alpha and occupation_beta
  // no longer need the spin dof, the occupation vector stores the occ mo
  occupation_alpha.Resize(nfunc, ndet); 
  occupation_beta.Resize(nfunc, ndet); 
  occupation_orig_alpha.Resize(nfunc, ndet);
  occupation_orig_beta.Resize(nfunc, ndet);



  for(int i=0; i< nfunc; i++) {
    for(int det=0; det < ndet; det++)  {
      //for(int s=0; s<2; s++) {
      //cout << det  << endl;
      occupation_alpha(i,det).Resize(nelectrons);
      occupation_beta(i,det).Resize(nelectrons);
      occupation_orig_alpha(i,det).Resize(nelectrons);
      occupation_orig_beta(i,det).Resize(nelectrons);
      //}
    }
  }


  //Calculation helpers
  for(int i=0; i< nfunc; i++) {
    //cout << "i=" << i << endl;
    int counter=0;
    for(int det=0; det<ndet; det++)  {
      //cout << "det=" << det << endl;
    //  for(int s=0; s<2; s++) {
        //cout << "s=" << s << endl;
        //cout << "nelectrons " << nelectrons(s) << endl;
      for(int e=0; e<nelectrons; e++) {
        //cout << "e=" << e << endl;
        //cout << statevec_alpha[i][counter].c_str() <<  endl;
        //cout << statevec_beta[i][counter].c_str() << endl;
        occupation_orig_alpha(i,det)(e)=atoi(statevec_alpha[i][counter].c_str())-1;
        occupation_orig_beta(i,det)(e)=atoi(statevec_beta[i][counter].c_str())-1;
        //cout << "occupation_orig_alpha(" << i << "," << det << ")(" << e << ")=" 
        //     << occupation_orig_alpha(i,det)(e) << endl;
        //cout << "occupation_orig_beta(" << i << "," << det << ")(" << e << ")=" 
        //     << occupation_orig_beta(i,det)(e) << endl;

        counter++;
      }     
    }
  }

  //Find what MO's are necessary for each spin
  //for(int s=0; s<2; s++) {
  vector <int> totocctemp_a;
  vector <int> totocctemp_b;

  for(int f=0; f< nfunc; f++)  {
    for(int det=0; det<ndet; det++) {

      //first do the alpha occupations
      for(int mo=0; mo < nelectrons; mo++) {
        int place_a=-1;
        int ntot=totocctemp_a.size();
        for(int i=0; i< ntot; i++) {
          //cout << "i " << i << endl;
          if(occupation_orig_alpha(f,det)(mo)==totocctemp_a[i]) { //TODO: check this!
            place_a=i; // 
            break;
          }
        }
        //cout << "decide " << endl;
        if(place_a==-1) { //if we didn't find the MO
          //cout << "didn't find it " << endl;
          occupation_alpha(f,det)(mo)=totocctemp_a.size();
          totocctemp_a.push_back(occupation_orig_alpha(f,det)(mo));
        }
  
        else {
          //cout << "found it" << endl;
          occupation_alpha(f,det)(mo)=place_a;
          //cout << "occupation_alpha(" << f << "," << det << ")(" << mo << ")=" 
          //<< occupation_alpha(f,det)(mo) << endl;
        }

      }//alpha loop

      //now do the beta occupations
      for(int mo=0; mo < nelectrons; mo++) {
        int place_b=-1;
        int ntot=totocctemp_b.size();
        for(int i=0; i< ntot; i++) {
          //cout << "i " << i << endl;
          if(occupation_orig_beta(f,det)(mo)==totocctemp_b[i]) { //TODO: check this!
            place_b=i; // 
            break;
          }
        }
        //cout << "decide " << endl;
        if(place_b==-1) { //if we didn't find the MO
          //cout << "didn't find it " << endl;
          occupation_beta(f,det)(mo)=totocctemp_b.size();
          totocctemp_b.push_back(occupation_orig_beta(f,det)(mo));
        }
  
        else {
          //cout << "found it" << endl;
          occupation_beta(f,det)(mo)=place_b;
          cout << "occupation_beta(" << f << "," << det << ")(" << mo << ")=" 
          << occupation_beta(f,det)(mo) << endl;
        }

      }//beta loop

    } //det loop
  }//func loop
  totoccupation_alpha.Resize(totocctemp_a.size());
  totoccupation_beta.Resize(totocctemp_b.size());
  for(int i=0; i<totoccupation_alpha.GetDim(0); i++)
    totoccupation_alpha(i) = totocctemp_a[i];

  for(int i=0; i<totoccupation_beta.GetDim(0); i++)
    totoccupation_beta(i) = totocctemp_b[i];
    //cout << "total occupation for " << s<< " : "
    // << totoccupation(s)(i) << endl;
 // } //TODO: why is this here?
//} //no need to loop over spin

  //now build an occupation list from occupation_alpha and occupation_beta
  Array3 <Array1 <int> > occupation(nfunc,ndet,2);
  for(int f=0; f<nfunc; f++){
    for(int d=0; d<ndet; d++){
      occupation(f,d,0) = occupation_alpha(f,d);
      occupation(f,d,1) = occupation_beta(f,d);
    }
  }

  excitations.build_excitation_list(occupation,0); //TODO fix the clark_updates
  //Orbital_rotation
  pos=0;
  vector <string> optstring;
  Array1<Array1<int> >activeSpace;
  Array1<Array1<int> >lists;
  activeSpace.Resize(2);
  lists.Resize(4);
  if(optimize_mo){
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    /*
    if(!readsection(words,pos,optstring,"OPTIMIZE_DATA")){
      error("Require section OPTIMIZE_DATA with option optimize_mo");
    }
    //Read and initialize Orbital_rotation
    orbrot=new Orbital_rotation;
    orbrot->read(optstring, detwt.GetDim(0), occupation_orig, occupation, totoccupation);
    nmo=orbrot->getnmo();
    */
  }

  //Decide on the updating scheme:
 
  if(ndet > 1 && tote > 10) {
    use_clark_updates=true;
  }
  else { use_clark_updates=false; }
  if(haskeyword(words,pos=startpos,"CLARK_UPDATES")) {
    use_clark_updates=true;
  }
  else if(haskeyword(words,pos=startpos,"SHERMAN_MORRISON_UPDATES")) {
    use_clark_updates=false;
  }


  //molecorb->buildLists(totoccupation);
  if(genmolecorb) init_mo();


  //cout << "finish reading ncol wf" << endl;

}

//----------------------------------------------------------------------
//TODO check this
int ncol_Slat_wf_data::supports(wf_support_type support) {
  switch(support) {
    case laplacian_update:
      return 1;
    case density:
      return 1;
    case parameter_derivatives:
      return 1;
    default:
      return 0;
  }
}

//----------------------------------------------------------------------

void ncol_Slat_wf_data::init_mo() {
  //Array1 <int> nmospin(2);
  int nmospin_alpha=0;
  int nmospin_beta=0;

  if(optimize_mo){
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //nmo=molecorb->getNmo();
  }
  for(int i=0; i< nfunc; i++)
  {
    //cout << "i=" << i << endl;
    int counter=0; // don't know what this is for.
    for(int det=0; det<ndet; det++)
    {
      //cout << "det=" << det << endl;
      //for(int s=0; s<2; s++)   // no need to loop over spins
      //{
        //cout << "s=" << s << endl;
        //cout << "nelectrons " << nelectrons(s) << endl;
      for(int e=0; e<nelectrons; e++)
      {
        if(nmospin_alpha < occupation_orig_alpha(i,det)(e)+1)
          nmospin_alpha=occupation_orig_alpha(i,det)(e)+1;

        counter++;
      }

      for(int e=0; e<nelectrons; e++)
      {
        if(nmospin_beta < occupation_orig_beta(i,det)(e)+1)
          nmospin_beta=occupation_orig_beta(i,det)(e)+1;

        counter++;
      }


      //}
    }
  }

  if(nmospin_alpha > nmo)
    error("The alpha spin channel contains an orbital higher"
        " than requested NMO's.");
  if(nmospin_beta > nmo)
    error("The beta spin channel contains an orbital higher"
        " than requested NMO's.");

  if(optimize_mo){
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //nmo=orbrot->getnmo();
  }
  

  totoccupation.Resize(2); // 0 for alpha, 1 for beta
  totoccupation(0).Resize(totoccupation_alpha.GetDim(0));
  totoccupation(1).Resize(totoccupation_beta.GetDim(0));
   
  for(int i=0;i<totoccupation_alpha.GetDim(0);i++)  
    totoccupation(0)(i) = totoccupation_alpha(i);
  for(int i=0;i<totoccupation_beta.GetDim(0);i++)  
    totoccupation(1)(i) = totoccupation_beta(i);
  



  genmolecorb->buildLists(totoccupation);

}


//----------------------------------------------------------------------

void ncol_Slat_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  if(!genmolecorb)
    error("Slat_wf_data::Need to allocate molecular orbital before generating any Wavefunction objects");

  if(use_complexmo) {
    //wf=new Cslat_wf;
    //Cslat_wf * slatwf;
    wf=new ncol_Slat_wf<dcomplex>;
    ncol_Slat_wf<dcomplex> * slatwf;
    recast(wf,slatwf);

    slatwf->init(this,cmolecorb);
    //slatwf->init(this);
    attachObserver(slatwf);
  }
 else {
    wf=new ncol_Slat_wf<doublevar>;
    ncol_Slat_wf<doublevar> * slatwf;
    recast(wf, slatwf);

  
    slatwf->init(this,molecorb);
    attachObserver(slatwf);
  }
}


//----------------------------------------------------------------------
int ncol_Slat_wf_data::showinfo(ostream & os)
{
  //cout << "ncol_Slat_wf_data::showinfo" << endl;
  if(!genmolecorb)
    error("ncol_Slat_wf_data::showinfo() : Molecular orbital not allocated");
  os << "Slater Determinant" << endl;

  if(optimize_mo)
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //os << "Optimizing molecular orbital coefficients" << endl;
  if(ndet > 1)
    os << ndet << " determinants\n";
  else
    os << "1 determinant" << endl;

  if(use_clark_updates) {
    os << "Using fast updates for multideterminants.  Reference: \n";
    os << "Clark, Morales, McMinis, Kim, and Scuseria. J. Chem. Phys. 135 244105 (2011)\n";
  }
  for(int f=0; f< nfunc; f++) {
    if(nfunc > 1)
      os << "For function " << f << endl;
    for(int det=0; det<ndet; det++) {
      if(ndet > 1) {
        os << "Determinant " << det << ":\n";
        os << "Weight: " << detwt(det).val() << endl;
      }
      os << "State: \n";
      for(int e=0; e<nelectrons;e++) {
        os << "spin directions:\n";
        os << "  ";
        os << spin_dir_phi(f,det,e) << "  " << spin_dir_theta(f,det,e) << ", "; 
        if((e+1)%10 == 0) 
          os << endl << "  ";
      }
      /*
      for(int s=0; s<2; s++) {
        if(s==0) os << "spin up:\n";
        if(s==1) os << "spin down: \n";
        os << "  ";
        for(int e=0; e<nelectrons(s); e++) {
          os << occupation_orig(f,det,s)(e)+1 << " ";
          if((e+1)%10 == 0)
            os << endl << "  ";
        }
        os << endl;
      }
      */ 
    }
  }

  os << "Molecular Orbital object : ";
  genmolecorb->showinfo(os);

  return 1;


}

//----------------------------------------------------------------------

int ncol_Slat_wf_data::writeinput(string & indent, ostream & os) {

  //cout << "start ncol_Slat_wf:: writeinput " << endl;

  if(!genmolecorb)
    error("Slat_wf_data::writeinput() : Molecular orbital not allocated");

  os << indent << "SLATER" << endl;

  if(optimize_det)
    os << indent << "OPTIMIZE_DET" << endl;
  if(use_clark_updates)
    os << indent << "CLARK_UPDATES" << endl;
  else
    os << indent << "SHERMAN_MORRISON_UPDATES" << endl;
  if(!sort)
    os << indent << "NOSORT" << endl;
  Array1 <Array1 <dcomplex> > CSF_print(ncsf);
  Array1 <dcomplex> detwt_print(ndet);
  Array2 < Array1 <int> > occupation_orig_print(nfunc,ndet);

  if(sort){
    Array1 <doublevar> csf_tmp(ncsf);
    Array1 <int> list;
    for(int csf=0;csf<ncsf;csf++)
      csf_tmp(csf)=sqrt(CSF(csf)(0).real()*CSF(csf)(0).real()+CSF(csf)(0).imag()*CSF(csf)(0).imag());
    sort_abs_values_descending(csf_tmp,csf_tmp,list);


    Array1 < Array1 <int> > det_pos(ncsf);
    int counterr=0;
    for(int csf=0;csf<ncsf;csf++){
      det_pos(csf).Resize(CSF(csf).GetDim(0)-1);
      for(int j=1;j<CSF(csf).GetDim(0);j++){
        det_pos(csf)(j-1)=counterr++;
      }
    }
    //preparing CSF_print
    int counter_new=0;
    for(int csf=0;csf<ncsf;csf++){
      CSF_print(csf).Resize(CSF(list(csf)).GetDim(0));
      for(int j=0;j<CSF(list(csf)).GetDim(0);j++){
        CSF_print(csf)(j)=CSF(list(csf))(j);
        if(j>0){
          detwt_print(counter_new++)=CSF_print(csf)(0)*CSF_print(csf)(j);
        }
      }
    }

    Array1 <int> detlist(ndet);
    //cout <<"CSF list"<<endl;
    int counter_old=0;
    for(int csf=0;csf<ncsf;csf++){
      //cout << csf<<"  ";
      for(int j=0;j<CSF(list(csf)).GetDim(0)-1;j++){
        //cout <<det_pos(list(csf))(j) <<"  ";
        detlist(counter_old++)=det_pos(list(csf))(j);
      }
      //cout <<endl;
    }

    //preparing occupation_orig_print
    for(int f=0; f< nfunc; f++)
      for(int det=0; det < ndet; det++){
        //for(int s=0; s<2; s++){  // no need to loop over spin
        occupation_orig_print(f,det).Resize(nelectrons);
        occupation_orig_print(f,det)=occupation_orig_alpha(f,detlist(det)); //TODO: fix the printing of occup
        //}
      }
  }
  else{ //no sorting 
    int counter=0;
    for(int csf=0;csf<ncsf;csf++){
      CSF_print(csf).Resize(CSF(csf).GetDim(0));
      for(int j=0;j<CSF(csf).GetDim(0);j++){
        CSF_print(csf)(j)=CSF(csf)(j);
        if(j>0)
          detwt_print(counter++)=CSF(csf)(0)*CSF(csf)(j);
      }
    }
    for(int f=0; f< nfunc; f++)
      for(int det=0; det < ndet; det++){  
       // for(int s=0; s<2; s++){   // no need to loop over spin
        occupation_orig_print(f,det).Resize(nelectrons);
        occupation_orig_print(f,det)=occupation_orig_alpha(f,det);  //TODO: same here
       // }
      }
  }
  // do printout

  for(int csf=0;csf<ncsf;csf++){
    os << indent<<" CSF { ";
    for(int j=0;j<CSF_print(csf).GetDim(0);j++){
      os <<CSF_print(csf)(j)<<"  ";
    }
    os <<"} "<<endl;
  }
  for(int f=0; f< nfunc; f++){
    os << indent << "STATES { " << endl << indent <<"  ";
    for(int det=0; det < ndet; det++){
      if(ndet>1)
        os <<"#  Determinant "<<det+1<<": weight: "<<detwt_print(det)<<endl<< indent <<"  ";
      //for(int s=0; s<2; s++){  // no need to loop over spin
      for(int e=0; e< nelectrons; e++){
        os << occupation_orig_print(f,det)(e)+1 << " ";
        if((e+1)%20 ==0)
          os << endl << indent << "  ";
      }
      os << endl << indent << "  ";
      //}
    }
    os << "}" << endl;
  }
  if(optimize_mo) {
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //orbrot->writeinput(indent,os);
  }

  if(use_complexmo) os << indent << "CORBITALS { \n";
  else os << indent << "ORBITALS {\n";
  string indent2=indent+"  ";
  genmolecorb->writeinput(indent2, os);
  os << indent << "}\n";

  return 1;
}

//------------------------------------------------------------------------
void ncol_Slat_wf_data::getVarParms(Array1 <doublevar> & parms)
{ 
  //cout <<"start getVarParms"<<endl;
  if(optimize_mo) {
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //orbrot->getParms(parms);
  }
  else if(optimize_det) {
    parms.Resize(ncsf-1); 
    for(int i=1; i< ncsf; i++) {
      //parms(i-1)=CSF(i)(0);
      parms(i-1)=sqrt(CSF(i)(0).real()*CSF(i)(0).real()+CSF(i)(0).imag()*CSF(i)(0).imag());  
      //TODO: not sure what this parms is doing if we only take the CSF amp.
    }
  }
  else {
    parms.Resize(0);
  }
}

//----------------------------------------------------------------------
void ncol_Slat_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  //cout <<"start setVarParms"<<endl;
  if(optimize_mo) {
    error("The optimize_mo option is not yet supported for noncollinear-spin wave functions.");
    //orbrot->setParms(parms);
  }
  else if(optimize_det) {
    for(int csf=1; csf< ncsf; csf++)
      CSF(csf)(0)=parms(csf-1);
    int counter=0;
    for(int csf=0; csf< ncsf; csf++) {
      for(int j=1;j<CSF(csf).GetDim(0);j++){
        detwt(counter++)=CSF(csf)(0)*CSF(csf)(j);
      }
    }
    assert(counter==ndet);
  }
  else {
    parms.Resize(0);
  } 
  
  int max=wfObserver.size();
  //cout << "slatmax " << max << endl;
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  } 
} 
//----------------------------------------------------------------------

void ncol_Slat_wf_data::linearParms(Array1 <bool> & is_linear) {
  if(optimize_det) {
    is_linear.Resize(nparms());
    is_linear=true;
  }
  else {
    is_linear.Resize(nparms());
    is_linear=false;
  }
}

//----------------------------------------------------------------------

void ncol_Slat_wf_data::renormalize(){
}

//----------------------------------------------------------------------



