
#include "Average_ncol_so.h"
#include "Pseudopotential_so.h"
#include "ulec.h"
//----------------------------------------------------------------------
void Average_ncol_so::read(System * sys, Wavefunction_data * wfdata, 
                      vector <string> & words){
  unsigned int pos=0;

  vector <string> psp_file;
  psp_file.resize(0);
  vector <string> words_pspfile;
  words_pspfile.reserve(50);
  vector <string> pseudosec;
  vector < vector <string> > psp_txt;
  while(readsection(words,pos,psp_file,"READPSP_SO")!=0){
    //ifstream psp_soFile(psp_file[0]);
    ifstream psp_soFile(psp_file[0].c_str());
    parsefile(psp_soFile,words_pspfile);
  }
  pos=0;
  while( readsection(words_pspfile, pos, pseudosec, "PSEUDO") != 0) {
    psp_txt.push_back(pseudosec);
  } 
  psp_so = new Pseudopotential_so;
  psp_so->read(psp_txt,sys);

  //init_grid.Resize(2);
  //del_grid.Resize(2);  // \delta_\theta, \delta_\phi
  //num_grid.Resize(2);  // num_\theta, num_\phi
  
  int nions=sys->nIons();

  // assign the spin directions:
  // They should be able to be read off directly from the ncol_Slat_wf_data 
  //but currenly these arrays are private there, so we have to redefine them here

  vector <string> strspindir_phi;
  vector <string> strspindir_theta;
  vector <vector <string> > spindir_phi;
  vector <vector <string> > spindir_theta;
  //cout << "start reading average_ncol_so" << endl;

  pos = 0; 
  while(readsection(words, pos, strspindir_phi, "SPIN_DIR_PHI"))
  {
    spindir_phi.push_back(strspindir_phi);
  }
  //pos=startpos;

  pos = 0; 
  //cout << "spindir_phi" << spindir_phi[0][0] << endl;
  while(readsection(words, pos, strspindir_theta, "SPIN_DIR_THETA"))
  {
    spindir_theta.push_back(strspindir_theta);
  }

  //cout << "spindir_theta" << spindir_theta[0][0] << endl;
  spin_dir_phi.Resize(spindir_phi[0].size());
  spin_dir_theta.Resize(spindir_theta[0].size());
  for(int e=0; e< spindir_phi[0].size(); e++){
    spin_dir_phi(e) = atof(spindir_phi[0][e].c_str());
    spin_dir_theta(e) = atof(spindir_theta[0][e].c_str());
  }
 
  //cout << "finish reading input for Average_ncol_so" << endl; 
 
}
//----------------------------------------------------------------------
void Average_ncol_so::read(vector <string> & words){
}
//----------------------------------------------------------------------
void Average_ncol_so::write_init(string & indent, ostream & os){
  os << indent << "psp_ncol_so\n";
  os << "\n";
}
//----------------------------------------------------------------------
void Average_ncol_so::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="psp_ncol_so";

  //cout << "start evaluate" << endl;

  int nwf=wf->nfunc();
  avg.vals.Resize(nwf); // how many soi energies do we compute
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
 
  //Array1 <doublevar> totalv(nwf); //nwf=1
  //totalv = 0.0;
  
  int nrandvar=psp_so->nTest();
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++)
    rand_num(i)=rng.ulec();
 
  vector<Tmove>  tmoves; //TODO: include tmoves
  Array1 <doublevar> parm_deriv;
  Array3 <doublevar> totalv_alle(nelectrons, wf->nfunc(),3);

  psp_so->calcNonlocWithAllvariables(wfdata,sys, sample, wf, rand_num,totalv_alle, Tmoves::no_tmove, tmoves,false,parm_deriv);

  //psp_so->calcNonlocWithTest(wfdata,sys,sample,wf,rand_num,totalv);
  for (int w = 0; w < wf->nfunc(); w++)
  {
    for (int e=0; e<nelectrons; e++)
    {
      doublevar phi = spin_dir_phi(e);
      doublevar theta = spin_dir_theta(e);

      avg.vals(w) += sin(phi)*cos(theta)*totalv_alle(e, w, 0)
                   + sin(phi)*sin(theta)*totalv_alle(e, w, 1)
                   + cos(phi)*totalv_alle(e, w, 2);
    }
  }
  
}

//----------------------------------------------------------------------

void Average_ncol_so::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) {
  psp_so->randomize();
}

//----------------------------------------------------------------------

void Average_ncol_so::write_summary(Average_return & avg, Average_return & err, 
                               ostream & os) {
  int nwf = avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
//  psp_so->showinfo(os);
  os << "First order energy correction due to spin-orbit interaction\n";
 
  for(int i=0;i<nwf;i++)
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;

 /*
  for (int i=0; i<num_grid(0); i++){
    doublevar spin_theta = 0.0+i*del_grid(0);
    for (int j=0; j< num_grid(1); j++){
      doublevar spin_phi = 0.0+j*del_grid(1);
      for (int iwf=0; iwf<1; iwf++){
        doublevar spin_x = sin(spin_theta)*cos(spin_phi);
        doublevar spin_y = sin(spin_theta)*sin(spin_phi);
        doublevar spin_z = cos(spin_theta);
        os << "spin direction" << spin_theta << "  " << spin_phi <<  "   ";
        os << spin_x << "  " << spin_y << "  " << spin_z << endl; 
        os << avg.vals(iwf*1+i*num_grid(0)+j) << " +/- " << err.vals(i) << endl;
      }
    }
  }

 */

  
}

//----------------------------------------------------------------------

void Average_ncol_so::jsonOutput(Average_return & avg, Average_return & err,
                           ostream & os){
  
  os << "\"" << avg.type << "\:{" << endl;

  for(int i=0;i<avg.vals.GetDim(0);i++){
    os << "\"value\":[" << "  " <<  avg.vals(i) << " ]," << endl;
    os << "\"error\":[" << "  " << err.vals(i) << " ]," << endl; 
  }
}



