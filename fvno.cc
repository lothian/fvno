#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

INIT_PLUGIN

using namespace boost;

namespace psi { namespace fvno {

void title(void);
void get_moinfo(boost::shared_ptr<Wavefunction> wfn, boost::shared_ptr<Chkpt> chkpt);
void integrals(void);
void denom(void);
void cleanup(void);

extern "C" 
int read_options(std::string name, Options& options)
{
  if(name == "FVNO" || options.read_globals()) {
    options.add_str("REFERENCE", "RHF");
    options.add_double("THRESHOLD", 1e-5);
  }

  return true;
}

extern "C" 
PsiReturnType fvno(Options& options)
{
  title();
  params.ref = options.get_str("REFERENCE");
  params.thresh = options.get_double("THRESHOLD");

  fprintf(outfile, "\tReference      = %s\n", params.ref.c_str());
  fprintf(outfile, "\tThreshold      = %3.1e\n", params.thresh);
  fflush(outfile);

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  get_moinfo(wfn, chkpt);
  integrals();
  denom();

  int no = moinfo.no;
  int nv = moinfo.nv;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double ****D2 = moinfo.D2;

  // MP2 energy as a test
  double emp2 = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          emp2 += ints[i][j][a+no][b+no]*L[i][j][a+no][b+no]/D2[i][j][a][b];

  fprintf(outfile, "\tMP2 correlation energy      = %20.15f\n", emp2);
  fprintf(outfile, "\tMP2 energy                  = %20.15f\n", moinfo.escf+emp2);

  //double **D = block_matrix(nv,nv);
  SharedMatrix D(new Matrix("VV One-Electron Density", nv, nv));
  double **Dp = D->pointer();
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            Dp[a][b] += (ints[m][n][b+no][e+no]/D2[m][n][b][e]) * 
                       (L[m][n][a+no][e+no]/D2[m][n][a][e]);

  D->print_out();
  SharedVector OCC(new Vector("Virtual Occupation Numbers", nv));
  SharedMatrix NO(new Matrix("Virtual NOs", nv, nv));
  D->diagonalize(NO, OCC, descending);
  NO->print_out();
  OCC->print_out();

  // Test diagonalization
  // NB: eigenvector dimensions = orig x new (i.e. ordered by cols)
  SharedMatrix TMP(new Matrix(nv, nv));
  TMP->gemm(1, 0, 1.0, NO, D, 0);
  D->gemm(0, 0, 1.0, TMP, NO, 0);
  D->print_out();

  // Determine FVNO dimension
  int nfv = 0;
  for(int i=0; i < nv; i++) {
    if(fabs(OCC->get(i)) > params.thresh) nfv++;
    else break;
  }

  // Collect the active virtual NOs
  SharedMatrix FVNO(new Matrix("FVNOs", nv, nfv));
  for(int i=0; i < nv; i++)
    for(int j=0; j < nfv; j++)
      FVNO->set(i, j, NO->get(i,j));

  FVNO->print_out();

  // Transform the VV block of the Fock matrix to the FVNO basis
  SharedMatrix Fock_vv(new Matrix("VV Fock Matrix", nv, nv));
  for(int a=0; a < nv; a++) Fock_vv->set(a, a, moinfo.fock[a+no][a+no]);
  SharedMatrix Fock_fvno(new Matrix("FVNO Fock Matrix", nfv, nfv));
  SharedMatrix TMP1(new Matrix(nfv, nv));
  TMP1->gemm(1, 0, 1, FVNO, Fock_vv, 0);
  Fock_fvno->gemm(0, 0, 1, TMP1, FVNO, 0);

  Fock_fvno->print_out();

  // Generate semicanonical FVNOs
  // Eigenvectors are (old_nfv x new_nfv)
  SharedVector F_fvno_evals(new Vector("FVNO Fock Matrix Eigenvalues", nfv));
  SharedMatrix F_fvno_evecs(new Matrix("FVNO Fock Matrix Eigenvectors", nfv, nfv));
  Fock_fvno->diagonalize(F_fvno_evecs, F_fvno_evals, ascending);

  SharedMatrix C = Process::environment.wavefunction()->Ca();
  C->print_out();

  // Extract virtual MOs
  SharedMatrix Cv(new Matrix("Virtual MOs", moinfo.nso, nv));
  for(int i=0; i < nv; i++)
    for(int j=0; j < moinfo.nso; j++)
      Cv->set(j, i, C->get(j, i+no+moinfo.nfzc));

  Cv->print_out();

  // Transform virtual MOs to semicanonical FVNO basis
  // Cv(so, mo_v) * FVNO(mo_v x nfv) * F_fvno_evecs(nfv x nfv) 
  SharedMatrix TMP2(new Matrix(moinfo.nso, nfv));
  TMP2->gemm(0, 0, 1, Cv, FVNO, 0);
  SharedMatrix Cv_new(new Matrix(moinfo.nso, nfv));
  Cv_new->gemm(0, 0, 1, TMP2, F_fvno_evecs, 0);

  // Insert the new virtual MOs and dump to both wfn and chkpt
  for(int i=0; i < nfv; i++)
    for(int j=0; j < moinfo.nso; j++)
      C->set(j, i+no+moinfo.nfzc, Cv_new->get(j, i));

  chkpt->wt_scf(C->pointer());

  // Set the number of frozen virtuals
  fprintf(outfile, "My nfzv = %d; wfn nfzv = %d\n", moinfo.nfzv, wfn->frzvpi()[0]);
  wfn->frzvpi()[0] = nfv;

  cleanup();

  return Success;
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          FVNO          *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

}} // End namespaces

