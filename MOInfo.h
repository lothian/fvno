#ifndef _psi_psi_fvno_moinfo_h
#define _psi_psi_fvno_moinfo_h

namespace psi { namespace fvno {

struct MOInfo {
  int nmo;          /* # molecular orbitals */
  int nso;          /* # symmetry orbitals */
  int nact;         /* # active orbitals */
  int no;           /* # occupied orbitals */
  int nv;           /* # unoccupied orbitals */
  int nfzc;         /* # frozen core orbitals */
  int nfzv;         /* # frozen virtual orbitals */
  double enuc;      /* nuclear repulsion energy */
  double escf;      /* SCF energy */
  double efzc;      /* frozen core energy */
  double **h;       /* h(p,q) */
  double **fock;    /* f(p,q) */
  double ****ints;  /* <pq|rs> */
  double ****L;     /* 2<pq|rs> - <pq|sr> */
  double ****D2;    /* two-electron denominators */
};

}} // namespace psi::fvno

#endif // _psi_psi_fvno_moinfo_h
