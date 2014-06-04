#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace fvno {

void cleanup(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  int nact = moinfo.nact;

  free_block(moinfo.fock);
  free_4d_array(moinfo.ints, nact, nact, nact);
  free_4d_array(moinfo.L, nact, nact, nact);
  free_4d_array(moinfo.D2, no, no, nv);
}

}} // namespace psi::fvno
