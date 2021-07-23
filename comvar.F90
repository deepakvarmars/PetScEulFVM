module comvar
#include <petsc/finclude/petsc.h>
   use petsc
   implicit none

   integer                   :: comm
   PetscInt                  :: nx, ny, nxl, nyl
   PetscReal                 :: dx, dy, xmax, xmin, ymax, ymin
   PetscInt                  :: itmax, itplot
   PetscReal                 :: gasR, gasGam
   PetscInt                  :: imin, imax, jmin, jmax
   PetscInt                  :: limin, limax, ljmin, ljmax
   PetscInt                  :: id, nprocs
   PetscInt                  :: gr = 1, gm = 2, dt = 3, cfl = 4, ft =5
   PetscInt                  :: solfile, resfile
      
end module comvar

