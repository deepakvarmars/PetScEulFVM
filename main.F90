!-------------------------------------
!MAIN PROGRAM - START
!-------------------------------------

program main
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscts.h>
   use petscts
   use comvar
   implicit none
   !Define all the variables
   TS  :: ts
   DM  :: da 
   Vec :: ug, rg
   PetscReal :: final_time, time1, tot_time
   PetscErrorCode :: ierr
   PetscBool ::  flg
   PetscReal :: user(5)
   !Define the external subroutines to be  used
   external :: initconds, RHSFunction, Monitor


   call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   if(ierr /= 0) stop "Unable to initialise Petsc"

   comm = PETSC_COMM_WORLD

   call MPI_COMM_RANK(comm, id, ierr); CHKERRA(ierr)
   call MPI_COMM_SIZE(comm, nprocs, ierr); CHKERRA(ierr)
   
   if(id == 0) print*, "Number of processors = ", nprocs
   
   !Read the parameters
   nx = 51
   call PetscOptionsGetInt(Petsc_Null_Options, Petsc_Null_Character, &
                           '-nx', nx, flg, ierr); CHKERRQ(ierr)
   
   ny = 51
   call PetscOptionsGetInt(Petsc_Null_Options, Petsc_Null_Character, &
                           '-ny', ny, flg, ierr); CHKERRQ(ierr)
   user(dt) = 0.000001d0
   user(cfl) = 0.2d0
   user(ft) = 150.0d0

   itmax = 1000000
   call PetscOptionsGetInt(Petsc_Null_Options, Petsc_Null_Character, &
                           '-itmax', itmax, flg, ierr); CHKERRQ(ierr)
   itplot = 1000
   call PetscOptionsGetInt(Petsc_Null_Options, Petsc_Null_Character, &
                           '-itplot', itplot, flg, ierr); CHKERRQ(ierr)

   user(gr) = 1.0d0
   user(gm) = 1.4d0
   !-------------------------------------
   !Doman boundaries
   !-------------------------------------
   xmin =  0.0d0; xmax = 10.d0
   ymin = -5.0d0; ymax = 5.d0

   dx = (xmax-xmin)/(nx-1)
   dy = (ymax-ymin)/(ny-1)

   !-------------------------------------
   !Start the timer
   !-------------------------------------
   time1 = MPI_WTime()


   !----------------------------------------------------
   !Set up the grid - specifically for isentropic vortex
   !----------------------------------------------------
   call DMDACreate2d(comm, DM_BOUNDARY_MIRROR, DM_BOUNDARY_MIRROR, &    
            DMDA_STENCIL_STAR, nx, ny, PETSC_DECIDE, PETSC_DECIDE, &
            4, 2, Petsc_Null_Integer, Petsc_Null_Integer, &
            da, ierr); CHKERRA(ierr)
   call DMSetFromOptions(da,ierr); CHKERRA(ierr)
   call DMSetUp(da,ierr); CHKERRA(ierr)

   !Create global coordinates
   call DMDASetUniformCoordinates(da, xmin, &
             xmax, ymin, ymax,&
             0.0,max(dx,dy),ierr); CHKERRA(ierr)
   
   !Create conservative variables global vector
   call DMCreateGlobalVector(da, ug, ierr); CHKERRA(ierr)

   !Generate a residual vector
   call VecDuplicate(ug, rg, ierr); CHKERRA(ierr)
  
   !Find the grid indices
   call DMDAGetCorners(da, imin, jmin, &
       PETSC_NULL_INTEGER, nxl, nyl, &
      PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)

   imax = imin + nxl - 1
   jmax = jmin + nyl - 1

   !Find the local grid indices
   call DMDAGetGhostCorners(da, limin, ljmin, &
           PETSC_NULL_INTEGER, nxl, nyl, &
           PETSC_NULL_INTEGER, ierr); CHKERRA(ierr)

   limax = limin + nxl - 1
   ljmax = ljmin + nyl - 1

   !print*, "limin, limax, ljmin, limax", limin, limax, ljmin, ljmax
   !call MPI_BARRIER(comm, ierr)
   !call PETSCFINALIZE(ierr)
   !call exit(0)

   !Set the initial conditions
   call initconds(da, ug, user, ierr)

   !----------------------------------------------------
   !Set up the time-stepping schemes 
   !----------------------------------------------------
   call TSCreate(comm, ts, ierr); CHKERRA(ierr)
   call TSSetDM(ts, da, ierr); CHKERRA(ierr)
   call TSSetProblemType(ts, TS_NONLINEAR, ierr); CHKERRA(ierr)
   call TSSetRHSFunction(ts, rg, RHSFunction, user, ierr); CHKERRA(ierr)
   call TSSetTime(ts, 0.0, ierr); CHKERRA(ierr)
   call TSSetTimeStep(ts, user(dt), ierr); CHKERRA(ierr)
   call TSSetType(ts, TSRK, ierr); CHKERRA(ierr)
   call TSRKSetType(ts, TSRK3BS, ierr); CHKERRA(ierr)
   call TSSetMaxSteps(ts, itmax, ierr); CHKERRA(ierr)
   call TSSetMaxTime(ts, user(ft), ierr); CHKERRA(ierr)
   call TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP, ierr); CHKERRA(ierr)
   call TSSetSolution(ts, ug, ierr); CHKERRA(ierr)
   call TSMonitorSet(ts, Monitor, PETSC_NULL_FUNCTION, PETSC_NULL_FUNCTION, ierr); CHKERRA(ierr)
   call TSSetFromOptions(ts, ierr); CHKERRA(ierr)
   call TSSetUp(ts, ierr); CHKERRA(ierr)
   call TSSolve(ts, ug, ierr); CHKERRA(ierr)  

   !----------------------------------------------------
   !Destroy the global data structures 
   !----------------------------------------------------
   call VecDestroy(ug, ierr)
   call VecDestroy(rg, ierr)
   call TSDestroy(ts, ierr)
   call DMDestroy(da, ierr); CHKERRA(ierr)
    
   !-------------------------------------
   !Finish the timer
   !-------------------------------------
   tot_time = MPI_WTime() - time1

   if(id == 0) then
      print*, nprocs, tot_time
   endif

   !-------------------------------------
   !Finalize
   !-------------------------------------
   call PetscFinalize(ierr)

end program main

!-------------------------------------
!Set up the initial conditions
!-------------------------------------
subroutine initconds(da, ug, user, ierr)
#include <petsc/finclude/petscdmda.h>
   use petsc
   use comvar
   implicit none
   DM         :: da
   Vec        :: ug
   PetscReal, intent(in) :: user(5)
   PetscErrorCode :: ierr
   PetscReal  :: mach, beta, aoa
   PetscReal  :: x0, y0
   PetscReal  :: xr, yr, rr, temp
   PetscInt   :: i, j, n
   PetscReal, parameter :: pi = 4.0d0*atan(1.0d0)
   PetscScalar, pointer :: arrc(:,:,:)
   PetscScalar :: primg(4)
   PetscScalar, allocatable :: cong(:,:,:)
   
   gasR = user(gr)
   gasGam = user(gm)

   mach = dsqrt(2.0d0)
   aoa = 45.0d0

   aoa = aoa*(pi/180.0d0)

   solfile = 0
   resfile = 0

   x0 = 5.0d0
   y0 = 0.0d0

   beta = 5.0d0

   allocate(cong(4, imin:imax, jmin:jmax))
   
   call DMDAVecGetArrayF90(da, ug, arrc, ierr); CHKERRA(ierr)
   cong(:,:,:) = arrc(:,:,:)

   do j = jmin, jmax
      do i = imin, imax
         xr = xmin + (i)*dx
         yr = ymin + (j)*dy
         rr = (xr - x0)**2 + (yr - y0)**2
         temp = 1.0d0 - (((gasGam - 1.0d0)*beta*beta)/(8.0d0*gasGam*pi*pi))*dexp(1.0d0 - rr)
         !Create primitive variables
         primg(1) = temp**(1.0d0/(gasGam - 1.0d0))
         primg(2) = mach*cos(aoa) - (yr - y0)*(0.5d0*beta/pi)*dexp(0.5d0*(1.0d0 - rr))
         primg(3) = mach*sin(aoa) + (xr - x0)*(0.5d0*beta/pi)*dexp(0.5d0*(1.0d0 - rr))
         primg(4) = primg(1)**gasGam
         !Create conserved variables
         cong(1,i,j) = primg(1)
         cong(2,i,j) = primg(2)*primg(1)
         cong(3,i,j) = primg(3)*primg(1)
         cong(4,i,j) = primg(4)/(gasGam - 1.0d0) + &
                     0.5d0*primg(1)*(primg(2)**2 + primg(3)**2)
      enddo
   enddo

   arrc = cong

   call DMDAVecRestoreArrayF90(da, ug, arrc, ierr); CHKERRQ(ierr) 
   deallocate(cong)

   call VecAssemblyBegin(ug, ierr); CHKERRQ(ierr)
   call VecAssemblyEnd(ug, ierr); CHKERRQ(ierr)

   return

end subroutine initconds


!-------------------------------------
!Compute Residual Function
!-------------------------------------
subroutine RHSFunction(ts, time, ug, rg, user, ierr)
#include <petsc/finclude/petscts.h>
   use petscts
   use comvar
   implicit none
   TS               :: ts
   Vec              :: ug, rg
   PetscReal, intent(in) :: user(5)
   PetscErrorCode   :: ierr
   PetscReal        :: time
   DM               :: da
   Vec              :: uloc, rloc
   PetscScalar, pointer :: arrc(:,:,:)
   PetscScalar, pointer :: arrr(:,:,:)
   PetscScalar, allocatable :: cl(:,:,:)
   PetscScalar, allocatable :: rl(:,:,:)
   PetscScalar :: pl(4), pr(4),flux(4)
   PetscInt         :: i, j, step
   PetscViewer          :: viewer
   character(len = 100) :: fname

   call TSGetDM(ts, da, ierr); CHKERRQ(ierr)
   call TSGetTime(ts, time, ierr); CHKERRQ(ierr)
   call TSGetTimeStepNumber(ts, step, ierr); CHKERRQ(ierr)

   allocate(cl(4,limin:limax,ljmin:ljmax))
   allocate(rl(4,limin:limax,ljmin:ljmax))
   
   !Create local vectors
   call DMGetLocalVector(da, uloc, ierr); CHKERRQ(ierr)
   call VecDuplicate(uloc, rloc, ierr); CHKERRQ(ierr)
   call DMGlobaltoLocalBegin(da, ug, INSERT_VALUES, uloc, ierr); CHKERRQ(ierr)
   call DMGlobaltoLocalEnd(da, ug, INSERT_VALUES, uloc, ierr); CHKERRQ(ierr)

   call DMDAVecGetArrayF90(da, uloc, arrc, ierr); CHKERRA(ierr)
   call DMDAVecGetArrayF90(da, rloc, arrr, ierr); CHKERRA(ierr)

   cl(:,:,:) = arrc(:,:,:)
   rl(:,:,:) = arrr(:,:,:)

   rl = 0.0d0

   !X-flux
   do j=ljmin+2,ljmax-2
      do i = limin+1,limax-2
         call reconstruct(cl(:,i-1,j), cl(:,i,j), &
                          cl(:,i+1,j), cl(:,i+2,j), pl, pr)
         call numflux(1.0,0.0,pl,pr,flux)
         rl(:,i,j) = rl(:,i,j) - flux(:)*dy
         rl(:,i+1,j) = rl(:,i+1,j)+ flux(:)*dy        
      enddo
   enddo

   !Y-flux
   do j=ljmin+1,ljmax-2
      do i = limin+2,limax-2
         call reconstruct(cl(:,i,j-1), cl(:,i,j), &
                          cl(:,i,j+1), cl(:,i,j+2), pl, pr)
         call numflux(0.0,1.0,pl,pr,flux)
         rl(:,i,j) = rl(:,i,j) - flux(:)*dx
         rl(:,i,j+1) = rl(:,i,j+1)+ flux(:)*dx
      enddo
   enddo

   !Set Ghost Values to zero
   rl(:,limin+1,ljmin:ljmax) = 0.0d0      !Left
   rl(:,limax-1,ljmin:ljmax) = 0.0d0      !Right
   rl(:,limin:limax,ljmin+1) = 0.0d0      !Bottom
   rl(:,limin:limax,ljmax-1) = 0.0d0      !Top

   !Restore the local residual vector
   arrr = rl
   arrc = cl

   call DMDAVecRestoreArrayF90(da, rloc, arrr, ierr); CHKERRQ(ierr)
   call DMDAVecRestoreArrayF90(da, uloc, arrc, ierr); CHKERRQ(ierr)

   !Make rg zero
   call VecSet(rg, 0.0, ierr); CHKERRQ(ierr)

   call DMLocaltoGlobalBegin(da, rloc, ADD_VALUES, rg, ierr); CHKERRQ(ierr)
   call DMLocaltoGlobalEnd(da, rloc, ADD_VALUES, rg, ierr); CHKERRQ(ierr)

   call DMRestoreLocalVector(da, rloc, ierr); CHKERRQ(ierr)
   call DMRestoreLocalVector(da, uloc, ierr); CHKERRQ(ierr)

   deallocate(cl, rl)

end subroutine  RHSFunction

!-------------------------------------
!Reconstruction subroutine
!-------------------------------------
subroutine reconstruct(conjm1, conj, conjp1, conjp2, pl, pr)
   use comvar
   implicit none

   real(8), intent(in)    :: conjm1(4), conj(4), conjp1(4), conjp2(4)
   real(8), intent(out)   :: pl(4), pr(4)
   real(8)                :: cl(4), cr(4)
   integer                :: i

   ! reconstructed states
    do i = 1, 4
       cl(i) = conj(i) + ((conj(i) - conjm1(i))/6.0d0 + &
                          (conjp1(i) - conj(i))/3.0d0)
       cr(i) = conjp1(i) - ((conjp1(i) - conj(i))/3.0d0 + &
                          (conjp2(i) - conjp1(i))/6.0d0)
    enddo

    pl(1) = cl(1)
    pr(1) = cr(1)
    
    pl(2) = cl(2)/cl(1)
    pr(2) = cr(2)/cr(1)

    pl(3) = cl(3)/cl(1)
    pr(3) = cr(3)/cr(1)

    pl(4) = (gasGam - 1.0d0)*(cl(4) - &
             0.5d0*pl(1)*(pl(2)**2 + pl(3)**2))
    pr(4) = (gasGam - 1.0d0)*(cr(4) - &
             0.5d0*pr(1)*(pr(2)**2 + pr(3)**2))

end subroutine reconstruct

!-------------------------------------
!HLLC Flux scheme
!-------------------------------------
subroutine numflux(ct, st, qinl, qinr, flux)
   use comvar
   implicit none
   real(8),intent(in)  :: qinl(4), qinr(4), ct, st
   real(8),intent(out) :: flux(4)
   ! Local variables
   real(8) :: ql(4), qr(4), rho_l_sqrt, rho_r_sqrt, fact_l, fact_r, v2_l, v2_r
   real(8) :: v_l_normal, v_r_normal
   real(8) :: p_l, p_r, h_l, h_r, c_l, c_r, e_l, e_r
   real(8) :: velocity(2), unorm(2), vel_normal, v2
   real(8) :: h, c, s_l, s_r, s_m, pStar
   real(8) :: invSLmSs, sLmuL, rhoSL, rhouSL(2), eSL
   real(8) :: invSRmSs, sRmuR, rhoSR, rhouSR(2), eSR

   unorm(1) = ct
   unorm(2) = st

   ql(1) = qinl(2)
   ql(2) = qinl(3)
   ql(3) = qinl(1)
   ql(4) = qinl(4)

   qr(1) = qinr(2)
   qr(2) = qinr(3)
   qr(3) = qinr(1)
   qr(4) = qinr(4)

   rho_l_sqrt = sqrt(ql(3))
   rho_r_sqrt = sqrt(qr(3))
   fact_l = rho_l_sqrt / (rho_l_sqrt + rho_r_sqrt)
   fact_r = 1.0 - fact_l
   
   v2_l = ql(1)**2 + ql(2)**2 
   v2_r = qr(1)**2 + qr(2)**2 
   v_l_normal = ql(1)*unorm(1) + ql(2)*unorm(2)
   v_r_normal = qr(1)*unorm(1) + qr(2)*unorm(2)

   ! Roe average velocity
   velocity = ql(1:2) * fact_l + qr(1:2) * fact_r
   vel_normal = velocity(1)*unorm(1) + velocity(2)*unorm(2)
   v2 = velocity(1)**2 + velocity(2)**2 
   
   ! pressure
   p_l = ql(4)
   p_r = qr(4)
   
   ! sound speed
   c_l = sqrt(gasGam * p_l / ql(3))
   c_r = sqrt(gasGam * p_r / qr(3))
   
   ! enthalpy
   h_l = c_l/(gasGam - 1.0) + 0.5*(ql(1)**2 + ql(2)**2)
   h_r = c_r/(gasGam - 1.0) + 0.5*(qr(1)**2 + qr(2)**2)

   ! energy per unit mass
   e_l = p_l/((gasGam-1.0)*ql(3)) + 0.5*v2_l
   e_r = p_r/((gasGam-1.0)*qr(3)) + 0.5*v2_r
   
   ! roe average enthalpy and sound speed
   h = h_l * fact_l + h_r * fact_r
   c = sqrt( (gasGam-1.0) * (h - 0.5*v2) )
   
   ! speed of sound at l and r
   s_l = min(vel_normal-c, v_l_normal-c_l)
   s_r = max(vel_normal+c, v_r_normal+c_r)

   ! speed of contact
   s_m = (p_l - p_r &
          - ql(3) * v_l_normal * (s_l-v_l_normal) &
          + qr(3) * v_r_normal * (s_r-v_r_normal)) &
         /(qr(3)*(s_r-v_r_normal) - ql(3)*(s_l-v_l_normal))
   
   ! Pressure at right and left (Pressure_j=Pressure_i) side of contact surface
   pStar = qr(3) * (v_r_normal-s_r)*(v_r_normal-s_m) + p_r

   if (s_m >= 0.0) then
      if (s_l > 0.0) then
         flux(1)   = ql(3)*v_l_normal
         flux(2:3) = (ql(3)*v_l_normal)*ql(1:2) + p_l*unorm(1:2)
         flux(4)   = e_l*ql(3)*v_l_normal + p_l*v_l_normal
      else
         invSLmSs = 1.0/(s_l-s_m)
         sLmuL = s_l - v_l_normal
         rhoSL = ql(3)*sLmuL*invSLmSs
         rhouSL(1:2) = ((ql(3)*sLmuL)*ql(1:2) + (pStar-p_l)*unorm(1:2))*invSLmSs
         eSL = (sLmuL*e_l*ql(3) - p_l*v_l_normal + pStar*s_m)*invSLmSs
         
         flux(1)   = rhoSL*s_m
         flux(2:3) = rhouSL(1:2)*s_m + pStar*unorm(1:2)
         flux(4)   = (eSL + pStar)*s_m
      endif
   else
      if (s_r >= 0.0) then
         invSRmSs = 1.0/(s_r-s_m)
         sRmuR = s_r - v_r_normal
         rhoSR = qr(3)*sRmuR*invSRmSs
         rhouSR(1:2) = ((qr(3)*sRmuR)*qr(1:2) + (pStar-p_r)*unorm(1:2))*invSRmSs
         eSR = (sRmuR*e_r*qr(3) - p_r*v_r_normal + pStar*s_m)*invSRmSs
         
         flux(1)   = rhoSR*s_m
         flux(2:3) = rhouSR(1:2)*s_m + pStar*unorm(1:2)
         flux(4)   = (eSR + pStar)*s_m
      else
         flux(1)   = qr(3)*v_r_normal
         flux(2:3) = (qr(3)*v_r_normal)*qr(1:2) + p_r*unorm(1:2)
         flux(4)   = e_r*qr(3)*v_r_normal + p_r*v_r_normal
      endif
   endif
   
end subroutine numflux

!-----------------------------------------------
!Monitor Function - Called after each time-step
!-----------------------------------------------
subroutine Monitor(ts, step, time, ug, user, ierr)
#include <petsc/finclude/petscts.h>
   use petscts
   use comvar
   implicit none
   TS                   :: ts
   Vec                  :: ug
   PetscReal, intent(in):: user(5)
   PetscInt             :: step
   PetscErrorCode       :: ierr
   PetscReal            :: time
   PetscViewer          :: viewer
   character(len = 100) :: fname

   call TSGetTime(ts, time, ierr); CHKERRQ(ierr)
   call TSGetTimeStepNumber(ts, step, ierr); CHKERRQ(ierr)
   call TSSetTimeStep(ts, 0.001, ierr); CHKERRQ(ierr)

   if(id == 0) print*, "Steps, time = ", step, time

   if(mod(step, itplot) == 0) then
      solfile = solfile + 1
      fname = 'Output/Solution'
      write(unit=fname, fmt = '(A15, I0.6)') trim(fname), solfile
      fname = trim(fname)//'.xml'
      call PetscViewerCreate(comm, viewer, ierr)
      call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr)
      call PetscViewerFileSetMode(viewer, PETSCVIEWERASCII, ierr)
      call PetscViewerASCIIOpen(comm, fname, viewer, ierr)
      call PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB, ierr)
      call VecView(ug, viewer, ierr); CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
   endif

   return
end subroutine Monitor
