SUBROUTINE readfields
  
   USE mod_precdef
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_traj
   USE mod_getfile
   use mod_seed
   use mod_tempsalt
   
   IMPLICIT none
   ! ==========================================================================
   ! === Read velocity, temperature and salinity for ORCA0083 configuration ===
   ! ==========================================================================
   ! Subroutine to read the ocean state from ORCA0083 config
   ! Run each time step
   ! --------------------------------------------------------------------------
   ! The following arrays will be populated:
   !
   ! uflux    - Zonal volume flux (U point)
   ! vflux    - Meridional volume flux (V point)
   !
   ! If run with tempsalt option, the following are also set
   ! tem      - Temperature (T point) 
   ! sal      - Salinity (T point)
   ! rho      - Potential density (T point)
   !
   ! --------------------------------------------------------------------------
   
   ! = Loop variables
   INTEGER                                       :: i, j, k ,kk, im, ip
   INTEGER                                       :: jm, jp, imm, ii
   INTEGER                                       :: jmm, jpp, l
   INTEGER                                       :: kbot,ktop
   INTEGER, SAVE                                 :: ntempus=0,ntempusb=0,nread
   ! = Variables used for getfield procedures
   CHARACTER (len=200)                           :: filename, dataprefix
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: xxx
   REAL*4                                      :: dd,hu,hv,uint,vint,zint,hh,h0
  
   LOGICAL                                       :: around
 
!---------------------------------------------------------------   
   !
   ! Allocate variables 
   !
   !alloCondUVW: if(.not. allocated (zstot)) then
   !   allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
   !   allocate ( xxx(imt,jmt,km))
   !endif alloCondUVW
 
   call datasetswap ! Swap between current and previous step
   call updateClock 

   ! Create filename for the current timestep
   dataprefix='mercator_phy_XXXXXX.nc'
   write(dataprefix(14:19),'(I6)') int(loopJD)
   filename = trim(indatadir)//trim(dataprefix)

   ! Read velocity fields
   uvel(:,:,km:1:-1) = get3DfieldNC(trim(filename), 'uo')
   vvel(:,:,km:1:-1) = get3DfieldNC(trim(filename), 'vo')
   where (uvel==-32767)
      uvel = 0
   end where
   where (vvel==-32767)
      vvel = 0
   end where
   uvel = uvel * 0.000610370188951
   vvel = vvel * 0.000610370188951
   ! Calculate volume fluxes
   do k = 1, km
      uflux(:,:,k,nsp)  = uvel(:imt,:,k) * dyu(:imt,:) * dzu(:,:,k,nsp)
      vflux(:,1:,k,nsp) = vvel(:imt,:,k) * dxv(:imt,:) * dzv(:,:,k,nsp)
   enddo

   ! Check that volume fluxes are zero below sea floor
   !do i=1,IMT
   !do j=1,JMT
   !do k=1,KM
   !   if(k > kmv(i,j) .and. vflux(i,j,km+1-k,nsp) /= 0.) then
   !      print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
   !      stop 4966
   !   endif
   !   if(k > kmu(i,j) .and. uflux(i,j,km+1-k,nsp) /= 0.) then
   !      print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
   !      stop 4967
   !   endif
   !enddo
   !enddo
   !enddo
   

   return
   
end subroutine readfields



