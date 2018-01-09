SUBROUTINE setupgrid
  
   USE mod_precdef
   USE netcdf
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_getfile
   
   IMPLICIT none
   ! =============================================================
   !    ===  Set up the grid for ORCA0083 configuration ===
   ! =============================================================
   ! Subroutine for defining the grid of the ORCA0083 config. 
   ! Run once before the loop starts.
   ! -------------------------------------------------------------
   ! The following arrays will be populated:
   !
   !  dxdy - Horizontal area of cells (T points)
   !  dz   - Thickness of standard level (T point) 
   !  dzt  - Time-invariant thickness of level (T point)
   !  dzu  - Time-invariant thickness of level (U point)
   !  dzv  - Time-invariant thickness of level (V point)
   !  kmt  - Number of levels from surface to seafloor (T point)
   !  kmu  - Number of levels from surface to seafloor (U point)
   !  kmv  - Number of levels from surface to seafloor (V point)
   !
   ! -------------------------------------------------------------
    
   ! === Init local variables for the subroutine ===
   INTEGER                                      :: i ,j ,k, n, kk, ii
   INTEGER                                      :: ip, jp, im, jm
   REAL(DP), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t    !! dx, dy [m]
   CHARACTER (len=200)                          :: gridFile 
   
   map2D    = [3, 4,  1, 1 ]
   map3D    = [2, 3,  4, 1 ]
   ncTpos   = 1
   !
   ! --- Read dx, dy at T points --- 
   !
   allocate ( e1t(imt,jmt) , e2t(imt,jmt) )
   gridFile = trim(inDataDir)//'GLO-MFC_001_024_coordinates.nc'
   e1t  = get2DfieldNC(gridFile, 'e1t')
   e2t  = get2DfieldNC(gridFile, 'e2t')
   dxdy(1:imt,1:jmt) = e1t(1:imt,1:jmt) * e2t(1:imt,1:jmt)
   !deallocate ( e1t, e2t )
  
   !
   ! --- Read dy at U points and dx at V points --- 
   !
   dyu  = get2DfieldNC(gridFile, 'e2t')
   where (dyu > 10000)
      dyu = 0
   end where
   dxv  = get2DfieldNC(gridFile, 'e1t')
   where (dxv > 10000)
      dxv = 0
   end where
   dx   = dxv(imt/2, jmt/2)
   dy   = dyu(imt/2, jmt/2)
   
   !
   ! Read dz at T points without considering 
   !
   allocate ( dzu(imt,jmt,km,2),dzv(imt,jmt,km,2), dzt0(imt,jmt,km) )
   dzt0(:,:,km:1:-1) = get3DfieldNC(gridFile, 'e3t')
   where (dzt0 > 10000)
      dzt0 = 0
   end where
   dzt(:,:,:,1) = dzt0
   dzt(:,:,:,2) = dzt0
   dzu = dzt
   dzv = dzt
   
   kmt = 0
   do k=1,km
      where (dzt0(:,:,k) > 0)
         kmt = kmt + 1
      end where
   end do

   !kmt = get2DfieldNC(gridFile, 'mbathy')
   allocate ( kmu(imt,jmt), kmv(imt,jmt) )
   kmu=0 ; kmv=0
   do j=1,jmt
      jp=j+1
      if(jp == jmt+1) jp=jmt
      do i=1,imt
         ip=i+1
         if(ip == imt+1) ip=1
         kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
         kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
      enddo
   enddo

   do i=4, imt
      ii = imt + 4 - i
      kmv(i,jmt) = kmv(ii,jmt-3)
   enddo
  
   !dzu(:,:,:,1) = get3DfieldNC(gridFile, 'e3t')
   !dzv(:,:,:,1) = get3DfieldNC(gridFile, 'e3t')
         
end SUBROUTINE setupgrid
