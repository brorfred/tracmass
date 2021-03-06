SUBROUTINE readfields

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

#ifdef tempsalt
  USE mod_dens
  USE mod_stat
#endif
  IMPLICIT none
  
  
  ! = Loop variables
  INTEGER                                      :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
  INTEGER                                      :: kbot,ktop
  INTEGER, SAVE                                :: ntempus=0,ntempusb=0,nread
  ! = Variables used for getfield procedures
  CHARACTER (len=200)                          :: fieldFile
  ! = Variables for filename generation
  CHARACTER (len=200)                          :: dataprefix
  REAL*8, ALLOCATABLE, DIMENSION(:,:)          :: zstot,zstou,zstov
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:)        :: xxx
  REAL*4 dd,hu,hv,uint,vint,zint,hh,h0
  
#ifdef initxyt
  INTEGER, PARAMETER :: NTID=73
  INTEGER, PARAMETER :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
  REAL*4 , SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
#endif

  
#ifdef tempsalt
 REAL*4, ALLOCATABLE, DIMENSION(:)             :: rhozvec, depthzvec, latvec
  REAL*4 ,       ALLOCATABLE, DIMENSION(:)     :: tmpzvec, salzvec
#endif

 LOGICAL around

!---------------------------------------------------------------

#ifdef initxyt
  alloCondGrid: if ( .not. allocated (ntimask) ) then
     allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
  endif alloCondGrid
#endif
 alloCondUVW: if(.not. allocated (zstot)) then
    allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
    allocate ( xxx(imt,jmt,km))
#ifdef tempsalt
   allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
#endif
 endif alloCondUVW
 
 call datasetswap
 call updateClock

! === Initialising fields ===
 initFieldcond: if(ints.eq.intstart) then

#ifdef initxyt
    ! Time for individual start positions
    if(IJKMAX2.eq.7392) open(84,file=trim(inDataDir)//'topo/masktime_32_025', &
         form='unformatted')
    read(84) trajinit
    close(84)
    j=0
    do k=1,NTID
       do i=1,IJKMAX2
          if(trajinit(k,i,3).ne.0.) then
             j=j+1
#if orca025l75h6
             trajinit(k,i,3)=float(kst2)-0.5
             !   print *,j,trajinit(k,i,:)
#endif
          endif
       enddo
       ! print *,k,j
    enddo
    ijkst=0
    ! print *,'ijkmax=',j,IJKMAX2,ijkmax
    if(j.ne.IJKMAX2) then
       stop 4396
    endif
#endif

 endif initFieldcond

dataprefix='xxxx/ORCA0083-N01_xxxxxxxx'
write(dataprefix(1:4),'(i4)')   currYear
write(dataprefix(19:26),'(i4,i2.2,i2.2)') currYear,currMon,currDay
fieldFile = trim(inDataDir)//'fields/'//trim(dataprefix)//'d05'

hs(:,     :, nsp) = get2DfieldNC(trim(fieldFile)//'T.nc', 'sossheig')
hs(imt+1, :, nsp) = hs(1,:,nsp)



 where (abyst /= 0)
   zstot = hs(:imt,:jmt,nsp)/abyst + 1
 elsewhere
   zstot = 0.d0
 end where
 
  where (abysu /= 0)
   zstou = 0.5*(hs(:imt,:jmt,nsp)+hs(2:imt+1,:jmt,nsp))/abysu + 1
 elsewhere
   zstou = 0.d0
 end where
 
 where (abysv /= 0)
   zstov = 0.5*(hs(:imt,:jmt,nsp)+hs(:imt,2:jmt+1,nsp))/abysv + 1
 elsewhere
   zstov = 0.d0
 end where

!do j=1,JMT
! print *,j,zstot(1,j),zstou(1,j),zstov(1,j)
!enddo
!stop 34965




 do k=1,km
   dzt(:,:,k,nsp) = dztb(:,:,km+1-k) * zstot
 end do


!   do i=1,IMT
!    do j=1,JMT
!    
!     if(kmt(i,j)/=0) then
!      h0= abyst(i,j)   ! total depth
!      hh= h0 + hs(i,j,nsp)                      ! total thickness of water column on T-points
!      do k=1,KM
!       kk = KM+1-k
!       if (kmt(i,j) >= k ) then ! for the levels above the bottom box
!        dzt(i,j,kk,nsp) = dztb(i,j,k)    * hh / h0
!       else
!        dzt(i,j,kk,nsp) = 0.
!       endif
!      enddo
!     endif
!     
!     if(kmu(i,j)/=0) then
!      h0= abysu(i,j)   ! total depth
!      hh= h0 + 0.5*(hs(i,j,nsp)+hs(i+1,j,nsp))  ! total thickness of water column
!      do k=1,KM
!       kk = KM+1-k
!       if (kmu(i,j) >= k ) then ! for the levels above the bottom box
!        dzt(i,j,kk,nsp) = dztb(i,j,k)    * hh / h0
!       else
!        dzt(i,j,kk,nsp) = 0.
!       endif
!      enddo
!     endif
!
!     
!    enddo
!   enddo






#if defined tempsalt 
 xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'votemper')
 tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
 
 xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'T.nc', 'vosaline')
 sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
 
 depthzvec = 0.
 do j=1,JMT
    latvec=-80+1./12.*float(j+subGridJmin-1)
    do i=1,IMT
       tmpzvec = tem(i,j,:,nsp)
       salzvec = sal(i,j,:,nsp)
       call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)
       rho(i,j,:,nsp)=rhozvec - 1000.
    end do
 end do
#endif     

 uvel = get3DfieldNC(trim(fieldFile)//'U.nc', 'vozocrtx')
 vvel = get3DfieldNC(trim(fieldFile)//'V.nc', 'vomecrty')
! do k=1,km
!!    uflux(:,:,km+1-k,nsp) = uvel(:,:,k) !* dyu * dzu(:,:,k) * zstou
!    vflux(:,:,km+1-k,nsp) = vvel(1:imt,:,k) !* dxv * dzv(:,:,k) * zstov
! end do
 
 
 
 do i=1,IMT
 do j=1,JMT
 do k=1,KM
    uflux(i,j,km+1-k,nsp) = uvel(i,j,k) * dyu(i,j) * dzu(i,j,k) * zstou(i,j)
    vflux(i,j,km+1-k,nsp) = vvel(i,j,k) * dxv(i,j) * dzv(i,j,k) * zstov(i,j)
 enddo
 enddo
 enddo


 do i=1,IMT
 do j=1,JMT
 do k=1,KM
  if(k>kmv(i,j) .and. vflux(i,j,km+1-k,nsp)/=0. ) then
   print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
   stop 4966
  endif
  if(k>kmu(i,j) .and. uflux(i,j,km+1-k,nsp)/=0. ) then
   print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
   stop 4967
  endif
 enddo
 enddo
 enddo



#ifdef drifter
! average velocity/transport to simulate drifter trajectories
kbot=65 ; ktop=66 ! number of surface layers to integrate over 
uint=0. ; vint=0. ; zint=0.
do k=kbot,ktop
   uint = uint + uflux(:,:,k,nsp) ! integrated transport
   vint = vint + vflux(:,:,k,nsp)
   zint = zint + dz(k)          ! total depth of drougued drifter
end do
! weighted transport for each layer
do k=kbot,KM
   uflux(:,:,k,nsp) = uint*dz(k)/zint 
   vflux(:,:,k,nsp) = vint*dz(k)/zint
enddo
#endif

 
 !---
 
! Sea surface height
!ierr=NF90_OPEN(trim(fieldFile)//'T.nc',NF90_NOWRITE,ncid)
!if(ierr.ne.0) stop 3766
!ierr=NF90_INQ_VARID(ncid,'sossheig',varid)
!if(ierr.ne.0) stop 3768
!ierr=NF90_GET_VAR(ncid,varid,temp2d_simp,start2d,count2d)
!if(ierr.ne.0) stop 3799
!ierr=NF90_CLOSE(ncid)

!do j=1,JMT
! do i=1,IMT+1
!  ii=i
!  if(ii.eq.IMT+1) ii=1
!  hs(i,j,nsp)=temp2d_simp(ii,j)
! enddo
!enddo



!---------------------------------------------------------------
! Compute the level thickness of all boxes by z-star coordinates
! dz* = dz (H+ssh)/H and the variable bottom box
!do i=1,IMT
!      do j=1,JMT
!         do k=1,KM
!            kk = KM+1-k
!            if (kmt(i,j) == kk) then ! for the bottom box
!               dzt(i,j,k,nsp) = dztb(i,j,1) * (zw(kmt(i,j)) + hs(i,j,nsp)) / zw(kmt(i,j))
!               if(dzt(i,j,k,nsp).eq.0.) then
!                print *,i,j,kmt(i,j),dztb(i,j,1),hs(i,j,nsp),zw(kmt(i,j))
!                stop 4967
!               endif
!            elseif (kmt(i,j) /= 0) then ! for the levels above the bottom box
!               dzt(i,j,k,nsp) = dz(k)       * (zw(kmt(i,j)) + hs(i,j,nsp)) / zw(kmt(i,j))
!               if(dzt(i,j,k,nsp).eq.0.) then
!                print *,i,j,kmt(i,j),dz(k),hs(i,j,nsp),zw(kmt(i,j))
!                stop 4968
!               endif
!            else
!               dzt(i,j,k,nsp) = 0.
!            endif
!         enddo
!      enddo
!   enddo
   
!!$#ifdef tempsalt 
!!$! Temperature
!!$gridFile = trim(fieldFile)//'T.nc'
!!$ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!!$if(ierr.ne.0) stop 5751
!!$ierr=NF90_INQ_VARID(ncid,'votemper',varid) 
!!$if(ierr.ne.0) stop 3769
!!$ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
!!$if(ierr.ne.0) stop 3799
!!$ierr=NF90_CLOSE(ncid)
!!$!do i=1,IMT
!!$! do j=1,JMT
!!$  do k=1,KM 
!!$   kk=KM+1-k
!!$   tem(:,:,kk,nsp)=temp3d_simp(:,:,k)
!!$  enddo
!!$! enddo
!!$!enddo
!!$   
!!$! Salinity
!!$ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!!$if(ierr.ne.0) stop 5751
!!$ierr=NF90_INQ_VARID(ncid,'vosaline',varid) 
!!$if(ierr.ne.0) stop 3769
!!$ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
!!$if(ierr.ne.0) stop 3799
!!$ierr=NF90_CLOSE(ncid)
!!$!do i=1,IMT
!!$! do j=1,JMT
!!$  do k=1,KM
!!$   kk=KM+1-k
!!$   sal(:,:,kk,nsp)=temp3d_simp(:,:,k)
!!$  enddo
!!$! enddo
!!$!enddo
!!$
!!$depthb=0.
!!$do j=1,JMT
!!$ latb=-80+1./12.*float(j+subGridJmin-1)
!!$ do i=1,IMT
!!$  do k=1,kmt(i,j)
!!$   kk=KM+1-k
!!$   tempb(k)=tem(i,j,kk,nsp)
!!$   saltb(k)=sal(i,j,kk,nsp)
!!$  enddo
!!$  call statvd(tempb, saltb, rhob ,KM ,depthb ,latb)
!!$  do k=1,kmt(i,j)
!!$   kk=KM+1-k
!!$   rho(i,j,kk,nsp)=rhob(k)-1000.
!!$  enddo
!!$ enddo
!!$enddo
!!$
!!$#endif     

!!$! u velocity
!!$gridFile = trim(fieldFile)//'U.nc'
!!$ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!!$if(ierr.ne.0) stop 5751
!!$ierr=NF90_INQ_VARID(ncid,'vozocrtx',varid) 
!!$if(ierr.ne.0) stop 3769
!!$ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
!!$if(ierr.ne.0) stop 3799
!!$ierr=NF90_CLOSE(ncid)
!!$
!!$   do i=1,IMT
!!$   	  ip=i+1
!!$   	  if(i.eq.IMT) ip=1
!!$      do j=1,JMT
!!$       hu=min(zw(kmt(i,j)),zw(kmt(ip,j))) ! total depth at u point
!!$         do kk=1,kmu(i,j)
!!$            k=KM+1-kk
!!$            dd = dz(k)
!!$            if (kk == kmu(i,j)) THEN
!!$               dd = MIN (dztb(i,j,nsp),dztb(ip,j,nsp))
!!$            endif
!!$            if (kmu(i,j) <= 0) THEN
!!$               dd = 0.
!!$            endif
!!$            ! thickness of the wall at the u-point
!!$            if (kmu(i,j) /= 0) then
!!$             dd = dd * ( hu + 0.5*(hs(i,j,2) + hs(ip,j,2)) ) / hu 
!!$             uflux(i,j,k,2) = temp3d_simp(i,j,kk) * dyu(i,j) * dd 
!!$            endif
!!$         enddo
!!$      enddo
!!$   enddo


!!$! v velocity
!!$gridFile = trim(fieldFile)//'V.nc'
!!$ierr=NF90_OPEN(trim(gridFile),NF90_NOWRITE,ncid)
!!$if(ierr.ne.0) stop 5751
!!$ierr=NF90_INQ_VARID(ncid,'vomecrty',varid) ! kmt field
!!$if(ierr.ne.0) stop 3770
!!$ierr=NF90_GET_VAR(ncid,varid,temp3d_simp,start3d,count3d)
!!$if(ierr.ne.0) stop 3799
!!$ierr=NF90_CLOSE(ncid)
!!$
!!$
!!$   do i=1,IMT
!!$      do j=1,JMT-1
!!$       hv=min(zw(kmt(i,j)),zw(kmt(i,j+1))) ! total depth at v point
!!$         do kk=1,kmv(i,j)
!!$            k = KM+1-kk
!!$            dd = dz(k)
!!$            if (kk == kmv(i,j)) THEN
!!$               dd = MIN (dztb(i,j,2),dztb(i,j+1,2))
!!$            endif
!!$            if (kmv(i,j) <= 0) then
!!$               dd = 0.
!!$            endif
!!$            ! thickness of the wall at the u-point
!!$            if (kmv(i,j) /= 0) then
!!$             dd = dd * ( hv + 0.5*(hs(i,j,2) + hs(i,j+1,2)) ) / hv 
!!$             vflux(i,j,k,2) = temp3d_simp(i,j,kk) * dxv(i,j) * dd 
!!$            endif
!!$         enddo
!!$      enddo
!!$   enddo
!!$


!!$!  north fold 
!!$do i=4,IMT
!!$ ii=IMT+4-i
!!$! vflux(i,JMT,:,2)=-vflux(ii,JMT-3,:,2)
!!$enddo
!!$
!!$#ifdef drifter
!!$! average velocity/transport over surface drifter drogue depth to simulate drifter trajectories
!!$kbot=65 ; ktop=66 ! number of surface layers to integrat over
!!$do i=1,imt
!!$ do j=1,jmt
!!$ 
!!$  uint=0. ; vint=0. ; zint=0.
!!$  do k=kbot,ktop
!!$   uint=uint+uflux(i,j,k,2) ! integrated transport
!!$   vint=vint+vflux(i,j,k,2)
!!$   zint=zint+dz(k)          ! total depth of drougued drifter
!!$  enddo
!!$  ! weighted transport for each layer
!!$  do k=kbot,KM
!!$    uflux(i,j,k,2)=uint*dz(k)/zint 
!!$    vflux(i,j,k,2)=vint*dz(k)/zint
!!$  enddo
!!$
!!$ enddo
!!$enddo
!!$
!!$#endif


#ifdef initxyt
! Set the initial trajectory positions
!ijkst(:,5)=ntimask(ntempus,:)
#ifdef orca025l75h6
if( mod(ints,24/ngcm*5).eq.1 .or. ints.le.2) ntempus=ntempus+1
if(ntempus.ne.ntempusb .and. ntempus.le.NTID) then
ntempusb=ntempus
!print *,'ints=',ints,' ntempus=',ntempus,' ntempusb=',ntempusb
#else
if(ints.le.NTID) then
#endif

   
do ntrac=1,ijkmax
 if(trajinit(ntempus,ntrac,3).ne.0.) then
  ijkst(ntrac,4)=0
  ijkst(ntrac,5)=5
  ijkst(ntrac,6)=ijkst(ntrac,6)+1
  do l=1,3
   ijkst(ntrac,l)=trajinit(ntempus,ntrac,l)+1
   trj(ntrac,l)=trajinit(ntempus,ntrac,l)
!   if(l.eq.1) print *,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
   if(trj(ntrac,l).gt.float(ijkst(ntrac,l)) .or. trj(ntrac,l).lt.float(ijkst(ntrac,l)-1)) then
    print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
    stop 3946
   endif
  enddo
 else
  ijkst(ntrac,5)=0
  ijkst(ntrac,6)=0
 endif
enddo
endif
#ifdef orca025l75h6

#endif
if( mod(ints,24/ngcm*5).ne.1 .and. ints.gt.2) then
 ijkst=0 
endif
#endif


! deallocate ( temp3d_simp, temp2d_simp )
     
!#ifdef tempsalt
! deallocate ( tempb, saltb, rhob, depthb, latb )
!#endif

  return
end subroutine readfields



