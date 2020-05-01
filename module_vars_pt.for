  !##########################################################################
      module vars_pt
!##########################################################################
	  SAVE
      integer np,tsteps_pt,tsnr,cnt_pt,np_loc,ptnr
      logical DF,PSIcell,random,random_dp	
      integer,allocatable,dimension(:)::	ptsinproc
!	integer,allocatable,dimension(:)::  ipux,jpvy,kpwz
      integer,allocatable,dimension(:)::  id
      integer,allocatable,dimension(:):: ibs_u,ibe_u
      integer,allocatable,dimension(:):: jbs_u,jbe_u
      integer,allocatable,dimension(:):: kbs_u,kbe_u
      integer,allocatable,dimension(:):: ibs_v,ibe_v
      integer,allocatable,dimension(:):: jbs_v,jbe_v
      integer,allocatable,dimension(:):: kbs_v,kbe_v
      integer,allocatable,dimension(:):: ibs_w,ibe_w
      integer,allocatable,dimension(:):: jbs_w,jbe_w
      integer,allocatable,dimension(:):: kbs_w,kbe_w
      double precision, allocatable, dimension(:):: xp_pt,yp_pt,zp_pt
      double precision, allocatable, dimension(:):: uoi_pt,voi_pt,woi_pt
      double precision, allocatable, dimension(:):: ui_pt,vi_pt,wi_pt
      double precision, allocatable, dimension(:):: xp_loc,yp_loc,zp_loc
      double precision, allocatable, dimension(:):: uop_pt,vop_pt,wop_pt
      double precision, allocatable, dimension(:):: uop_loc,vop_loc
      double precision, allocatable, dimension(:):: wop_loc,dp_loc
      double precision, allocatable, dimension(:):: Fpu,Fpv,Fpw
      double precision, allocatable, dimension(:):: xpold,ypold,zpold
      double precision, allocatable, dimension(:):: uopold,vopold,wopold
      double precision, allocatable, dimension(:):: uoiold,voiold,woiold
      double precision, allocatable, dimension(:):: dp_pt,dp_old
      double precision xp,yp,zp,uop,vop,wop,div,Dp,sigma


      end module

