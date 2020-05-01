!##########################################################################
      module imbub
!##########################################################################
	  SAVE
      integer np,tsteps_pt,tsnr,cnt_pt,np_loc,ptnr
      logical random,random_dp	
      integer,allocatable,dimension(:):: imbsinproc,imbinblk,imb_block
      INTEGER,allocatable,dimension(:) :: imb_proc!,bub_proc
      INTEGER,allocatable,dimension(:) :: lag_bod_loc,lag_bod
      INTEGER,allocatable,dimension(:) :: imb_block_loc
!      integer,allocatable,dimension(:):: bub_block,bubsinproc,bubinblk

!      integer,allocatable,dimension(:):: ibs_u,ibe_u
!      integer,allocatable,dimension(:):: jbs_u,jbe_u
!      integer,allocatable,dimension(:):: kbs_u,kbe_u
!      integer,allocatable,dimension(:):: ibs_v,ibe_v
!      integer,allocatable,dimension(:):: jbs_v,jbe_v
!      integer,allocatable,dimension(:):: kbs_v,kbe_v
!      integer,allocatable,dimension(:):: ibs_w,ibe_w
!      integer,allocatable,dimension(:):: jbs_w,jbe_w
!      integer,allocatable,dimension(:):: kbs_w,kbe_w
      double precision, allocatable, dimension(:):: xp_pt,yp_pt,zp_pt
      double precision, allocatable, dimension(:):: uoi_pt,voi_pt,woi_pt
      double precision, allocatable, dimension(:):: ui_pt,vi_pt,wi_pt
      double precision, allocatable, dimension(:):: xp_loc,yp_loc,zp_loc
      double precision, allocatable, dimension(:):: uop_pt,vop_pt,wop_pt
!      double precision, allocatable, dimension(:):: uop_loc,vop_loc
!      double precision, allocatable, dimension(:):: wop_loc!,dp_loc
      double precision, allocatable, dimension(:):: Fpu,Fpv,Fpw
      double precision, allocatable, dimension(:):: xpold,ypold,zpold
      double precision, allocatable, dimension(:):: uopold,vopold,wopold
      double precision, allocatable, dimension(:):: uoiold,voiold,woiold
      double precision, allocatable, dimension(:):: dp_pt,dp_old
      DOUBLE PRECISION,allocatable,dimension (:,:) :: pre_nodex 
      DOUBLE PRECISION,allocatable,dimension (:,:) :: pre_nodey    
      DOUBLE PRECISION,allocatable,dimension (:,:) :: pre_nodez   
      double precision xp,yp,zp,uop,vop,wop,div,Dp,sigma
	DOUBLE PRECISION::dxm,dym,dzm,nxl
	

      INTEGER :: maxnode,maxnodeIBS,nibs_loc!,nbubs_loc

      DOUBLE PRECISION,allocatable,dimension(:)::nodex_loc,nodey_loc
	DOUBLE PRECISION,allocatable,dimension(:)::nodez_loc
      DOUBLE PRECISION,allocatable,dimension(:)::FX1_loc,FX2_loc
      DOUBLE PRECISION,allocatable,dimension(:)::FX3_loc

      DOUBLE PRECISION,allocatable,dimension(:)::U_Beta1_loc
	DOUBLE PRECISION,allocatable,dimension(:)::U_Beta2_loc
      DOUBLE PRECISION,allocatable,dimension(:)::U_Beta3_loc
      DOUBLE PRECISION,allocatable,dimension(:,:)::dh1_loc,dh2_loc
      DOUBLE PRECISION,allocatable,dimension(:,:)::dh3_loc

      DOUBLE PRECISION,allocatable,dimension(:,:)::nodex,nodey,nodez

      DOUBLE PRECISION,allocatable,dimension(:)::Fu_bub,Fv_bub,Fw_bub   	
      DOUBLE PRECISION,allocatable,dimension(:)::U_bub,V_bub,W_bub
       
	INTEGER,allocatable,dimension(:,:) :: I_nr_V,J_nr_V,K_nr_V
	INTEGER,allocatable,dimension(:,:) :: I_nr_U,J_nr_U,K_nr_U
	INTEGER,allocatable,dimension(:,:) :: I_nr_W,J_nr_W,K_nr_W

	INTEGER,allocatable,dimension(:) :: kmaxU,kmaxV,kmaxW,nodes
              
	INTEGER :: nl,cmax     


      end module
