!=======================================================================
!    			ImBubble
!			Bruño Fraga Bugallo
!			Birmingham - CTR Stanford 2018
!=======================================================================
!##########################################################################      
      subroutine imb_bubbles

!	Creates gas bubbles whose interface is defined by IBs

!##########################################################################
        use multidata
        use mpi
        use vars   
	use imbub
        use omp_lib, only : omp_get_num_threads,
     &                    omp_get_thread_num,
     &                    omp_set_num_threads
      implicit none

      	integer L
      	CHARACTER*31 :: gridfile
      	double precision REp,rho_p,PI
      	double precision aslip,bslip,cslip,ao,bo,co,upmod
	double precision epsylon,gamma_p
	double precision Vcell,Vp!,Vball
      	double precision, allocatable, dimension(:):: up_pt,vp_pt,wp_pt
	double precision, allocatable, dimension(:):: Fgw,Fau,Fav,Faw
     

	allocate (up_pt(np),vp_pt(np),wp_pt(np))
	!allocate (Fgw(np),Fau(np),Fav(np),Faw(np))

	rho_p = 1.4d0
	dens = 1000.d0
        PI = 4.D0*DATAN(1.D0)	
	Vcell = g_dx*g_dy*g_dz	

!loop in bubbles
      do L=1,np
		!write(44,*)wi_pt(L),Fpw(L)

	Vp = (PI*dp_pt(L)**3.d0)/6.d0

!	if (LENERGY) then
!		dom(ib)%dens(ip,jp,kp) = 				Can't do stratification yet, I guess this should be implemented in IBM
!     &999.8/(1.+0.000088*(dom(ib)%T(ip,jp,kp)+20.))
!		dom(ib)%mu(ip,jp,kp) = 
!     &2.414d-5*10.d0**(-25.2/(dom(ib)%T(ip,jp,kp)+20.-413.d0))

!	Re=dom(ib)%dens(ip,jp,kp)/dom(ib)%mu(ip,jp,kp)
!	endif

	ui_pt(L)= ui_pt(L)*0.5*Vcell/Vp					!from IMB
	vi_pt(L)= vi_pt(L)*0.5*Vcell/Vp					!from IMB
	wi_pt(L)= wi_pt(L)*0.5*Vcell/Vp					!from IMB

	Fpu(L) = Fpu(L)*0.5*rho_p/dens								!from IMB
	Fpv(L) = Fpv(L)*0.5*rho_p/dens									!from IMB
	Fpw(L) = Fpw(L)*0.5*rho_p/dens		

	!Old slip vel components 
     	ao = uop_pt(L)-uoi_pt(L)
      	bo = vop_pt(L)-voi_pt(L)
	co = wop_pt(L)-woi_pt(L)

	up_pt(L) = uop_pt(L)
	vp_pt(L) = vop_pt(L)
	wp_pt(L) = wop_pt(L)

	
	epsylon=1	
	upmod=0

	do while (epsylon.ge.1d-3)
		!New slip vel components 
	      	aslip = up_pt(L)-ui_pt(L)
	     	bslip = vp_pt(L)-vi_pt(L)
		cslip = wp_pt(L)-wi_pt(L)

		!Added mass
!		Fau(L) = -dens*Vp*0.5*(aslip-ao)/dt
!		Fav(L) = -dens*Vp*0.5*(bslip-bo)/dt
!		Faw(L) = -dens*Vp*0.5*(cslip-co)/dt

		!Buoyancy
!		Fgw(L)= Vp*(dens-rho_p)*9.81

      	up_pt(L) = uop_pt(L)  + dt* (1.d0-rho_p/dens)*gz 
     &	- 0.5 * (aslip-ao)						!Added Mass
     &	- Fpu(L)							!Drag & Lift

!     	up_pt(L) = uop_pt(L) 		 				!Schwarz, Kempe, Frolich
!     &	- (1.d0/rho_p) * 0.5 * (aslip-ao)					
!     &	- (1.d0/(rho_p*Vp*dens)) * Fpu(L)
			
      	vp_pt(L) = vop_pt(L)  
     &	- 0.5 * (bslip-bo)						!Added Mass
     &	- Fpv(L) 							!Drag & Lift	

!     	vp_pt(L) = vop_pt(L)  		 				!Schwarz, Kempe, Frolich		 		
 !    &	- (1.d0/rho_p) * 0.5 * (bslip-bo)						
 !    &	- (1.d0/(rho_p*Vp*dens)) * Fpv(L)

     	wp_pt(L) = wop_pt(L) 						!Buoyancy
     &	- 0.5 * (cslip-co)						!Added Mass
     &	- Fpw(L)


!      	wp_pt(L) = wop_pt(L) + dt * (2.d0 * 9.81) -			!Schwarz et al
!     &	Fpw(L) 									


!	      up_pt(L) = uop_pt(L) + dt* ((-Fpu(L) + Fau(L))/(rho_p*Vp))				
!	      vp_pt(L) = vop_pt(L) + dt* ((-Fpv(L) + Fav(L))/(rho_p*Vp))
!	      wp_pt(L) = wop_pt(L) + dt* ((-Fpw(L) + Faw(L) + Fgw(L))/
 !    &			(rho_p*Vp))


	      epsylon = 
     &	abs((upmod-sqrt(up_pt(L)**2.d0+vp_pt(L)**2.d0+wp_pt(L)**2.d0))
     &	/sqrt(up_pt(L)**2.d0+vp_pt(L)**2.d0+wp_pt(L)**2.d0))

	      upmod=sqrt(up_pt(L)**2.d0+vp_pt(L)**2.d0+wp_pt(L)**2.d0)
	enddo	

!	Actualizar velocidad paso previo
	uop_pt(L) = up_pt(L)
      	vop_pt(L) = vp_pt(L)
      	wop_pt(L) = wp_pt(L)

	uoi_pt(L) = ui_pt(L)
      	voi_pt(L) = vi_pt(L)
      	woi_pt(L) = wi_pt(L)

!	Actualizar posicion de particula
      	xp_pt(L)=xp_pt(L)+up_pt(L)*dt
      	yp_pt(L)=yp_pt(L)+vp_pt(L)*dt
      	zp_pt(L)=zp_pt(L)+wp_pt(L)*dt

!     		write(myrank+700,*)'xp_loc',nodex_loc(L),nodey_loc(L),nodez_loc(L)

	if (L.eq.1) write(44,*)itime,xp_pt(L),up_pt(L),vp_pt(L)
     &	,wp_pt(L),Fpu(L),Fpv(L),Fpw(L)
      end do  								!end of loop in bubbles

!      deallocate (Fgw,Fau,Fav,Faw)
!      deallocate (up_pt,vp_pt,wp_pt)

      return
      end subroutine
!######################################################################
      SUBROUTINE bub_algorithm
!######################################################################
      use vars
      use mpi
      use multidata
      use imbub
      implicit none
      INTEGER :: I,J,K,L,M,ib,nt,nnmls,KK,iii,s
      DOUBLE PRECISION :: dh!,dhtotal
      DOUBLE PRECISION :: PI,UIB_loc,VIB_loc,WIB_loc,fbeta

      	if (nibs_loc.gt.0) then

	U_Beta1_loc=0.d0 ; U_Beta2_loc=0.d0 ; U_Beta3_loc=0.d0	

	call OMP_SET_NUM_THREADS(OMP_threads)

	Do ib=1,nbp

!Using delta functions. Interpolation of U-velocities


      Do L = 1,nibs_loc
	 nl=0 
	IF(imb_block_loc(L).eq.dom_id(ib)) then 

          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%x(i) .le.(nodex_loc(L)+nxl*dom(ib)%dx) .and.
     &     dom(ib)%x(i) .ge.(nodex_loc(L)-nxl*dom(ib)%dx)) then
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).le.(nodey_loc(L)+nxl*dom(ib)%dy) .and.
     &     dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy)) then
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).le.(nodez_loc(L)+nxl*dom(ib)%dz) .and.
     &     dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz)) then 
	 nl=nl+1

        dh1_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),order)  
        U_Beta1_loc(L)=U_Beta1_loc(L)+dom(ib)%USTAR(I,J,K)*dh1_loc(L,nl)     
        KmaxU(L)=nl ; I_nr_U(L,nl)=I ; J_nr_U(L,nl)=J ; K_nr_U(L,nl)=K    

		if (order.eq.2 .and. nl.ge.125) GOTO 700
		if (order.eq.4 .and. nl.ge.125) GOTO 700
		if (order.eq.3 .and. nl.ge.64)  GOTO 700
		if (order.eq.6 .and. nl.ge.64)  GOTO 700
		if (order.eq.1 .and. nl.ge.27)  GOTO 700
		if (order.eq.5 .and. nl.ge.27)  GOTO 700
 
           endif
           END DO
           endif
           END DO
           endif
           END DO        
       	   if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'
        
700	  endif								!point belongs to block
      	  Enddo								!loop in ibs



      Do L = 1,nibs_loc
	IF(imb_block_loc(L).eq.dom_id(ib)) then
	 nl=0 

          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).le.(nodex_loc(L)+nxl*dom(ib)%dx) .and.
     &     dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx)) then
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%y(j) .le.(nodey_loc(L)+nxl*dom(ib)%dy) .and.
     &     dom(ib)%y(j) .ge.(nodey_loc(L)-nxl*dom(ib)%dy)) then 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).le.(nodez_loc(L)+nxl*dom(ib)%dz) .and.
     &     dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz)) then
	 nl=nl+1
        dh2_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),order) 		!June 2015
        U_Beta2_loc(L)=U_Beta2_loc(L)+dom(ib)%VSTAR(I,J,K)*dh2_loc(L,nl)         
        KmaxV(L)=nl ; I_nr_V(L,nl)=I ; J_nr_V(L,nl)=J ; K_nr_V(L,nl)=K


	if (order.eq.2 .and. nl.ge.125) GOTO 701
	if (order.eq.4 .and. nl.ge.125) GOTO 701
	if (order.eq.3 .and. nl.ge.64)  GOTO 701
	if (order.eq.6 .and. nl.ge.64)  GOTO 701
	if (order.eq.1 .and. nl.ge.27)  GOTO 701
	if (order.eq.5 .and. nl.ge.27)  GOTO 701
 
	  endif
          END DO
          endif
          END DO
          endif
          END DO
       if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'         

701	endif

      Enddo

      Do L = 1,nibs_loc
	IF(imb_block_loc(L).eq.dom_id(ib)) then 
	 nl=0 

          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).le.(nodex_loc(L)+nxl*dom(ib)%dx) .and.
     &     dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx)) then 
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).le.(nodey_loc(L)+nxl*dom(ib)%dy) .and.
     &     dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy)) then 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%z(k) .le.(nodez_loc(L)+nxl*dom(ib)%dz) .and.
     &     dom(ib)%z(k) .ge.(nodez_loc(L)-nxl*dom(ib)%dz)) then
	 nl=nl+1
        dh3_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%YC(J),dom(ib)%Z(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),order)  		
        U_Beta3_loc(L)=U_Beta3_loc(L)+dom(ib)%WSTAR(I,J,K)*dh3_loc(L,nl)
        KmaxW(L)=nl ;  I_nr_W(L,nl)=I ; J_nr_W(L,nl)=J ; K_nr_W(L,nl)=K


	if (order.eq.2 .and. nl.ge.125) GOTO 702		
	if (order.eq.4 .and. nl.ge.125) GOTO 702
	if (order.eq.3 .and. nl.ge.64)  GOTO 702
	if (order.eq.6 .and. nl.ge.64)  GOTO 702
	if (order.eq.1 .and. nl.ge.27)  GOTO 702
	if (order.eq.5 .and. nl.ge.27)  GOTO 702
 
	    endif
            END DO
	   endif
           END DO
	  endif
          END DO
        if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'

702	endif
      Enddo

	
	Enddo 								!mpi domains loop
        endif								!domain has ibs

!#################   SUBROUTINE calfl   #################################
	Fu_bub = 0.d0      ; Fv_bub = 0.d0      ; Fw_bub = 0.d0
	U_bub = 0.d0	   ; V_bub = 0.d0       ; W_bub = 0.d0

      	if (nibs_loc.gt.0) then

	FX1_loc = 0.d0     ; FX2_loc = 0.d0	; FX3_loc = 0.d0

	DO ib=1,nbp
        Do L = 1,nibs_loc
	  IF(imb_block_loc(L).eq.dom_id(ib)) then

		M=lag_bod_loc(L)

	   	UIB_loc = uop_pt(M)
	   	VIB_loc = vop_pt(M)
	   	WIB_loc = wop_pt(M)

		FX1_loc(L)=UIB_loc-U_Beta1_loc(L)
		FX2_loc(L)=VIB_loc-U_Beta2_loc(L)
		FX3_loc(L)=WIB_loc-U_Beta3_loc(L)  

!		if (lag_bod_loc(M).eq.dom_id(ib)) then			!ibs only contribute to bubbles in the same domain

		Fu_bub(M)=FX1_loc(L)+Fu_bub(M)
		Fv_bub(M)=FX2_loc(L)+Fv_bub(M)
		Fw_bub(M)=FX3_loc(L)+Fw_bub(M)
		U_bub(M)=U_Beta1_loc(L)+U_bub(M)
		V_bub(M)=U_Beta2_loc(L)+V_bub(M)
		W_bub(M)=U_Beta3_loc(L)+W_bub(M)

!		endif

	  endif
	  ENDDO    
	Enddo !ib-loop	                   
	endif


!!################   SUBROUTINE distfbeta   ###################
      	if (nibs_loc.gt.0) then
	Do ib=1,nbp
        Do L = 1,nibs_loc

	IF(imb_block_loc(L).eq.dom_id(ib)) then
	 Do nl=1,KmaxU(L)

          I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
          fbeta = FX1_loc(L)*dh1_loc(L,nl)   
          dom(ib)%USTAR(I,J,K) = dom(ib)%USTAR(I,J,K) + fbeta*0.5

  	 Enddo	
	 Do nl=1,KmaxV(L)
	   I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl);  K=K_nr_V(L,nl)
          fbeta = FX2_loc(L)*dh2_loc(L,nl) 
          dom(ib)%VSTAR(I,J,K) = dom(ib)%VSTAR(I,J,K) + fbeta*0.5

  	 Enddo
	 Do nl=1,KmaxW(L)
	   I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl);  K=K_nr_W(L,nl)
          fbeta = FX3_loc(L)*dh3_loc(L,nl)     
          dom(ib)%WSTAR(I,J,K) = dom(ib)%WSTAR(I,J,K) + fbeta*0.5
  	 Enddo

	endif
        End do 

	Enddo !ib-loop 
        endif

! Send these forces and velocities to the master processor			(Brunho 2018)

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_REDUCE (Fu_bub,Fpu,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE (Fv_bub,Fpv,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE (Fw_bub,Fpw,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        call MPI_REDUCE (U_bub,ui_pt,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE (V_bub,vi_pt,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE (W_bub,wi_pt,np,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

	if (nibs_loc.gt.0) then
	       deallocate (kmaxU,kmaxV,kmaxW)
	       deallocate (nodex_loc,U_Beta1_loc,FX1_loc)
	       deallocate (nodey_loc,U_Beta2_loc,FX2_loc)
	       deallocate (nodez_loc,U_Beta3_loc,FX3_loc)
	       deallocate (imb_block_loc,lag_bod_loc)
	       deallocate (dh1_loc,dh2_loc,dh3_loc)
	       deallocate (I_nr_U,J_nr_U)
	       deallocate (I_nr_V,J_nr_V)
	       deallocate (I_nr_W,J_nr_W)
	       deallocate (K_nr_U,K_nr_V)
	       deallocate (K_nr_W)
	endif

      RETURN
      END
!######################################################################
      SUBROUTINE MPI_bub 
!######################################################################
      use vars
      use imbub
      use mpi
      use multidata
      implicit none

	INTEGER 	:: M,L,ii,nxdom,nydom,nzdom,ntot,s,procprev,o,jj
 	double precision :: lx,ly,lz
	integer,dimension(nprocs) :: strider_imb!,strider_bub
      	double precision, allocatable, dimension(:):: X_MPI,Y_MPI,Z_MPI
      	integer, allocatable, dimension(:):: lag_bod_MPI
      	integer, allocatable, dimension(:):: imb_block_MPI
!	call MPI_BCAST(maxnodeIBS,1,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)


       IF (myrank.eq.0) THEN

!	ntot=np*maxnodeIBS

	allocate(X_MPI(maxnodeIBS),Y_MPI(maxnodeIBS),Z_MPI(maxnodeIBS))
	allocate(lag_bod_MPI(maxnodeIBS),imb_block_MPI(maxnodeIBS))
	
	lx=xen-xst +2.d0 * pl *dxm								
	ly=yen-yst +2.d0 * pl *dym
	lz=zen-zst +2.d0 * pl *dzm

	imbinblk=0  							!imbs in each block
!	bubinblk=0
	imbsinproc=0							!Processor to which every imb belongs
!	bubsinproc=0
	ii=0  				
	do M=1,np  							!All bubbles
!		if (nodex(M,L).lt.xst.or.nodex(M,L).gt.xen) then
!			write(6,*)'ib out of domain in X',M,L,nodex(M,L)
!		endif
!		if (nodey(M,L).lt.yst.or.nodey(M,L).gt.yen) then
!			write(6,*)'ib out of domain in Y',M,L,nodey(M,L)
!		endif
!		if (nodez(M,L).lt.zst.or.nodez(M,L).gt.zen) then
!			write(6,*)'ib out of domain in Z',M,L,nodez(M,L)
!		endif

	    DO L=1,nodes(M)						!Analyze all IB points of the bubble
		ii=ii+1	
		nxdom=INT((nodex(M,L)-(xst-pl*dxm)-1d-12)/(lx/idom))				!Domain where point belongs in i
		nydom=INT((nodey(M,L)-(yst-pl*dym)-1d-12)/(ly/jdom))				!Domain where point belongs in j
		nzdom=INT((nodez(M,L)-(zst-pl*dzm)-1d-12)/(lz/kdom))				!Domain where point belongs in k

		imb_block(ii)=idom*jdom*nzdom+idom*nydom+nxdom		!MPI domain to which the ib belongs
	
!     		write(6,*)'part',L,ii,nxdom,nydom,nzdom,imb_block(ii)
!     		write(6,*)'coord',nodex(M,L),nodey(M,L),nodez(M,L)
!     		write(6,*)'domain',lx,ly,lz

		if (imb_block(ii).ge.num_domains)
     &		write(6,*)'allocation error',M,L,nxdom,nydom,nzdom

!		imbinblk(imb_block(ii)+1)=imbinblk(imb_block(ii)+1)+1	!number of ibs in the MPI dom
		imb_proc(ii)=dom_ad(imb_block(ii))+1 			!Processor to which the ib belongs (1 to num_procs+1)
		imbsinproc(imb_proc(ii))=imbsinproc(imb_proc(ii))+1	!number of ibs in each processor (for SCATTERV)
!	   	lag_bod(ii)=M 						!Bubble to which the ib belongs (number of points per bubble is nodes(M))
	    ENDDO
	   nxdom=INT((xp_pt(M)-xst-1d-12)/(lx/idom))				!Domain where point belongs in i
	   nydom=INT((yp_pt(M)-yst-1d-12)/(ly/jdom))				!Domain where point belongs in j
	   nzdom=INT((zp_pt(M)-zst-1d-12)/(lz/kdom))				!Domain where point belongs in k

!	   bub_block(M)=idom*jdom*nzdom+idom*nydom+nxdom		!MPI domain to which the bub belongs
     		!write(6,*)'part',M,nxdom,nydom,nzdom,bub_block(M)
     		!write(6,*)'coord',xp_pt(M),yp_pt(M),zp_pt(M)
     		!write(6,*)'dom',xst,xen,yst,yen,zst,zen
!		if (bub_block(M).ge.num_domains)
!     &		write(6,*)'allocation error',M,nxdom,nydom,nzdom
!	   bubinblk(bub_block(M)+1)=bubinblk(bub_block(M)+1)+1		!number of bubs in the MPI dom
!	   bub_proc(M)=dom_ad(bub_block(M))+1 				!Processor to which the bub belongs (1 to num_procs+1)
!	   bubsinproc(bub_proc(M))=bubsinproc(bub_proc(M))+1		!number of bus in each processor (for SCATTERV)
	enddo 

	strider_imb(1) = 0						!strider array for SCATTERV
!	strider_bub(1) = 0
	do s=2,nprocs
		strider_imb(s) = imbsinproc(s-1) + strider_imb(s-1)
!		strider_bub(s) = bubsinproc(s-1) + strider_bub(s-1)
	enddo

!create superarrays with the ibs ordered by processors
	do o=1,nprocs
	ii=0								!non-ordered
	jj=0								!ordered
	do M=1,np	;	DO L=1,nodes(M)	
	ii=ii+1		
	if (imb_proc(ii).eq.o) then
	jj=jj+1	
	X_MPI(jj + strider_imb(imb_proc(ii)))=nodex(M,L)				
	Y_MPI(jj + strider_imb(imb_proc(ii)))=nodey(M,L)
	Z_MPI(jj + strider_imb(imb_proc(ii)))=nodez(M,L)

	lag_bod_MPI(jj + strider_imb(imb_proc(ii))) = M

	imb_block_MPI(jj + strider_imb(imb_proc(ii)))=imb_block(ii)
	endif
	ENDDO	;	enddo	
	enddo

	ENDIF								!master

!	if (myrank.eq.0) then 
!	ii=0								!non-ordered
!	do M=1,np	;	DO L=1,nodes(M)	
!	ii=ii+1
!		write(33,*) lag_bod(ii),lag_bod_MPI,imb_block_MPI(ii)
!	ENDDO	;	enddo	
!	endif

!	FOR ALL DOMAINS.................................................
	

!	call MPI_BCAST(imbsinproc,nprocs,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)
!	call MPI_BCAST(bubsinproc,nprocs,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)
!	call MPI_BCAST(imbinblk,num_domains,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)
!	call MPI_BCAST(imb_proc,num_domains,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)
!	call MPI_BCAST(imb_block,maxnodeIBS,MPI_INTEGER,0,
!     &			MPI_COMM_WORLD,ierr)
!	call MPI_SCATTER(bubsinproc,1,MPI_INTEGER,nbubs_loc,1,		!local number of bubs in each processor
!     &		MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

	call MPI_SCATTER(imbsinproc,1,MPI_INTEGER,nibs_loc,1,		!local number of ibs in each processor
     &		MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

!	allocate (uop_loc(nbubs_loc),vop_loc(nbubs_loc))
!	allocate (wop_loc(nbubs_loc))

	if (nibs_loc.gt.0) then
	allocate (imb_block_loc(nibs_loc))
        allocate (kmaxU(nibs_loc),kmaxV(nibs_loc),kmaxW(nibs_loc))
        allocate (nodex_loc(nibs_loc),U_Beta1_loc(nibs_loc))
        allocate (FX1_loc(nibs_loc))
        allocate (nodey_loc(nibs_loc),U_Beta2_loc(nibs_loc))
        allocate (FX2_loc(nibs_loc))
        allocate (nodez_loc(nibs_loc),U_Beta3_loc(nibs_loc))
        allocate (FX3_loc(nibs_loc),lag_bod_loc(nibs_loc))


	kmaxU=0 	   ; kmaxV=0 		; kmaxW=0	         		  
	U_Beta1_loc=0.d0 ; U_Beta2_loc =0.d0 	; U_Beta3_loc =0.d0 
	FX1_loc =0.d0    ; FX2_loc =0.d0    	; FX3_loc =0.d0 
	imb_block_loc=0  	
	Fu_bub=0.d0
	Fv_bub=0.d0
	Fw_bub=0.d0
	U_bub=0.d0
	V_bub=0.d0
	W_bub=0.d0
	endif

	if (order.eq.2 .or. order.eq.4) then
	 nxl=2.499999d0 	; nl=125
	endif
	if (order.eq.3 .or. order.eq.6) then
	 nxl=1.999999d0  	; nl=64 
	endif
	if (order.eq.1 .or. order.eq.5) then
	 nxl=1.499999d0  	; nl=27
	endif

	if (nibs_loc.gt.0) then
       	allocate (dh1_loc(nibs_loc,nl),dh2_loc(nibs_loc,nl))
       	allocate (dh3_loc(nibs_loc,nl))
	 allocate (I_nr_U(nibs_loc,nl),J_nr_U(nibs_loc,nl))
	 allocate (I_nr_V(nibs_loc,nl),J_nr_V(nibs_loc,nl))
	 allocate (I_nr_W(nibs_loc,nl),J_nr_W(nibs_loc,nl))
	 allocate (K_nr_U(nibs_loc,nl),K_nr_V(nibs_loc,nl))
	 allocate (K_nr_W(nibs_loc,nl))
	

	  dh1_loc=0.d0		; dh2_loc=0.d0		; dh3_loc=0.d0
	  I_nr_U=0		; J_nr_u=0 		; K_nr_U=0
	  I_nr_V=0 		; J_nr_V=0 		; K_nr_V=0
	  I_nr_W=0 		; J_nr_W=0 		; K_nr_W=0
	endif

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

	call MPI_SCATTERV(imb_block_MPI,imbsinproc,strider_imb,
     &	MPI_INTEGER,imb_block_loc,nibs_loc,MPI_INTEGER,0
     &	,MPI_COMM_WORLD,ierr)

        call MPI_SCATTERV(X_MPI,imbsinproc,strider_imb,		
     &	MPI_DOUBLE_PRECISION,nodex_loc,nibs_loc,MPI_DOUBLE_PRECISION,0,
     &	MPI_COMM_WORLD,ierr)
        call MPI_SCATTERV(Y_MPI,imbsinproc,strider_imb,		
     &	MPI_DOUBLE_PRECISION,nodey_loc,nibs_loc,MPI_DOUBLE_PRECISION,0,
     &	MPI_COMM_WORLD,ierr)		
        call MPI_SCATTERV(Z_MPI,imbsinproc,strider_imb,		
     &	MPI_DOUBLE_PRECISION,nodez_loc,nibs_loc,MPI_DOUBLE_PRECISION,0,
     &	MPI_COMM_WORLD,ierr)

        call MPI_SCATTERV(lag_bod_MPI,imbsinproc,strider_imb,		
     &	MPI_INTEGER,lag_bod_loc,nibs_loc,MPI_INTEGER,
     &	0,MPI_COMM_WORLD,ierr)

      	call MPI_BCAST(uop_pt,np,MPI_DOUBLE_PRECISION
     &		,0,MPI_COMM_WORLD,ierr)
!	write(6,*)'vop'
      	call MPI_BCAST(vop_pt,np,MPI_DOUBLE_PRECISION
     &		,0,MPI_COMM_WORLD,ierr)
!	write(6,*)'wop'
      	call MPI_BCAST(wop_pt,np,MPI_DOUBLE_PRECISION
     &		,0,MPI_COMM_WORLD,ierr)	
!        call MPI_SCATTERV(uop_pt,bubsinproc,strider_bub,		
!     &	MPI_DOUBLE_PRECISION,uop_loc,nbubs_loc,MPI_DOUBLE_PRECISION,0,
!     &	MPI_COMM_WORLD,ierr)
!        call MPI_SCATTERV(vop_pt,bubsinproc,strider_bub,		
!     &	MPI_DOUBLE_PRECISION,vop_loc,nbubs_loc,MPI_DOUBLE_PRECISION,0,
!     &	MPI_COMM_WORLD,ierr)	
!        call MPI_SCATTERV(wop_pt,bubsinproc,strider_bub,		
!     &	MPI_DOUBLE_PRECISION,wop_loc,nbubs_loc,MPI_DOUBLE_PRECISION,0,
!     &	MPI_COMM_WORLD,ierr)	
	call MPI_BARRIER (MPI_COMM_WORLD,ierr)

	if (myrank.eq.0) then
		deallocate(X_MPI,Y_MPI,Z_MPI)
		deallocate(lag_bod_MPI,imb_block_MPI)
	endif

      RETURN
      END
!######################################################################
      SUBROUTINE imb_spheric_bub(M)
!######################################################################
      use vars
      use multidata
      use imbub
      use mpi
      implicit none
      DOUBLE PRECISION:: thc,PI,thz(np),zccos,zcsin,R
      INTEGER      :: M,L,K,c
      INTEGER      :: strlen2,strlen3,nzr(np),izr,maxnzr
      CHARACTER*8  :: char_block2,char_block3
      CHARACTER*31 :: gridfile
      DOUBLE PRECISION,allocatable,dimension (:,:) :: ztemp_layer
      DOUBLE PRECISION,allocatable,dimension (:,:,:) :: Rtemp_layer      
      INTEGER,allocatable,dimension (:,:) :: ctot_layer,nodes_layer
      INTEGER,allocatable,dimension (:,:,:) :: nodes_percyl_layer


      PI = 4.D0*DATAN(1.D0)

!	if ((mod(itime,tsteps_pt).eq.0).and.
!     &		(itime.gt.itime_start)) then 			

!        write(char_block2,'(I8)') M					!bubble number
!        write(char_block3,'(I8)') cnt_pt				!toutput
!        strlen2=LEN(TRIM(ADJUSTL(char_block2)))
!        strlen3=LEN(TRIM(ADJUSTL(char_block3)))
!        char_block2=REPEAT('0',(3-strlen2))//TRIM(ADJUSTL(char_block2))
!        char_block3=REPEAT('0',(4-strlen2))//TRIM(ADJUSTL(char_block3))
!        gridfile='Bubble_'//TRIM(ADJUSTL(char_block2))//'_'//		!create a file Bubble_block_num.dat
!     &	TRIM(ADJUSTL(char_block3))//'.dat'
!         open (unit=2, file=gridfile)

!	endif

	R=0.5*dp_pt(M)
	maxnzr=0 ;   maxnode = 0

	nzr(M) = nint((2.d0*PI*R/2.d0)/dxm) !Number of planes
	thz(M) = PI/(nzr(M))
	maxnzr=max(maxnzr,nzr(M))

!	if ((mod(itime,tsteps_pt).eq.0).and.
!     &		(itime.gt.itime_start)) then 			    
!		write (2,*) 'variables="x","y","z"'
!	endif

	allocate (ztemp_layer(np,maxnzr),ctot_layer(np,maxnzr))
        allocate (nodes_layer(np,100000),Rtemp_layer(np,maxnzr,100000))
	allocate (nodes_percyl_layer(np,maxnzr,100000))

         nodes_percyl_layer = 0 ;    nodes_layer = 0 ;nodes_layer=0
                  
         nodes(M) = 0

	 do izr = 1,nzr(M)	 
	      c =1
         if (izr.eq.1) then
	   Rtemp_layer(M,izr,c) = 0.d0
         else
	   Rtemp_layer(M,izr,c) = R*cos(thz(M)*(izr-1)-(PI/2.d0))
	   ztemp_layer(M,izr) = R*sin(thz(M)*(izr-1)-(PI/2.d0))+zp_pt(M)
         end if

	      do while (Rtemp_layer(M,izr,c).ge.0.)		!gt!!!
                nodes_layer(M,izr) = nodes_layer(M,izr) + 
     &  NINT(2.0*PI*Rtemp_layer(M,izr,c)/dxm)
                nodes_percyl_layer(M,izr,c) = 
     &  NINT(2.0*PI*Rtemp_layer(M,izr,c)/dxm)
               if (Rtemp_layer(M,izr,c).eq.0.) then
                nodes_layer(M,izr) = nodes_layer(M,izr) + 1
                nodes_percyl_layer(M,izr,c) = 1
               end if
               c = c + 1
               Rtemp_layer(M,izr,c) = Rtemp_layer(M,izr,c-1)- dxm
               
             if(c.ge.cmax) goto 555 
               
	      end do

555	CONTINUE
              ctot_layer(M,izr) = c - 1
	      if(ctot_layer(M,izr) .gt. 100) then
	      print*, 'allocate problem in Rtemp_layer'
	      end if
          end do      	    

	 izr=1
	    K = 1
	 do while(izr.le.nzr(M))
	   do c = 1,ctot_layer(M,izr)
	     thc = 2.d0*PI/nodes_percyl_layer(M,izr,c)
	     do L = 1,nodes_percyl_layer(M,izr,c)
	     zccos = cos(2.d0*PI*(L-1)/nodes_percyl_layer(M,izr,c)+thc)
	     zcsin = sin(2.d0*PI*(L-1)/nodes_percyl_layer(M,izr,c)+thc)

             pre_nodex(M,K) = xp_pt(M) + Rtemp_layer(M,izr,c)*zccos
	     pre_nodey(M,K) = yp_pt(M) + Rtemp_layer(M,izr,c)*zcsin
	     pre_nodez(M,K) = zp_pt(M) + R*sin(thz(M)*(izr-1)-(PI/2.d0))

  !	if ((mod(itime,tsteps_pt).eq.0).and.
!     &		(itime.gt.itime_start)) then 			
! 	      write(2,89) pre_nodex(M,K), pre_nodey(M,K), pre_nodez(M,K)
!       end if  

               K = K + 1   
!THIS SHOULD PREVENT IBs OUTSIDE THE DOMAIN
	      if (pre_nodex(M,K).lt.(xst-pl*dxm)) K = K -1
	      if (pre_nodex(M,K).gt.(xen+pl*dxm)) K = K -1
	      if (pre_nodey(M,K).lt.(yst-pl*dym)) K = K -1
	      if (pre_nodey(M,K).gt.(yen+pl*dym)) K = K -1
	      if (pre_nodez(M,K).lt.(zst-pl*dzm)) K = K -1
	      if (pre_nodez(M,K).gt.(zen+pl*dzm)) K = K -1


	     end do
	    end do
            izr = izr + 1     
          end do
 	      nodes(M) = K-1
  	      maxnode = max(maxnode,nodes(M))	          

	deallocate (ztemp_layer,ctot_layer,nodes_layer)
	deallocate (nodes_percyl_layer,Rtemp_layer)

!	if ((mod(itime,tsteps_pt).eq.0).and.
!     &		(itime.gt.itime_start)) then 			

!       		close (2)
!       end if        

   88 FORMAT (i15)
   89 FORMAT (3e20.6)

      RETURN
      end subroutine 
!#############################################################
      SUBROUTINE init_IMBub
!#############################################################
      use vars
      use multidata
      use imbub
      use mpi
      implicit none
      INTEGER      :: L,I,M,maxn

!	allocate (bubsinproc(nprocs))
	allocate (imbsinproc(nprocs))
	allocate (imbinblk(num_domains),imb_block(maxnodeIBS))

	imbsinproc=0	;	imbinblk=0
!	bubsinproc=0		
	imb_block=0

	IF (myrank.eq.0) then

	WRITE(6,*)' '
	WRITE(6,*)'~~~~~~~~~~  Starting Immersed Boundaries  ~~~~~~~~~~'
	WRITE(6,*)' '

	 cmax=1	
      	 allocate (nodes(np))
	 nodes=0 	
	 maxn=np*800.d0	!maximum number of Lagrangian allowed

      	 allocate (pre_nodex(np,maxn))
	 allocate (pre_nodey(np,maxn))
	 allocate (pre_nodez(np,maxn))
        pre_nodex = 0.d0 ; pre_nodey = 0.d0 ; pre_nodez = 0.d0  
	maxnodeIBS=0
	
        dxm=g_dx/rdivmax ; dym=g_dy/rdivmax ; dzm=g_dz/rdivmax		!Minimum grid sizes
	 
!	write(6,*)'Largest rdivmax  :',rdivmax	
!	write(6,'(a,3e12.4)')'Smallest gridsize: ',dxm,dym,dzm	 
		
	Do M=1,np
 	    call imb_spheric_bub(M)
   	     maxnodeIBS=maxnodeIBS+nodes(M)
	  IF (maxnodeIBS.gt.maxn) write(6,*)'Too many ib points'
	  IF (maxnodeIBS.gt.maxn) STOP
	Enddo
        WRITE(6,*)'Number of IBs per bubble:',nodes(1)
	WRITE(6,*)'Total number of IBs:',maxnodeIBS
	WRITE(6,*)' '

      	 allocate (nodex(np,maxnodeIBS))
	 allocate (nodey(np,maxnodeIBS))
	 allocate (nodez(np,maxnodeIBS))
        nodex = 0.d0	; nodey = 0.d0 	  ; nodez = 0.d0  

	do M=1,np 	; 	do L=1,nodes(M)
		nodex(M,L)=pre_nodex(M,L)
		nodey(M,L)=pre_nodey(M,L)
		nodez(M,L)=pre_nodez(M,L)
	enddo		;	enddo

	allocate (lag_bod(maxnodeIBS))
!	allocate (bub_block(np))
!	allocate (bub_proc(np))!,bubinblk(num_domains))
	allocate (imb_proc(maxnodeIBS))

	lag_bod=0;	!bub_block=0
!	bubinblk=0	; bub_proc=0
	imb_proc=0	

      	 deallocate (pre_nodex)
	 deallocate (pre_nodey)
	 deallocate (pre_nodez)

	ENDIF

      RETURN
      end subroutine
!#############################################################
      SUBROUTINE update_IMBub
!#############################################################
      use vars
      use multidata
      use imbub
      use mpi
      implicit none
      INTEGER      :: L,I,M,maxn     
	
	 maxn=np*800.d0	!maximum number of Lagrangian allowed

      	 allocate (pre_nodex(np,maxn))
	 allocate (pre_nodey(np,maxn))
	 allocate (pre_nodez(np,maxn))
        pre_nodex = 0.d0 ; pre_nodey = 0.d0 ; pre_nodez = 0.d0  
	maxnodeIBS=0
	
	Do M=1,np
 	    call imb_spheric_bub(M)
   	     maxnodeIBS=maxnodeIBS+nodes(M)
	  IF (maxnodeIBS.gt.maxn) write(6,*)'Too many ib points'
	  IF (maxnodeIBS.gt.maxn) STOP
	Enddo

      	 deallocate (nodex)
	 deallocate (nodey)
	 deallocate (nodez)

      	 allocate (nodex(np,maxnodeIBS))
	 allocate (nodey(np,maxnodeIBS))
	 allocate (nodez(np,maxnodeIBS))
        nodex = 0.d0	; nodey = 0.d0 	  ; nodez = 0.d0  

	do M=1,np 	; 	do L=1,nodes(M)
		nodex(M,L)=pre_nodex(M,L)
		nodey(M,L)=pre_nodey(M,L)
		nodez(M,L)=pre_nodez(M,L)
	enddo		;	enddo


	deallocate (imb_block,lag_bod,imb_proc)

	allocate (imb_block(maxnodeIBS),lag_bod(maxnodeIBS))
	allocate (imb_proc(maxnodeIBS))

	imb_block=0	; lag_bod=0
	imb_proc=0	


      	 deallocate (pre_nodex)
	 deallocate (pre_nodey)
	 deallocate (pre_nodez)

      RETURN
      end subroutine
