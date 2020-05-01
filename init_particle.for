!=======================================================================
!    			LAGRANGIAN PARTICLE TRACKING
!			Bruño Fraga Bugallo
!			Cardiff 2013-2015
!=======================================================================
C#############################################################
      SUBROUTINE INIT_PARTICLE
C#############################################################
        use multidata
        use mpi
        use vars   
 	use imbub

	implicit none	

	integer l,ib
	logical dummy
      	double precision X,Y,Z
	double precision :: Dx,Dy,Dz,random_number_normal
      	integer :: i,j,k,tti,ttj,ttk

20  	format(' Domain dimensions (m):',3x,e14.5)
21  	format(' Average diameter (mm):',F5.2)

	if (myrank.eq.0) write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       open(10,file='IMBub.cin')
	if (myrank.eq.0) write(6,*)'IMBub model for DNS of bubbly flows'
     &,' with immersed boundaries'
	if (myrank.eq.0) write(6,*)'////\\\\Bruño Fraga Bugallo////\\\\'
	if (myrank.eq.0) write(6,*)'========CTR Stanford 2018=========='
	read(10,*) order
	if (myrank.eq.0) write(6,*)'Delta function:',order
       	read(10,*) tsteps_pt,tsnr
	 if (myrank.eq.0) write(6,*)'Transient output:',tsteps_pt!,
!     &'Release frequency:',100/tsnr,'%'
       	read(10,*) ptnr
	if (myrank.eq.0) write(6,*)'Number of bubbles released each time
     &:',ptnr
      	read(10,*) Dp,random_dp,sigma
	if (myrank.eq.0) write(6,21)Dp*1000
       	read(10,*) Dx,Dy,Dz
	if (myrank.eq.0) write(6,20)Dx,Dy,Dz
	read(10,*) random
!	 div = 1./(Dw)
	if (random) read(10,*)xp,yp,zp,uop,vop,wop
	if (myrank.eq.0) write(6,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	 np=ptnr

       if (LRESTART) then
       	open(20,file='final_particle.dat')
       	read(20,*) np
       	read(20,*) cnt_pt
       end if

       	allocate (uop_pt(np),vop_pt(np),wop_pt(np))
	allocate (Fpu(np),Fpv(np),Fpw(np))
	allocate (ui_pt(np),vi_pt(np),wi_pt(np))
	allocate (Fu_bub(np),Fv_bub(np),Fw_bub(np))   	
	allocate (U_bub(np),V_bub(np),W_bub(np))

	IF (myrank.eq.0) then
       	allocate (xp_pt(np),yp_pt(np),zp_pt(np))
	allocate (dp_pt(np))
        allocate (uoi_pt(np),voi_pt(np),woi_pt(np))
	if (tsnr.gt.0) then
       	allocate (xpold(np),ypold(np),zpold(np))
       	allocate (uopold(np),vopold(np),wopold(np))
       	allocate (uoiold(np),voiold(np),woiold(np))
	allocate (dp_old(np))
	endif

!	 allocate(ptsinproc(nprocs))

       if (LRESTART) then

       do l=1,np
       	read(20,*) xp_pt(l),yp_pt(l),zp_pt(l)
       	read(20,*) uop_pt(l),vop_pt(l),wop_pt(l)
		read(20,*) dp_pt(l)
		uoi_pt(l)=0.d0
		voi_pt(l)=0.d0
		woi_pt(l)=0.d0
		if (tsnr.gt.0) then
		xpold(l)=xp_pt(l)
		ypold(l)=yp_pt(l)
		zpold(l)=zp_pt(l)
		uopold(l)=uop_pt(l)
		vopold(l)=vop_pt(l)
		wopold(l)=wop_pt(l)	
		dp_old(l)=dp_pt(l)
		uoiold(l)=uoi_pt(l)
		voiold(l)=voi_pt(l)
		woiold(l)=woi_pt(l)
		endif
       end do
	close (20)

       elseif (LRESTART.eq..false.) then
	if (tsnr.gt.0) then					!continuous release
	if (myrank.eq.0) then
	write(202,*) 'First release:',np,'new particles. Total:',np
	endif
	       if (random) then
	       do l=1,np
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(Z)
			zp_pt(l)=Z*Dz				!4cm CTR Stanford
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(Y)
			yp_pt(l)=Y*Dy				!4cm CTR Stanford

			xp_pt(l)=xp	

			if (random_dp) then						!random dimeter of particle

			dp_pt(l)= random_number_normal(Dp,sigma)

			else

			dp_pt(l) = Dp

			endif
			uop_pt(l)=uop
			vop_pt(l)=vop
			wop_pt(l)=wop
			uoi_pt(l)=0.d0
			voi_pt(l)=0.d0
			woi_pt(l)=0.d0
			if (tsnr.gt.0) then
			xpold(l)=xp_pt(l)
			ypold(l)=yp_pt(l)
			zpold(l)=zp_pt(l)
			uopold(l)=uop_pt(l)
			vopold(l)=vop_pt(l)
			wopold(l)=wop_pt(l)
			dp_old(l)=dp_pt(l)
			uoiold(l)=uoi_pt(l)
			voiold(l)=voi_pt(l)
			woiold(l)=woi_pt(l)
			endif
!			write(202,*)l,X,xp_pt(l),yp_pt(l),zp_pt(l)
	       end do
		else								!not random
       		do l=1,np
      			read(10,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &				uop_pt(l),vop_pt(l),wop_pt(l)
				dp_pt(l) = Dp				!if release isnt random, then Dp is constant
				uoi_pt(l)=0.d0
				voi_pt(l)=0.d0
				woi_pt(l)=0.d0
				if (tsnr.gt.0) then
				xpold(l)=xp_pt(l)
				ypold(l)=yp_pt(l)
				zpold(l)=zp_pt(l)
				uopold(l)=uop_pt(l)
				vopold(l)=vop_pt(l)
				wopold(l)=wop_pt(l)
				dp_old(l)=dp_pt(l)
				uoiold(l)=uoi_pt(l)
				voiold(l)=voi_pt(l)
				woiold(l)=woi_pt(l)
				endif
      		 end do
		endif							!random
 !		call TECPARTICLE(0)
	else								!Periodic domain, no further release
		if (random) then
		do l=1,np
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(X)
			xp_pt(l)=X*Dx					! CTR Stanford
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(Y)
			yp_pt(l)=Y*Dy					! CTR Stanford
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(Z)
			zp_pt(l)=Z*Dz					! CTR Stanford	

			if (random_dp) then						!random dimeter of particle

			dp_pt(l)= random_number_normal(Dp,sigma)

			else

			dp_pt(l) = Dp

			endif
			uop_pt(l)=uop
			vop_pt(l)=vop
			wop_pt(l)=wop
			uoi_pt(l)=0.d0
			voi_pt(l)=0.d0
			woi_pt(l)=0.d0
			if (tsnr.gt.0) then
			xpold(l)=xp_pt(l)
			ypold(l)=yp_pt(l)
			zpold(l)=zp_pt(l)
			uopold(l)=uop_pt(l)
			vopold(l)=vop_pt(l)
			wopold(l)=wop_pt(l)
			dp_old(l)=dp_pt(l)
			uoiold(l)=uoi_pt(l)
			voiold(l)=voi_pt(l)
			woiold(l)=woi_pt(l)
			endif
		enddo
		else								!not random
       		do l=1,np
      			read(10,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &				uop_pt(l),vop_pt(l),wop_pt(l)
				dp_pt(l) = Dp				!if release isnt random, then Dp is constant
				uoi_pt(l)=0.d0
				voi_pt(l)=0.d0
				woi_pt(l)=0.d0
				if (tsnr.gt.0) then
				xpold(l)=xp_pt(l)
				ypold(l)=yp_pt(l)
				zpold(l)=zp_pt(l)
				uopold(l)=uop_pt(l)
				vopold(l)=vop_pt(l)
				wopold(l)=wop_pt(l)
				dp_old(l)=dp_pt(l)
				uoiold(l)=uoi_pt(l)
				voiold(l)=voi_pt(l)
				woiold(l)=woi_pt(l)
				endif
      		 end do
		endif		
	endif									!tsnr=0
      end if									!restart
      ENDIF									!myrank
      close (10)

      RETURN
      END SUBROUTINE

C **********************************************************************
      SUBROUTINE TECPLOT(num_output)
C **********************************************************************

        use multidata
        use mpi
        use vars   
	use imbub

	implicit none	

	integer strlen,i,j,k,ib,ni,nj,nk,ii,idfile,num_output
	integer is,ie,js,je,ks,ke
      	character(LEN=19) filename
      	character(LEN=4) b_str
      	character(LEN=3) c_str
	double precision u_cn,v_cn,w_cn,p_cn,S_cn!,T_cn,k_cn,eps_cn,vis_cn
   		
!        if (LRESTART) KK1=KK+KK2
!        if (LRESTART.eq..false.) KK1=KK

	do ib=1,nbp

!	if ((dom_id(ib).ge.20.and.dom_id(ib).le.29).or.(dom_id(ib).ge.70
!     &	.and.dom_id(ib).le.79).or.(dom_id(ib).ge.120.and.dom_id(ib).le.
!     &  129).or.(dom_id(ib).ge.170.and.dom_id(ib).le.179).or.(dom_id(ib)
!     &  .ge.220.and.dom_id(ib).le.229)) 	then

	  idfile=600+dom_id(ib)

        write(b_str,'(I4)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(4-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"
        write(c_str,'(I3)') dom_id(ib)
        strlen=LEN(TRIM(ADJUSTL(c_str)))
        c_str=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(c_str)) ! e.g. "001"

      filename='tecout_'//b_str//'_'//c_str//'.dat'

      OPEN (UNIT=idfile, FILE=filename)

      WRITE (idfile,*) 'TITLE = ', '"Eulerian field"'
      WRITE (idfile,"(A)")'VARIABLES = "X","Y","Z","U","V","W","P"'!,"S",
!     &"k","eps","vis"'

        is=pl+1; ie=dom(ib)%ttc_i-pl+1
        js=pl+1; je=dom(ib)%ttc_j-pl+1
        ks=pl+1; ke=dom(ib)%ttc_k-pl+1
        ni=ie-(is-1)+1
        nj=je-(js-1)+1
        nk=ke-(ks-1)+1


      !WRITE(idfile,*)'ZONE T="','id:',dom_id(ib),'it:',ntime,'"'
      WRITE(idfile,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime
	WRITE(idfile,*)'I=',0.5*ni,', J=',0.5*nj,', K=',0.5*nk,'F=POINT'

		do k=ks-1,ke,2
		do j=js-1,je,2
		do i=is-1,ie,2

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))

                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 

                 p_cn  =0.125*(dom(ib)%p(i,j,k)+
     &dom(ib)%p(i+1,j,k)    +dom(ib)%p(i,j+1,k)+
     &dom(ib)%p(i+1,j+1,k)  +dom(ib)%p(i,j,k+1)+
     &dom(ib)%p(i+1,j,k+1)  +dom(ib)%p(i,j+1,k+1)+
     &dom(ib)%p(i+1,j+1,k+1))
!                 S_cn  =0.125*(dom(ib)%S(i,j,k)+
!     &dom(ib)%S(i+1,j,k)    +dom(ib)%S(i,j+1,k)+
!     &dom(ib)%S(i+1,j+1,k)  +dom(ib)%S(i,j,k+1)+
!     &dom(ib)%S(i+1,j,k+1)  +dom(ib)%S(i,j+1,k+1)+
!     &dom(ib)%S(i+1,j+1,k+1))
!                 k_cn  =0.125*(dom(ib)%ksgs(i,j,k)+
!     &dom(ib)%ksgs(i+1,j,k)    +dom(ib)%ksgs(i,j+1,k)+
!     &dom(ib)%ksgs(i+1,j+1,k)  +dom(ib)%ksgs(i,j,k+1)+
!     &dom(ib)%ksgs(i+1,j,k+1)  +dom(ib)%ksgs(i,j+1,k+1)+
!     &dom(ib)%ksgs(i+1,j+1,k+1))
!                 eps_cn  =0.125*(dom(ib)%eps(i,j,k)+
!     &dom(ib)%eps(i+1,j,k)    +dom(ib)%eps(i,j+1,k)+
!     &dom(ib)%eps(i+1,j+1,k)  +dom(ib)%eps(i,j,k+1)+
!     &dom(ib)%eps(i+1,j,k+1)  +dom(ib)%eps(i,j+1,k+1)+
!     &dom(ib)%eps(i+1,j+1,k+1))
!                 vis_cn  =0.125*(dom(ib)%vis(i,j,k)+
!     &dom(ib)%vis(i+1,j,k)    +dom(ib)%vis(i,j+1,k)+
!     &dom(ib)%vis(i+1,j+1,k)  +dom(ib)%vis(i,j,k+1)+
!     &dom(ib)%vis(i+1,j,k+1)  +dom(ib)%vis(i,j+1,k+1)+
!     &dom(ib)%vis(i+1,j+1,k+1))
!                 T_cn  =0.125*(dom(ib)%T(i,j,k)+
!     &dom(ib)%T(i+1,j,k)    +dom(ib)%T(i,j+1,k)+
!     &dom(ib)%T(i+1,j+1,k)  +dom(ib)%T(i,j,k+1)+
!     &dom(ib)%T(i+1,j,k+1)  +dom(ib)%T(i,j+1,k+1)+
!     &dom(ib)%T(i+1,j+1,k+1))

      write (idfile,'(11e14.6)') dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k)
     & ,u_cn,v_cn,w_cn,p_cn!,S_cn!T_cn,k_cn,eps_cn,vis_cn

	enddo
	enddo
	enddo
!		write (90,*) dom(ib)%isp,dom(ib)%iep,
!     & dom(ib)%jsp,dom(ib)%jep,dom(ib)%ksp,dom(ib)%kep

!	endif
	end do



      close (idfile)

!   88 FORMAT (10F15.8)

      END SUBROUTINE


C **********************************************************************
      SUBROUTINE TECPARTICLE(num_output)
C **********************************************************************
C
  
        use multidata
        use mpi
        use vars   
   	use imbub

	implicit none	

      integer l,strlen,ib,num_output
      character(LEN=80) filename,filename2
      character(LEN=4) b_str

        write(b_str,'(I4)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(4-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"

        filename='tecout_'//b_str//'_pt.dat'

        OPEN (UNIT=95, FILE=TRIM(ADJUSTL(filename)))

      WRITE (95,*) 'TITLE = ', '"Lagrangian field"'
      WRITE (95,"(A)")'VARIABLES = "X","Y","Z","U<sub>Lag<\sub>","V<sub>
     &Lag<\sub>","W<sub>Lag<\sub>"'!,"D<sub>Lag<\sub>"'
      WRITE(95,*)'zone ','STRANDID=', 2, 'SOLUTIONTIME=', ctime

        do l=1,np
      	WRITE (95,*) xp_pt(l),yp_pt(l),zp_pt(l)
     &			,uop_pt(l),vop_pt(l),wop_pt(l)
!     &			,dp_pt(l)
        end do
	  close (95)

!   88 FORMAT (10F15.8)

      END SUBROUTINE


	
