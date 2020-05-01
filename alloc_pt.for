!!=======================================================================
!    			LAGRANGIAN PARTICLE TRACKING
!			Bruño Fraga Bugallo
!			Cardiff 2013-2015
!			Birmingham-Stanford 2018
!=======================================================================
!##########################################################################
        subroutine alloc_pt
!##########################################################################
        use vars
        use mpi
        use multidata
	 use imbub
        implicit none

        integer l,out_cnt   
	  logical,allocatable,dimension(:):: out_pt			

!	wtime_refresh = MPI_WTIME ( ) 

	allocate(out_pt(np))

	out_pt=.false.

	out_cnt = 0

	do l=1,np

!	Comprobar si permanece en el dominio
      	if ((xp_pt(l).le.xst).or.(xp_pt(l).ge.xen)) then
       		out_pt(l) = .TRUE.
	 	out_cnt = out_cnt + 1
	 	goto 50					
      	elseif ((yp_pt(l).le.yst).or.(yp_pt(l)
     &.ge.yen)) then
       		out_pt(l) = .TRUE.
	 	out_cnt = out_cnt + 1
	 	goto 50
      	elseif ((zp_pt(l).lt.zst).or.(zp_pt(l)
     &.ge.zen)) then
       		out_pt(l) = .TRUE.
	 	out_cnt = out_cnt + 1
	 	goto 50
      	end if

50	continue

		if (out_pt(l).eq..false.) then
      			xpold(l-out_cnt)=xp_pt(l)
           		ypold(l-out_cnt)=yp_pt(l)
           		zpold(l-out_cnt)=zp_pt(l)
           		uopold(l-out_cnt)=uop_pt(l)
           		vopold(l-out_cnt)=vop_pt(l)
           		wopold(l-out_cnt)=wop_pt(l)

			dp_old(l-out_cnt)=dp_pt(l)

			uoiold(l-out_cnt)=uoi_pt(l)
           		voiold(l-out_cnt)=voi_pt(l)
           		woiold(l-out_cnt)=woi_pt(l)
		endif

      	enddo
!	enddo

	np = np - out_cnt

	if (np.eq.0) then
		if (myrank.eq.0) then
		write(202,*)'Time:',ntime,'No particles in the domain'
		endif
		write(6,*)'No bubbles left. EXIT'
		Stop
	endif

	if (out_cnt.gt.0) then
		if (myrank.eq.0) then
		write(202,*)'Time:',ntime,'Removing',out_cnt,'particles'
     &,'. Total:',np
		endif
	endif
      		deallocate (xp_pt,yp_pt,zp_pt)
		deallocate (Fpu,Fpv,Fpw)
	 	deallocate (ui_pt,vi_pt,wi_pt)
      		deallocate (uop_pt,vop_pt,wop_pt)

         	allocate (xp_pt(np),yp_pt(np),zp_pt(np))
         	allocate (uop_pt(np),vop_pt(np),wop_pt(np))
		allocate (Fpu(np),Fpv(np),Fpw(np))
	 	allocate (ui_pt(np),vi_pt(np),wi_pt(np))

		deallocate (out_pt)
      		deallocate (uoi_pt,voi_pt,woi_pt)
		deallocate (dp_pt)

         	allocate (uoi_pt(np),voi_pt(np),woi_pt(np))
		allocate (dp_pt(np))

	do l=1,np
			xp_pt(l) = xpold(l)
           		yp_pt(l) = ypold(l)
           		zp_pt(l) = zpold(l)
           		uop_pt(l)= uopold(l)
           		vop_pt(l)= vopold(l)
           		wop_pt(l)= wopold(l)
           		uoi_pt(l)= uoiold(l)
           		voi_pt(l)= voiold(l)
           		woi_pt(l)= woiold(l)
			dp_pt(l) = dp_old(l)
	enddo

	deallocate (xpold,ypold,zpold,uopold,vopold,wopold,dp_old)
	deallocate (uoiold,voiold,woiold)

        allocate(xpold(np),ypold(np),zpold(np))
    	allocate(uopold(np),vopold(np),wopold(np)) 
    	allocate(uoiold(np),voiold(np),woiold(np)) 
	allocate(dp_old(np))

	do l=1,np
			xpold(l)=xp_pt(l)
           		ypold(l)=yp_pt(l)
           		zpold(l)=zp_pt(l)
           		uopold(l)=uop_pt(l)
           		vopold(l)=vop_pt(l)
           		wopold(l)=wop_pt(l)
           		uoiold(l)=uoi_pt(l)
           		voiold(l)=voi_pt(l)
           		woiold(l)=woi_pt(l)
			dp_old(l)=dp_pt(l)
	enddo		

      return
      end subroutine


!##########################################################################
        subroutine release_pt
!##########################################################################
        use vars
        use mpi
        use multidata
 	use imbub
        implicit none

        integer l
        double precision X,Y	
	double precision random_number_normal

       	np = np + ptnr
		if (myrank.eq.0) then
	  	write(202,*) ntime,'Releasing',ptnr,'new particles.', 
     &			' Total:',np
		endif

      	deallocate (xp_pt,yp_pt,zp_pt)
      	deallocate (uop_pt,vop_pt,wop_pt)
      	deallocate (uoi_pt,voi_pt,woi_pt)
	deallocate (dp_pt)
	deallocate (Fpu,Fpv,Fpw)
	 deallocate (ui_pt,vi_pt,wi_pt)

         allocate (xp_pt(np),yp_pt(np),zp_pt(np))
         allocate (uop_pt(np),vop_pt(np),wop_pt(np))
         allocate (uoi_pt(np),voi_pt(np),woi_pt(np))
	allocate (dp_pt(np))
	allocate (Fpu(np),Fpv(np),Fpw(np))
	 allocate (ui_pt(np),vi_pt(np),wi_pt(np))

         do l=1,(np-ptnr)

         	xp_pt(l)=xpold(l)
         	yp_pt(l)=ypold(l)
         	zp_pt(l)=zpold(l)

         	uop_pt(l)=uopold(l)
         	vop_pt(l)=vopold(l)
         	wop_pt(l)=wopold(l)    

		uoi_pt(l)=0.d0
		voi_pt(l)=0.d0
		woi_pt(l)=0.d0
		uoiold(l)=uoi_pt(l)
		voiold(l)=voi_pt(l)
		woiold(l)=woi_pt(l)

		dp_pt(l)=dp_old(l)

         end do

         deallocate (xpold,ypold,zpold)
         allocate(xpold(np),ypold(np),zpold(np))
         deallocate (uopold,vopold,wopold)
         allocate (uopold(np),vopold(np),wopold(np))
         deallocate (uoiold,voiold,woiold)
         allocate (uoiold(np),voiold(np),woiold(np))
         deallocate (dp_old)
         allocate (dp_old(np))

	   if (random) then

!         	read(35,*)xp,yp,zp,uop,vop,wop

!	write(401,*)'A Release',
!     &	np-ptnr+1,xp_pt(np-ptnr+1),yp_pt(np-ptnr+1),
!     &	zp_pt(np-ptnr+1),out_pt(np-ptnr+1),' '

         	do l=(np-ptnr+1),np
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(X)
			zp_pt(l)=X/div				!10cm CTR Stanford
			CALL RANDOM_SEED
			CALL RANDOM_NUMBER(Y)
			yp_pt(l)=X/div				!10cm CTR Stanford

			xp_pt(l)=zp	
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
			uoiold(l)=uoi_pt(l)
			voiold(l)=voi_pt(l)
			woiold(l)=woi_pt(l)

         	end do
	
		else										!if release isnt random, then Dp is constant
		open(35,file='LPT.cin')
       		read(35,*)
       		read(35,*)
       		read(35,*) 
       		read(35,*) 
       		read(35,*) 
		read(35,*) 
       		do l=np-ptnr+1,np
      			read(35,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &				uop_pt(l),vop_pt(l),wop_pt(l)
				xpold(l)=xp_pt(l)
				ypold(l)=yp_pt(l)
				zpold(l)=zp_pt(l)
				uopold(l)=uop_pt(l)
				vopold(l)=vop_pt(l)
				wopold(l)=wop_pt(l)
				dp_pt(l) = Dp

				uoi_pt(l)=0.d0
				voi_pt(l)=0.d0
				woi_pt(l)=0.d0
				uoiold(l)=uoi_pt(l)
				voiold(l)=voi_pt(l)
				woiold(l)=woi_pt(l)
      		 end do
          close(35)

		endif	


      return
      end subroutine

!##########################################################################
        subroutine periodic_pt
!##########################################################################
        use vars
        use mpi
        use multidata
	use imbub
        implicit none

        integer l,out_cnt   

	do l=1,np

!	Comprobar si permanece en el dominio
      	if (xp_pt(l).lt.xst) then
		xp_pt(l)=xp_pt(l)+(xen-xst)				!bubble comes back at the top (unlikely)
	elseif(xp_pt(l).gt.xen) then
		xp_pt(l)=xp_pt(l)-(xen-xst)				!bubble comes back at the bottom (likely)
	endif
      	if (yp_pt(l).lt.yst) then
		yp_pt(l)=yp_pt(l)+(yen-yst)				!bubble comes back
	elseif (yp_pt(l).gt.yen) then
		yp_pt(l)=yp_pt(l)-(yen-yst)				!bubble comes back
	endif
      	if (zp_pt(l).lt.zst) then
		zp_pt(l)=zp_pt(l)+(zen-zst)				!bubble comes back
	elseif (zp_pt(l).gt.zen) then
		zp_pt(l)=zp_pt(l)-(zen-zst)				!bubble comes back
      	end if

	enddo


      return
      end subroutine
