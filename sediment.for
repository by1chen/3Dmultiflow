!=======================================================================
!    			Scalar transport
!			Yan Liu
!			Bruño Fraga Bugallo
!			Cardiff 2015-2016
!=======================================================================
!##########################################################################
        subroutine sediment_init
!##########################################################################

	use vars
	use mpi
	use multidata
	implicit none
	integer :: ib,i,j,k,tti,ttj,ttk
        integer :: is,ie,js,je,ks,ke
	double precision :: Vcell, temp1,temp2

	do ib=1,nbp

           tti=dom(ib)%ttc_i
           ttj=dom(ib)%ttc_j
           ttk=dom(ib)%ttc_k

		do k=1,ttk  !ks-1,ke
		   do j=1,ttj !js-1,je
		      do i=1,tti !is-1,ie

		if (dom(ib)%z(k).gt.0.625) then  !========= 1st layer
                        dom(ib)%Sp(i,j,k) = 0.0 
                        dom(ib)%S(i,j,k) = 0.435 
                        temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
			dom(ib)%dens(i,j,k) = temp2*1000
!		elseif (dom(ib)%z(k).gt.0.480.and.dom(ib)%z(k).le.0.600) then !========= 2nd layer
!			dom(ib)%Sp(i,j,k) = 0.0
!			dom(ib)%S(i,j,k) = 1.358
!			temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
!			dom(ib)%dens(i,j,k) = temp2*1000
!		elseif (dom(ib)%z(k).gt.0.360.and.dom(ib)%z(k).le.0.480) then !=======3rd layer 
!			dom(ib)%S(i,j,k) = 2.017   
!			dom(ib)%Sp(i,j,k) = 0.0  
!			temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
!			dom(ib)%dens(i,j,k) = temp2*1000
!		elseif (dom(ib)%z(k).gt.0.240.and.dom(ib)%z(k).le.0.360) then !=======4th layer
!			dom(ib)%Sp(i,j,k) = 0.0  
!			dom(ib)%S(i,j,k) = 2.676 
!			temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
!			dom(ib)%dens(i,j,k) = temp2*1000
!		elseif (dom(ib)%z(k).gt.0.120.and.dom(ib)%z(k).le.0.240) then ! =======5th layer
!                        dom(ib)%Sp(i,j,k) = 0.0  
!                        dom(ib)%S(i,j,k) = 3.33 
!                        temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
!                        dom(ib)%dens(i,j,k) = temp2*1000
		else							      ! ========6th layer
                        dom(ib)%Sp(i,j,k) = 0.0 
                        temp1 = dom(ib)%z(k)        ! correct the unit of depth cm
                        dom(ib)%S(i,j,k) = -5.69*temp1 + 3.994 
!                        dom(ib)%S(i,j,k) = 3.994 
                        temp2 = (0.007587*dom(ib)%S(i,j,k)+0.9947)
                        dom(ib)%dens(i,j,k) = temp2*1000
		endif
			  end do
		      end do
		    end do


!		if (dom_id(ib).lt.60) then 
!                        dom(ib)%S = 3.994 
!		elseif (dom_id(ib).lt.120) then 
!                        dom(ib)%S = 3.33 
!		elseif (dom_id(ib).lt.180) then 
!                        dom(ib)%S = 3.33 
!		elseif (dom_id(ib).lt.240) then 
!                        dom(ib)%S = 2.017
!		else
!                        dom(ib)%S = 1.358 
!		if (dom_id(ib).eq.18) then
!		Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	
!		dom(ib)%S(dom(ib)%isp,dom(ib)%jep,4)=
!     & dom(ib)%S(dom(ib)%isp,dom(ib)%jep,4) + (1d-10/Vcell)*dt
!		endif
	enddo

	end

!##########################################################################
        subroutine sediment_4thtest
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib,tti,ttj,ttk
        double precision :: dxx,dyy,dzz
        double precision :: conv,diff
        double precision :: duSdx,dvSdy,dwSdz,dwsSdz,dSdt
        double precision :: awS,aeS,asS,anS,ab_S,atS,apS
	double precision :: masout
	double precision :: ues,uws,vns,vss,wts,wbs,wsbs
	double precision :: masconv,masdiff,masws
	double precision :: kp,km,ku,kc,kd,b_r,ws,Vcell

        do ib=1,nbp

 	  Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	

!	  if (dom_id(ib).eq.9) then							!release point for bubble plume
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'1',dom(ib)%S(dom(ib)%iep,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.10) then
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'2',dom(ib)%S(dom(ib)%isp,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.17) then
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'3',dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.18) then
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'4',dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
!        endif


	if (itime .eq. itime_start) then
	dom(ib)%sfactor = 1.0
	end if 

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
!-------Convection
	if(dom(ib)%u(i-1,j,k).gt.0.0) then
	ku=dom(ib)%So(i-2,j,k)
	kc=dom(ib)%So(i-1,j,k)
	kd=dom(ib)%So(i,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else if(dom(ib)%u(i-1,j,k).lt.0.0) then
	ku=dom(ib)%So(i+1,j,k)
	kc=dom(ib)%So(i,j,k)
	kd=dom(ib)%So(i-1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else
	km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i-1,j,k))
	end if

	if(dom(ib)%u(i,j,k).gt.0.0) then
	ku=dom(ib)%So(i-1,j,k)
	kc=dom(ib)%So(i,j,k)
	kd=dom(ib)%So(i+1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	kp=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%u(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i+2,j,k)
           kc=dom(ib)%So(i+1,j,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i+1,j,k))
        end if
        duSdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!------
        if(dom(ib)%v(i,j-1,k).gt.0.0) then
           ku=dom(ib)%So(i,j-2,k)
           kc=dom(ib)%So(i,j-1,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%v(i,j-1,k).lt.0.0) then
           ku=dom(ib)%So(i,j+1,k)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j-1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j-1,k))
        end if

        if(dom(ib)%v(i,j,k).gt.0.0) then
           ku=dom(ib)%So(i,j-1,k)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j+1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%v(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i,j+2,k)
           kc=dom(ib)%So(i,j+1,k)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j+1,k))
        end if
        dvSdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!-------
        if(dom(ib)%w(i,j,k-1).gt.0.0) then
           ku=dom(ib)%So(i,j,k-2)
           kc=dom(ib)%So(i,j,k-1)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k-1).lt.0.0) then
           ku=dom(ib)%So(i,j,k+1)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j,k-1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k-1))
        end if
        if(dom(ib)%w(i,j,k).gt.0.0) then
           ku=dom(ib)%So(i,j,k-1)
           kc=dom(ib)%So(i,j,k)
           kd=dom(ib)%So(i,j,k+1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k).lt.0.0) then
           ku=dom(ib)%So(i,j,k+2)
           kc=dom(ib)%So(i,j,k+1)
           kd=dom(ib)%So(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k+1))
        end if
        dwSdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz

        conv=(duSdx+dvSdy+dwSdz)

!-------Diffusion
                 awS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                 aeS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                 anS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                 asS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                 atS=-dom(ib)%vis(i,j,k)/(dzz*Pr)
                 ab_S=-dom(ib)%vis(i,j,k)/(dzz*Pr)

                 apS = -1.0*(awS+aeS+asS+anS+ab_S+atS)

        	diff=(apS*dom(ib)%So(i,j,k)+
     & anS*dom(ib)%So(i,j+1,k) + asS*dom(ib)%So(i,j-1,k)+
     & aeS*dom(ib)%So(i+1,j,k) + awS*dom(ib)%So(i-1,j,k)+
     & atS*dom(ib)%So(i,j,k+1) + ab_S*dom(ib)%So(i,j,k-1))

		dom(ib)%S(i,j,k)=dom(ib)%So(i,j,k)-dt*(conv+diff)

!		dom(ib)%dens(i,j,k)=(0.007587*dom(ib)%S(i,j,k)+0.9947)*1000.0
!	if (dom(ib)%S(i,j,k) .lt. 0.0) then
!	write (81,*) dom(ib)%S(i,j,k), dom(ib)%So(i,j,k)
!	write (81,*) i,j,k
!	write (81,*) conv,diff
!	write (81,*) duSdx,dvSdy,dwSdz
!	end if

	end do
	end do
	end do



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Testing the mpi block 10/2019~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	! correct the index for the scalar 

	! calculate the density 
	dom(ib)%dens=(0.007587*dom(ib)%S+0.9947)*1000.0
	

        end do

	do ib =1, nbp
           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
	if (dom(ib)%S(i,j,k) .gt. 100) then
!	call tecplot_S(itime)
	write(6,*)'ERROR: scalar too big'
	stop
	end if
	end do
	end do
	end do
	end do

 	call exchange(8)
       	call exchange(10)

        call boundS

        return
        end subroutine sediment_4thtest
!##########################################################################
        subroutine boundS
!##########################################################################
	use mpi
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ni,nj,nk,ly
        integer :: is,ie,js,je,ks,ke
	double precision :: absz,absy

        if (PERIODIC) call exchange_bc(8,pl_ex)

        do ly=0,pl_ex

        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for S
!..............................................................................
!=== West ===>   ..  4=wall  ..    1=Inflow
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%Tbc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(is-1-ly,j,k)= dom(ib)%S(is+ly,j,k)
              end do; end do

           else if (dom(ib)%Tbc_west.eq.1) then					!CHANGE
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(is-1-ly,j,k)= 1.0
              end do; end do
	   else if (dom(ib)%Tbc_west.eq. 11) then
	      do k=ks-1,ke+1; do j=js-1,je+1
	      absz=abs(dom(ib)%zc(k)-5.22)
	      absy=abs(dom(ib)%yc(j)-4.32)
	      if (absz .le. 0.5*dom(ib)%dz) then 
	      	if (absy .le. 0.5*dom(ib)%dy) then
		   dom(ib)%S(is-1-ly,j,k)= 1.0
	      	else 
	           dom(ib)%S(is-1-ly,j,k)= 0.0
	      	end if
	      else 
	         dom(ib)%S(is-1-ly,j,k)= 0.0
	      end if	
	      end do; end do
           end if
        end if
!...............................................................................
!=== East ===>   ..  4=wall  ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%Tbc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)	
              end do; end do

           else if (dom(ib)%Tbc_east.eq.2) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)
              end do; end do
           end if
        end if
!...............................................................................
!=== South ===>  ..  4=wall  ..      
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%Tbc_south.eq.4) then 
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%S(i,js-1-ly,k)= dom(ib)%S(i,js+ly,k)	
              end do; end do

           end if
        end if
!.............................................................................
!=== North ===>  ..  4=wall  ..    
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%Tbc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%S(i,je+1+ly,k) = dom(ib)%S(i,je-ly,k) 	
              end do; end do

           end if
        end if
!...............................................................................
!=== Bottom ===>  ..  6=Net deposition  ..   7=Erosion
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%Tbc_bottom.eq.6) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ks-1-ly)= 0.0	
              end do; end do

           else if (dom(ib)%Tbc_bottom.eq.7) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ks-1-ly)= 1.5 * dom(ib)%S(i,j,ks+ly)	
              end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>  ..  8=free surface
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%Tbc_top.eq.8) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%S(i,j,ke+1+ly) = dom(ib)%S(i,j,ke-ly)	
              end do; end do

           end if
        end if

!==============================================================================
        end do

	end do ! ly

        end subroutine boundS
!#############################################################################
       subroutine sediment_4thtest_passive
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
	integer :: tti,ttj,ttk
        double precision :: dxx,dyy,dzz,disx,disy,disz,dis
        double precision :: conv,diff
        double precision :: duSdx,dvSdy,dwSdz,dwsSdz,dSdt
        double precision :: awS,aeS,asS,anS,ab_S,atS,apS
	double precision :: masout
	double precision :: ues,uws,vns,vss,wts,wbs,wsbs
	double precision :: masconv,masdiff,masws
	double precision :: kp,km,ku,kc,kd,b_r,ws,Vcell

        do ib=1,nbp

 	  Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	
	 
	tti=dom(ib)%ttc_i
	ttj=dom(ib)%ttc_j
	ttk=dom(ib)%ttc_k


        do k=1,ttk
           do j=1,ttj
              do i=1,tti
		disx = (dom(ib)%x(i)-1.76)**2
		disy = (dom(ib)%y(j)-0.45)**2
		disz = (dom(ib)%z(k)-0.095)**2
!		if (dom(ib)%z(k).ge.0.005.and.dom(ib)%z(k).le.0.01) then
		dis=sqrt(disx+disy+disz)
		if (dis.le.0.01) then 
			dom(ib)%Spo(i,j,k)=dom(ib)%Spo(i,j,k)+1.0*dt
		endif 
!		endif 	
	       end do
	   end do
	end do

!	  if (dom_id(ib).eq.9) then							!release point for bubble plume
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'1',dom(ib)%S(dom(ib)%iep,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.10) then
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'2',dom(ib)%S(dom(ib)%isp,dom(ib)%jep,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jep,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.17) then
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31)=
!     & 	dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'3',dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
!	  elseif (dom_id(ib).eq.18) then
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!	  	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31)=
!     & 	dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
!		write(6,*)'4',dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,29)
!		write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,28)
!		write(6,*)Vcell,dt
!        endif


	if (itime .eq. itime_start) then
	dom(ib)%sfactor = 1.0
	end if 

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
!-------Convection
	if(dom(ib)%u(i-1,j,k).gt.0.0) then
	ku=dom(ib)%Spo(i-2,j,k)
	kc=dom(ib)%Spo(i-1,j,k)
	kd=dom(ib)%Spo(i,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else if(dom(ib)%u(i-1,j,k).lt.0.0) then
	ku=dom(ib)%Spo(i+1,j,k)
	kc=dom(ib)%Spo(i,j,k)
	kd=dom(ib)%Spo(i-1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	km=(kc+0.5*b_r*(kc-ku))
	else
	km=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i-1,j,k))
	end if

	if(dom(ib)%u(i,j,k).gt.0.0) then
	ku=dom(ib)%Spo(i-1,j,k)
	kc=dom(ib)%Spo(i,j,k)
	kd=dom(ib)%Spo(i+1,j,k)
	b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
	kp=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%u(i,j,k).lt.0.0) then
           ku=dom(ib)%Spo(i+2,j,k)
           kc=dom(ib)%Spo(i+1,j,k)
           kd=dom(ib)%Spo(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i+1,j,k))
        end if
        duSdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!------
        if(dom(ib)%v(i,j-1,k).gt.0.0) then
           ku=dom(ib)%Spo(i,j-2,k)
           kc=dom(ib)%Spo(i,j-1,k)
           kd=dom(ib)%Spo(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else if(dom(ib)%v(i,j-1,k).lt.0.0) then
           ku=dom(ib)%Spo(i,j+1,k)
           kc=dom(ib)%Spo(i,j,k)
           kd=dom(ib)%Spo(i,j-1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i,j-1,k))
        end if

        if(dom(ib)%v(i,j,k).gt.0.0) then
           ku=dom(ib)%Spo(i,j-1,k)
           kc=dom(ib)%Spo(i,j,k)
           kd=dom(ib)%Spo(i,j+1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%v(i,j,k).lt.0.0) then
           ku=dom(ib)%Spo(i,j+2,k)
           kc=dom(ib)%Spo(i,j+1,k)
           kd=dom(ib)%Spo(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i,j+1,k))
        end if
        dvSdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!-------
        if(dom(ib)%w(i,j,k-1).gt.0.0) then
           ku=dom(ib)%Spo(i,j,k-2)
           kc=dom(ib)%Spo(i,j,k-1)
           kd=dom(ib)%Spo(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k-1).lt.0.0) then
           ku=dom(ib)%Spo(i,j,k+1)
           kc=dom(ib)%Spo(i,j,k)
           kd=dom(ib)%Spo(i,j,k-1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=(kc+0.5*b_r*(kc-ku))
        else
           km=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i,j,k-1))
        end if
        if(dom(ib)%w(i,j,k).gt.0.0) then
           ku=dom(ib)%Spo(i,j,k-1)
           kc=dom(ib)%Spo(i,j,k)
           kd=dom(ib)%Spo(i,j,k+1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))

        else if(dom(ib)%w(i,j,k).lt.0.0) then
           ku=dom(ib)%Spo(i,j,k+2)
           kc=dom(ib)%Spo(i,j,k+1)
           kd=dom(ib)%Spo(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=(kc+0.5*b_r*(kc-ku))
        else
           kp=0.5*(dom(ib)%Spo(i,j,k)+dom(ib)%Spo(i,j,k+1))
        end if
        dwSdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz

        conv=(duSdx+dvSdy+dwSdz)

!-------Diffusion
                 awS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                 aeS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                 anS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                 asS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                 atS=-dom(ib)%vis(i,j,k)/(dzz*Pr)
                 ab_S=-dom(ib)%vis(i,j,k)/(dzz*Pr)

                 apS = -1.0*(awS+aeS+asS+anS+ab_S+atS)

        	diff=(apS*dom(ib)%Spo(i,j,k)+
     & anS*dom(ib)%Spo(i,j+1,k) + asS*dom(ib)%Spo(i,j-1,k)+
     & aeS*dom(ib)%Spo(i+1,j,k) + awS*dom(ib)%Spo(i-1,j,k)+
     & atS*dom(ib)%Spo(i,j,k+1) + ab_S*dom(ib)%Spo(i,j,k-1))

		dom(ib)%Sp(i,j,k)=dom(ib)%Spo(i,j,k)-dt*(conv+diff)

!		dom(ib)%dens(i,j,k)=(0.007587*dom(ib)%S(i,j,k)+0.9947)*1000.0
!	if (dom(ib)%S(i,j,k) .lt. 0.0) then
!	write (81,*) dom(ib)%S(i,j,k), dom(ib)%So(i,j,k)
!	write (81,*) i,j,k
!	write (81,*) conv,diff
!	write (81,*) duSdx,dvSdy,dwSdz
!	end if

	end do
	end do
	end do

        end do

	do ib =1, nbp
           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep
!	if (dom(ib)%Sp(i,j,k) .gt. 100) then
!	call tecplot_S(itime)
!	write(6,*)'ERROR: scalar too big'
!	stop
!	end if
	end do
	end do
	end do
	end do

        call exchange(20)

        call boundSp

        return
        end subroutine sediment_4thtest_passive
!##########################################################################
        subroutine boundSp
!##########################################################################
	use mpi
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ni,nj,nk,ly
        integer :: is,ie,js,je,ks,ke
	double precision :: absz,absy

        if (PERIODIC) call exchange_bc(8,pl_ex)

        do ly=0,pl_ex

        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for S
!..............................................................................
!=== West ===>   ..  4=wall  ..    1=Inflow
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%Tbc_west.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%Sp(is-1-ly,j,k)= dom(ib)%Sp(is+ly,j,k)
              end do; end do

           else if (dom(ib)%Tbc_west.eq.1) then					!CHANGE
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%Sp(is-1-ly,j,k)= 1.0
              end do; end do
	   else if (dom(ib)%Tbc_west.eq. 11) then
	      do k=ks-1,ke+1; do j=js-1,je+1
	      absz=abs(dom(ib)%zc(k)-5.22)
	      absy=abs(dom(ib)%yc(j)-4.32)
	      if (absz .le. 0.5*dom(ib)%dz) then 
	      	if (absy .le. 0.5*dom(ib)%dy) then
		   dom(ib)%Sp(is-1-ly,j,k)= 1.0
	      	else 
	           dom(ib)%Sp(is-1-ly,j,k)= 0.0
	      	end if
	      else 
	         dom(ib)%Sp(is-1-ly,j,k)= 0.0
	      end if	
	      end do; end do
           end if
        end if
!...............................................................................
!=== East ===>   ..  4=wall  ..    2=Outflow
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%Tbc_east.eq.4) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%Sp(ie+1+ly,j,k)= dom(ib)%Sp(ie-ly,j,k)	
              end do; end do

           else if (dom(ib)%Tbc_east.eq.2) then
              do k=ks-1,ke+1; do j=js-1,je+1
                 dom(ib)%Sp(ie+1+ly,j,k)= dom(ib)%Sp(ie-ly,j,k)
              end do; end do
           end if
        end if
!...............................................................................
!=== South ===>  ..  4=wall  ..      
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%Tbc_south.eq.4) then 
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%Sp(i,js-1-ly,k)= dom(ib)%Sp(i,js+ly,k)	
              end do; end do

           end if
        end if
!.............................................................................
!=== North ===>  ..  4=wall  ..    
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%Tbc_north.eq.4) then
              do k=ks-1,ke+1; do i=is-1,ie+1
                 dom(ib)%Sp(i,je+1+ly,k) = dom(ib)%Sp(i,je-ly,k) 	
              end do; end do

           end if
        end if
!...............................................................................
!=== Bottom ===>  ..  6=Net deposition  ..   7=Erosion
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%Tbc_bottom.eq.6) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%Sp(i,j,ks-1-ly)= 0.0	
              end do; end do

           else if (dom(ib)%Tbc_bottom.eq.7) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%Sp(i,j,ks-1-ly)= 1.5 * dom(ib)%Sp(i,j,ks+ly)	
              end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>  ..  8=free surface
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%Tbc_top.eq.8) then
              do j=js-1,je+1; do i=is-1,ie+1
                 dom(ib)%Sp(i,j,ke+1+ly) = dom(ib)%Sp(i,j,ke-ly)	
              end do; end do

           end if
        end if

!==============================================================================
        end do

	end do ! ly

        end subroutine boundSp
!#############################################################################
