﻿program $4sitemodel
    use basic_fnc
    use pot_4site
    use verlet
    use termo_barostat
    implicit none
    !osnovne spremenljivke:
    real(8)::r,deltay,deltax,tempread(3),kon
    logical:: prekritje
    integer:: i,j,time_step,wr=0 !za štetje v do loop
    real(8),allocatable :: x(:,:),v(:,:),a(:,:) !x-položaji atomov, v- hitrosti atomov, a-vsota sil na delec 
    real(8),allocatable :: xx(:),delec(:,:)  !xx-placeholder za 3 random števila,delec(i,(m,q))- masa,naboj delca
    !RDF sprem.
    real(8)::del_r=0.02,rdf_max=8!,n_pict=0 !del_r ->delta r, rdf_max -> do kje raèuna RDF
    integer::n_points !število toèk glede na del_r in rdf_max
    real(8),allocatable::RDF_n(:),rdf_oh(:),rdf_hh(:)
    !parametri modela:
    real(8):: masa_o, naboj_o, masa_h, naboj_h, r_hh, r_oh,dt,angle0,r_om,masa_m,naboj_m
    real(8)::a1,a2 !a1- potencialna jama LJ, a2 razdalja LJ
    !fizikalne konstante:
    real(8)::atomska_masa=1.66053904020d-27,angstrem=1.d-10,osnovninaboj=1.602176634d-19,avog_st=6.02214076d23,k_boltz=1.380649d-23,epsilon0=8.85418782d-12,eps0,pi,kb
    real(8)::cas_ref,temp_ref,tlak_ref,hitrost_ref,sila_ref,ener_ref
    !parametri simulacije:
    integer::K,sim !parameter ewaldove vsote k_max
    real(8)::G !parameter ewaldove vsote Ewald convergence parameter
    integer:: N,samplingst=10!N-število molekul vode, samplingst-> na koliko èasovnih korakov se zapišejo vrednosti
    real(8)::l ! l- dolžina škatle
    !ni_bar -> barostat scaling parameter, tau bar èasovna konstanta berendsenovega barostata, izoterm_st-konstanta pri barostatu, pressure0-nastavljen tlak
    !termo/barostati:
    real(8)::temp0,tau,th_eta,th_q !temp0-> nastavljena temperatura, tau- èasovni parameter berendsenovega termostata
    real(8)::ni_bar,tau_bar,pressure0,mtk_nib,mtk_w
    
    !rezultati/kolièine v simulaciji:
    real(8):: rho,temp,potencial_el,potencial_lj,virial,pressure,self_int !rho-gostota delcev,temp-temperatura,potencial_el-columb potencial,potencial_lj,virial,pressure,self_int-self interaction potencial
    real(8)::temp_avg=0,pressure_avg=0
 


    call random_seed()
    pi=atan(1.)*4
    
    !dolocitev parametrov
    
    !ewald parameter K_max
    K=5
    
    !masi v atomskih enotah
    masa_o=15.999
    masa_h=1.007
    masa_m=0
    
    !naboja, enota je osnovni naboj
    naboj_o=0
    naboj_h=0.5564
    naboj_m=-0.5564*2
    
    !dolžini vezi, enota je 1A
    r_oh=0.9572
    angle0=104.52 !v stopinjah
    angle0=angle0*pi/180 !v radianih
    r_hh=2*r_oh*sin(angle0/2)
    r_om=0.1546
    
    !lj parametri, a1 je enota za energijo, tu jo napišemo v j/mol
    a1=774.26832
    ener_ref=1.0/avog_st
    a2=3.1589
    
    !èasovni korak v sekundah
    dt=2.d-15
    cas_ref=sqrt(atomska_masa*angstrem**2/ener_ref)
    dt=dt/cas_ref
    
    !temp v k
    temp_ref=ener_ref/k_boltz
    temp0=298
    temp0=temp0/temp_ref
    temp=temp0
    
    !tau v s
    tau=2.d-15
    tau=dt
     
    !tlak v bar
     tlak_ref=ener_ref/(angstrem**3.)
    pressure0=1
    pressure0=pressure0 * 10**5 / tlak_ref
    !še ostale ref konstante
    hitrost_ref=sqrt(ener_ref/atomska_masa)
    sila_ref=ener_ref/angstrem
    
    

    
    !molekule vode
    N=256
    !postavitev dolžine škatle
    rho=0.03327 !0.994 g/cm^3
    l=(N/rho)**(1./3.)
    G=5.6/l
    eps0=epsilon0 * angstrem * ener_ref /osnovninaboj**2
    !mth in nh konstante:
    tau=cas_Ref/(1e-12)
    th_q=(6*N-3)*tau**2
    tau=cas_Ref/(2e-12)
    mtk_w=(6*N-3)*tau**2
    th_eta=0
    mtk_nib=0
    

    allocate(x(4*N,3),v(3*N,3),a(3*N,3),xx(3),delec(4*N,2))

    n_points=floor(rdf_max/del_r)+1
    allocate(RDF_n(n_points),RDF_oh(n_points),RDF_hh(n_points))
    RDF_n=0
    RDF_oh=0
    RDF_hh=0

    open(unit=211,status='replace',file='dens.txt')
    open(unit=201,status='replace',file='rdf.txt')
    open(unit=2011,status='replace',file='rdfOH.txt')
    open(unit=2012,status='replace',file='rdfHH.txt')
    open(unit=202,status='replace',file='results.txt') !time_step, tlak, p*V/NkT
    open(unit=203,status='replace',file='ener_pot.txt') !time_step, potencial_lj, potencial_el
    open(unit=204,status='replace',file='temp.txt') !timestep,temp
    
    
    
    
    call random_number(xx)
    x(1,:)=(xx(:))*l
    delec(1,1)=masa_o
    delec(1,2)=naboj_o
    do i=2,N !postavitev kisikov
        prekritje=.true.
        do while (prekritje)
            call random_number(xx)
            xx=(xx)*l
            prekritje=.false.
            do j=1,i-1
                r=dsqrt(image(x(j,1),xx(1),l)**2+image(x(j,2),xx(2),l)**2+image(x(j,3),xx(3),l)**2)
                if (r<2) prekritje=.true.
            enddo
        enddo
        x(i,:)=xx(:)
        delec(i,1)=masa_o
        delec(i,2)=naboj_o
    enddo
    
    !postavitev vodikov
    do i=1,N 
        call random_number(xx)
        xx(:)=floor(xx(:)*3)
        deltax=r_hh/2
        deltay=r_oh*cos(angle0/2)
        x(i+1*N,1)=x(i,1)+deltax
        x(i+2*N,1)=x(i,1)-deltax
        x(i+1*N,2)=x(i,2)+deltay
        x(i+2*N,2)=x(i,2)+deltay
        x(i+1*N,3)=x(i,3)
        x(i+2*N,3)=x(i,3)
        x(i,:)=periodic_bound(x(i,:),l)
        x(i+N,:)=periodic_bound(x(i+N,:),l)
        x(i+2*N,:)=periodic_bound(x(i+2*N,:),l)
        
        delec(i+1*N,1)=masa_h
        delec(i+1*N,2)=naboj_h
        delec(i+2*N,1)=masa_h
        delec(i+2*N,2)=naboj_h
    end do
    
    
    ! open (unit=99, file='starting_config_xyz.txt', status='old', action='read')   !<--- če beremo položaje iz datoteke
    !!read(99,*)
    !!read(99,*)
    !do i=1,N
    !    read(99,*)tempread
    !    x(i,:)=tempread
    !    read(99,*)tempread
    !    x(i+N,:)=tempread
    !    read(99,*)tempread
    !    x(i+2*N,:)=tempread
    !    
    !end do
    !x=x+12.43864445
    !delec(1,1)=masa_o
    !delec(1,2)=naboj_o
    !do i=1,N !nastavitev mase
    !    delec(i,1)=masa_o
    !    delec(i,2)=naboj_o
    !enddo
    !do i=1,N 
    !    delec(i+1*N,1)=masa_h
    !    delec(i+1*N,2)=naboj_h
    !    delec(i+2*N,1)=masa_h
    !    delec(i+2*N,2)=naboj_h
    !    x(i,:)=periodic_bound(x(i,:),l)
    !    x(i+N,:)=periodic_bound(x(i+N,:),l)
    !    x(i+2*N,:)=periodic_bound(x(i+2*N,:),l)
    !end do
    
    
    
    
    !zaèetna postavitev M-atomov
    call postavitev_m(x,N,r_om,l)
    do i=1,N
        delec(i+3*N,1)=masa_m
        delec(i+3*N,2)=naboj_m
    end do
    
    !izracun self-interaction potenciala
    self_int=0
    do i=1,4*N
        self_int=self_int+delec(i,2)**2
    end do
    kon=g/(4*eps0*pi**(3.0/2.0))
    self_int=self_int*kon
    
    kb=k_boltz * temp_ref/(ener_ref)
    do i=1,N
        call random_number(xx)
            v(i,1) = dsqrt(temp0/(delec(i,1)+delec(i+n,1)*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
            v(i,2) = dsqrt(temp0/(delec(i,1)+delec(i+n,1)*2))*sqrt(-2*log(xx(1)))*sin(2.*pi*xx(2))
            call random_number(xx)
            v(i,3) = dsqrt(temp0/(delec(i,1)+delec(i+n,1)*2))*sqrt(-2*log(xx(1)))*cos(2.*pi*xx(2))
            
            v(i+n,:) = v(i,:)
            v(i+2*n,:) = v(i,:)
    enddo
    
    v=0
    
    virial=0
    a(:,:)=0
    potencial_el=0
    potencial_lj=0
    call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int,r_om)
    call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
    !!ekvilibracija
    do time_step=0,10000
        tau=dt
        virial=0
        !uporabimo rattle algoritem
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)
        call postavitev_m(x,N,r_om,l) !posodobi hitrosti na halfstep, zraèuna g-je iterativno, zracuna nove popravljene položaje
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
        call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int,r_om)
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
  
        
            !nove sile -> nove hitrosti
        call rattle_hit(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)
        call calc_pressure(virial,v,l,N,pressure,delec,x,a)
        
        call remove_comm(x,v,N,delec) !odstranimo center-of-mass motion
        call berendsen_termostat(v,temp0,temp,delec,N,tau,dt)
            !zapišemo rezultate
        if (mod(time_step,samplingst)==0) then
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5
        end if
    end do
    

    
    wr=0
    th_q=(6*N-3)
    !NVT DEL
    do time_step=0,60000
        virial=0
        !NH termostat
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        call nh_v_adjust(v,dt,th_eta,N)
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        !uporabimo rattle algoritem
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)
        call postavitev_m(x,N,r_om,l)!posodobi hitrosti na halfstep, zraèuna g-je iterativno, zraèuna nove popravljene položaje
                !zracunamo sile na novo (z novimi položaji)
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
        call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int,r_om)
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
            !nove sile -> nove hitrosti
        call rattle_hit(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)  
            !zraèunamo tlak
        call calc_pressure(virial,v,l,N,pressure,delec,x,a)
        call calc_temp(v,N,delec,temp)
        
         call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        call  nh_v_adjust(v,dt,th_eta,N)
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
            !zapišemo rezultate
        if (mod(time_step,samplingst)==0) then
            call calc_rdf_pict(x,N,rdf_n,rdf_oh,rdf_hh,del_r,rdf_max,n_points,l)
            wr=wr+1 !šteje število zapisov za izraèun RDF
            call remove_comm(x,v,N,delec)
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5,sim
            write(202,*)time_step,pressure*tlak_ref/10**5,l
            write(203,*)time_step,potencial_lj*ener_ref,potencial_el*ener_ref
            write(204,*)time_step,temp*temp_ref
            write(211,*)time_step,N/(l*angstrem)**3 * 18/avog_st
        end if
        
          
    end do
    
    
    tau=cas_Ref/(1e-12)
    th_q=(6*N-3)*tau**2
    mtk_nib=0
    th_eta=0
    wr=0
    tau=cas_Ref/(2e-12)
    mtk_w=(6*N-3)*tau**2*10
    !NPT
    do time_step=0,300000
        
        !MTK termo/barostat
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_v_adjust(v,th_eta,N,dt,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        call mtk_box(l,mtk_nib,dt)
        virial=0
        !uporabimo rattle algoritem
        call shift_coms(x,dt,mtk_nib,N,delec,l)
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)
        call postavitev_m(x,N,r_om,l)!posodobi hitrosti na halfstep, zraèuna g-je iterativno, zraèuna nove popravljene položaje
                !zraèunamo sile na novo (z novimi položaji)
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
        
        call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int,r_om)
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
            !nove sile -> nove hitrosti
        call rattle_hit(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)  
            !zraèunamo tlak
        
        call calc_pressure(virial,v,l,N,pressure,delec,x,a)
        call calc_temp(v,N,delec,temp)
        
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_v_adjust(v,th_eta,N,dt,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        
        
            !zapišemo rezultate
        if (mod(time_step,samplingst)==0) then
            call calc_rdf_pict(x,N,rdf_n,rdf_oh,rdf_hh,del_r,rdf_max,n_points,l)
            wr=wr+1 !šteje število zapisov za izraèun RDF
            call remove_comm(x,v,N,delec)
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5,sim
            write(202,*)time_step,pressure*tlak_ref/10**5,l
            write(203,*)time_step,potencial_lj*ener_ref,potencial_el*ener_ref
            write(204,*)time_step,temp*temp_ref
            write(211,*)time_step,N/(l*angstrem)**3 * 18/avog_st
        end if
     
        
    !end do
          
    end do
    rdf_n=rdf_n/wr
    rdf_oh=rdf_oh/wr
    rdf_hh=rdf_hh/wr
    do i=1,n_points
        write(201,*)i*del_r,rdf_n(i)
        if (i*del_r .lt. 1.5) then
             write(2011,*)i*del_r,0.0
             write(2012,*)i*del_r,0.0
        else
            write(2011,*)i*del_r,rdf_oh(i)
            write(2012,*)i*del_r,rdf_hh(i)
        end if
    end do
    
    
 end program $4sitemodel

