  
    program Voda_main
    use basic_fnc
    use SPME_voda_3site
    use verlet
    use termo_barostat
    implicit none
    !osnovne spremenljivke:
    real(8)::r,deltay,deltax,tempread(3),kon
    logical:: prekritje
    integer:: i,j,time_step,wr=0,samplingst=10,N !za štetje v do loop
    real(8),allocatable :: x(:,:),v(:,:),a(:,:) !x-položaji atomov, v- hitrosti atomov, a-vsota sil na delec 
    real(8),allocatable :: xx(:),delec(:,:)  !xx-placeholder za 3 random števila,delec(i,(m,q))- masa,naboj delca
    !RDF sprem.
    real(8)::del_r=0.02,rdf_max=7.5!,n_pict=0 !del_r ->delta r, rdf_max -> do kje računa RDF
    integer::n_points !število točk glede na del_r in rdf_max
    real(8),allocatable::RDF_n(:),rdf_oh(:),rdf_hh(:)
    !parametri modela:
    real(8):: masa_o, naboj_o, masa_h, naboj_h, r_hh, r_oh,dt,angle0
    real(8)::a1,a2 !a1- potencialna jama LJ, a2 razdalja LJ
    !fizikalne konstante:
    real(8)::atomska_masa=1.66053904020d-27,angstrem=1.d-10,osnovninaboj=1.602176634d-19,avog_st=6.02214076d23,k_boltz=1.380649d-23,epsilon0=8.85418782d-12,eps0,pi,kb
    real(8)::cas_ref,temp_ref,tlak_ref,hitrost_ref,sila_ref,ener_ref
    !parametri simulacije:
    integer::K,sim !parameter ewaldove vsote k_max
    real(8)::G !parameter ewaldove vsote Ewald convergence parameter
    real(8)::l ! l- dolžina škatle
    !ni_bar -> barostat scaling parameter, tau bar èasovna konstanta berendsenovega barostata, izoterm_st-konstanta pri barostatu, pressure0-nastavljen tlak
    !termo/barostati:
    real(8)::temp0,tau,th_eta,th_q !temp0-> nastavljena temperatura, tau- èasovni parameter berendsenovega termostata
    real(8)::ni_bar,tau_bar,pressure0,mtk_nib,mtk_w
    !rezultati/količine v simulaciji:
    real(8):: rho,temp,potencial_el,potencial_lj,virial,pressure,self_int !rho-gostota delcev,temp-temperatura,potencial_el-columb potencial,potencial_lj,virial,pressure,self_int-self interaction potencial
    real(8)::temp_avg=0,pressure_avg=0
 


    call random_seed()
    pi=atan(1.)*4
    
    !določitev parametrov
    
    !ewald parameter K_max
    K=5
    
    !masi v atomskih enotah
    masa_o=15.999
    masa_h=1.007
    
    !naboja, enota je osnovni naboj
    naboj_o=-0.41*2
    naboj_h=0.41
    !naboj_o=-0.42380*2.0
    !naboj_h=0.42380
    
    !dolžini vezi, enota je 1A
    r_oh=1
    angle0=109.47 !v stopinjah
    angle0=angle0*pi/180 !v radianih
    r_hh=2*r_oh*sin(angle0/2)
    
    !lj parametri, a1 je enota za energijo, tu jo napišemo v j/mol
    a1=649.65
    !a1=650.1695808
    ener_ref=1.0/avog_st
    a2=3.1656
    
    !časovni korak v sekundah
    dt=2.d-15
    !dt=5.d-16
    cas_ref=sqrt(atomska_masa*angstrem**2/ener_ref)
    dt=dt/cas_ref
    
    !temp v k
    temp_ref=ener_ref/k_boltz
    temp0=298
    temp0=temp0/temp_ref
    temp=temp0
    
    !tau v s
    tau=2.d-15
    tau=tau /cas_ref
    !tau_bar v s
    tau_bar=1.d-14
    tau_bar=tau_bar/cas_ref
     
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
        !rho=0.03327
        !rho=0.033416597
    rho=0.03337 !1g/cm^3
    rho=0.03327 !0.997 g/cm^3
    !rho=0.03317 !0.994 g/cm^3
    !rho=0.03307 !0.991 g/cm^3
    !rho=0.03297 !0.989 g/cm^3
    l=(N/rho)**(1./3.)
    !l=30
    G=5.6/l
    eps0=epsilon0 * angstrem * ener_ref /osnovninaboj**2
    !mth in nh konstante:
    tau=cas_Ref/(1e-12)
    th_q=(6*N-3)*tau**2
    tau=cas_Ref/(2e-12)
    mtk_w=(6*N-3)*tau**2
    th_eta=0
    mtk_nib=0

    allocate(x(3*N,3),v(3*N,3),a(3*N,3),xx(3),delec(3*N,2))

    n_points=floor(rdf_max/del_r)+1
    allocate(RDF_n(n_points),RDF_oh(n_points),RDF_hh(n_points))
    
    open(unit=211,status='replace',file='dens.txt')
    open(unit=201,status='replace',file='rdf.txt')
    open(unit=2011,status='replace',file='rdfOH.txt')
    open(unit=2012,status='replace',file='rdfHH.txt')
    open(unit=202,status='replace',file='results.txt') !time_step, tlak, p*V/NkT
    open(unit=203,status='replace',file='ener_pot.txt') !time_step, potencial_lj, potencial_el
    open(unit=204,status='replace',file='temp.txt') !timestep,temp
    
    
    !temp0=273+sim*1
    !temp0=temp0/temp_ref
    !temp=temp0
    !
    !RDF_n=0
    !RDF_oh=0
    !RDF_hh=0
    
    
    
    
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
    
    ! open (unit=99, file='positions128ek.txt', status='old', action='read')   !<--- če beremo položaje iz datoteke
    !!read(99,*)
    !!read(99,*)
    !do i=1,N
    !    read(99,*)tempread
    !    x(i,:)=tempread
    !    read(99,*)tempread
    !    x(i+N,:)=tempread
    !    read(99,*)tempread
    !    x(i+2*N,:)=tempread
    !end do
    !!x=x+10
    !delec(1,1)=masa_o
    !delec(1,2)=naboj_o
    !do i=2,N !nastavitev mase
    !    delec(i,1)=masa_o
    !    delec(i,2)=naboj_o
    !enddo
    !do i=1,N 
    !    delec(i+1*N,1)=masa_h
    !    delec(i+1*N,2)=naboj_h
    !    delec(i+2*N,1)=masa_h
    !    delec(i+2*N,2)=naboj_h
    !end do
    

    !izraćčun self-interaction potenciala
    self_int=0
    do i=1,3*N
        self_int=self_int+delec(i,2)**2
    end do
    kon=g/(4*eps0*pi**(3.0/2.0))
    self_int=self_int*kon

    !zacetne hitrosti, Maxwellova porazdelitev, vsak delec dobi svojo hitrost
    
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

    
    
    virial=0
    a(:,:)=0
    potencial_el=0
    potencial_lj=0
    call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int)
    call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
    !ekvilibracija
    do time_step=0,10000
        tau=dt
        mtk_nib=0 !če je ni_bar=1 se velikost škatle ne spreminja-> NVT
        virial=0
        !uporabimo rattle algoritem
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial) !posodobi hitrosti na halfstep, zračuna g-je iterativno, zračuna nove popravljene položaje
                !zračunamo sile na novo (z novimi položaji)
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
        call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int)
        
            !nove sile -> nove hitrosti
        call rattle_hit(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)  
            !zračunamo tlak
        call calc_pressure(virial,v,l,N,pressure,delec,x,a)
        call remove_comm(x,v,N,delec) !odstranimo center-of-mass motion
        call berendsen_termostat(v,temp0,temp,delec,N,tau,dt)
            !zapišemo rezultate
        if (mod(time_step,samplingst)==0) then
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5
        end if
        
          
    end do
    
    
    do sim=0,9
    temp0=237+10*sim
    temp0=temp0/temp_ref
    temp=temp0
    
    
    th_q=10
    !!NVT DEL
    do time_step=0,60000
        mtk_nib=0  
        virial=0
        !NH termostat
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        call nh_v_adjust(v,dt,th_eta,N)
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        !uporabimo rattle algoritem
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial) !posodobi hitrosti na halfstep, zračuna g-je iterativno, zračuna nove popravljene položaje
                !zračunamo sile na novo (z novimi položaji)
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
       
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
         call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int)
            !nove sile -> nove hitrosti
        call rattle_hit(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)  
            !zračunamo tlak
        call calc_pressure(virial,v,l,N,pressure,delec,x,a)
        call calc_temp(v,N,delec,temp)
        
         call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
        call  nh_v_adjust(v,dt,th_eta,N)
        call calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
         
            !zapišemo rezultate
        if (mod(time_step,samplingst)==0) then
            call calc_rdf_pict(x,N,rdf_n,rdf_oh,rdf_hh,del_r,rdf_max,n_points,l)
            wr=wr+1 !šteje število zapisov za izračun RDF
            call remove_comm(x,v,N,delec) !odstranimo center-of-mass motion
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5,sim
            write(202,*)time_step,pressure*tlak_ref/10**5,l
            write(203,*)time_step,potencial_lj*ener_ref,potencial_el*ener_ref
            write(204,*)time_step,temp*temp_ref
            write(211,*)time_step,N/(l*angstrem)**3 * 18/avog_st /1000
        end if
        
          
    end do

    
    
    
        
     th_eta=0
    mtk_nib=0   
    wr=0
    rdf_n=0
    rdf_oh=0
    rdf_hh=0
        tau=cas_Ref/(1e-12)
    th_q=(6*N-3)*tau**2*0.9
    tau=cas_Ref/(2e-12)
    mtk_w=(6*N-3)*tau**2*0.9
    
    do time_step=0,60000
        
        !MTK termo/barostat
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_v_adjust(v,th_eta,N,dt,mtk_nib)
        call calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
        call mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
        call mtk_box(l,mtk_nib,dt)
        virial=0
        !uporabimo rattle algoritem
        call shift_coms(x,dt,mtk_nib,N,l,delec)
        call rattle_polozaji(N,x,a,v,delec,dt,l,r_OH,r_HH,virial)!posodobi hitrosti na halfstep, zraèuna g-je iterativno, zraèuna nove popravljene položaje
                !zraèunamo sile na novo (z novimi položaji)
        a(:,:)=0
        potencial_el=0
        potencial_lj=0
        call lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)
        call ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int)
        
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
            !call calc_rdf_pict(x,N,rdf_n,rdf_oh,rdf_hh,del_r,rdf_max,n_points,l)
            !wr=wr+1 !šteje število zapisov za izraèun RDF
            call remove_comm(x,v,N,delec) !odstranimo center-of-mass motion
            print*,time_step,temp*temp_ref, pressure*tlak_ref/10**5,sim
            write(202,*)time_step,pressure*tlak_ref/10**5,l
            write(203,*)time_step,potencial_lj*ener_ref,potencial_el*ener_ref
            write(204,*)time_step,temp*temp_ref
            write(211,*)time_step,N/(l*angstrem)**3 * 18/avog_st /1000
        end if
        
          
    end do

    end do
     rdf_n=rdf_n/wr
    rdf_oh=rdf_oh/wr
    rdf_hh=rdf_hh/wr
    do i=1,n_points !zapišemo točke za RDFs
        write(201,*)i*del_r,rdf_n(i)
            write(2011,*)i*del_r,rdf_oh(i)
            write(2012,*)i*del_r,rdf_hh(i)
    end do
    

    
    end program Voda_main
    
    

