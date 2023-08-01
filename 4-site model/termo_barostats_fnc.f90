module termo_barostat
    use basic_fnc
    implicit none
    
    
    contains
    
    subroutine berendsen_termostat(v,temp0,temp,delec,N,tau,dt)
        implicit none
        integer,intent(in)::N
        real(8),intent(in)::temp0,delec(4*N,2),tau,dt
        real(8),intent(inout)::v(3*N,3),temp
        real(8)::ksi
        integer::i
        
        temp=0
        do i=1,3*N
            temp=temp+delec(i,1)*len_sq(v(i,:)) !zračuna temp kot vsoto m*v**2
        end do
        temp=temp/(6*N-3) !prejšna vsota se deli z 3*(N-1)-N_c kjer je N_c število constraints-v tem primeru so 3 na molekulo
        ksi=sqrt(1-dt/tau *(1-temp0/temp)) !berendsenov faktor
        v(:,:)=ksi*v(:,:)
    end subroutine
    
    
    
    
    
    
    subroutine berendsen_barostat(pressure,pressure0,ni_bar,tau_bar,l,izoterm_st,dt)
        implicit none
        real(8),intent(in)::pressure,pressure0,tau_bar,izoterm_st,dt
        real(8),intent(inout)::ni_bar,l
        ni_bar=(1-izoterm_st * dt/tau_bar *(pressure0 - pressure))**(1./3.)
        l=l*ni_bar
    end subroutine
    
     subroutine remove_comm(x,v,N,delec)
        implicit none
        real(8),intent(in)::x(4*n,3),delec(4*N,2)
        integer,intent(in)::N
        real(8),intent(inout)::v(3*n,3)
        
        integer::i
        real(8)::comm(3),mass
        comm=0
        mass=0
        do i=1,3*N
            comm=comm+v(i,:)*delec(i,1)
            mass=mass+delec(i,1)
        end do
        comm=comm/mass
        do i=1,3*n
            v(i,:)=v(i,:)-comm
        end do
        
        
     end subroutine
     
     
     
     
     subroutine calc_th_eta(th_q,th_eta,v,dt,delec,N,temp0)
     integer,intent(in)::N
     real(8),intent(in)::dt,th_q,delec(4*N,2),temp0,v(3*N,3)
     real(8),intent(inout)::th_eta
     integer::i
     real(8)::sum=0
     sum=0
     do i=1,3*N
        sum=sum+delec(i,1)*len_sq(v(i,:))
     end do
     th_eta=th_eta + dt* (sum-(6*N-3)*temp0) /(4*th_q)
     end subroutine
     
     subroutine nh_v_adjust(v,dt,th_eta,N)
     real(8),intent(in)::dt,th_eta
     integer, intent(in):: N
     real(8),intent(inout):: v(3*N,3)
     
        v=v-v*dt*th_eta/2
     end subroutine
     
     
     subroutine calc_temp(v,N,delec,temp)
     integer,intent(in)::N
     real(8),intent(in)::delec(4*N,2),v(3*N,3)
     real(8),intent(out)::temp
     integer::i
     temp=0
     do i=1,3*N
         temp=temp+delec(i,1)*len_sq(v(i,:)) 
    end do
    temp=temp/(6*N-3)
     end subroutine
     
     
     subroutine mtk_calc_th_eta(th_eta,dt,th_q,mtk_w,N,v,delec,temp0,mtk_nib)
     implicit none
     integer,intent(in):: N
     real(8),intent(in)::dt,th_q,mtk_w,v(3*N,3),temp0,delec(4*n,2),mtk_nib
     real(8),intent(inout)::th_eta
     real(8)::vs,nf
     integer::i
     nf=6*n-3
     vs=0
     do i=1,3*n
         vs=vs+delec(i,1)*len_sq(v(i,:))
     end do
     
     th_eta=th_eta+ dt/(4*th_q) * (  vs - nf*temp0 + mtk_w * mtk_nib**2-temp0)
     
     end subroutine
     
     
     subroutine calc_mtk_nib(th_eta,dt,mtk_w,N,v,delec,virial,l,mtk_nib,pressure0,x,a)
     implicit none
     real(8),intent(inout)::mtk_nib
     integer,intent(in):: N
     real(8),intent(in):: virial,v(3*n,3),l,delec(4*n,2),x(4*n,2),a(3*n,3),dt,mtk_w,pressure0,th_eta
     real(8)::p,t   
     
     call calc_pressure(virial,v,l,N,p,delec,x,a)
    call calc_temp(v,N,delec,t)
    
     mtk_nib=mtk_nib + dt/4 * ( 3/mtk_w * l**3 * (p-pressure0) + 3/mtk_w * t - th_eta* mtk_nib)
     end subroutine
     
     
     subroutine mtk_v_adjust(v,th_eta,N,dt,mtk_nib)
     implicit none
     integer,intent(in)::n
     real(8),intent(in)::th_eta,dt,mtk_nib
     real(8),intent(inout)::v(3*N,3)
     integer::i
     do i=1,3*N
         v(i,:)=v(i,:) - dt*(th_eta +(1+3/(6*N-3))* mtk_nib)*v(i,:)/2
     end do
     
     
     end subroutine
     
     
     subroutine mtk_box(l,mtk_nib,dt)
     implicit none
     real(8),intent(in)::mtk_nib,dt
     real(8),intent(inout)::l
     
     l=l * dexp(dt*mtk_nib)
     
     
     
     end subroutine
     
     
     
     subroutine shift_coms(x,dt,mtk_nib,N,delec,l)
     real(8),intent(in)::dt,mtk_nib,delec(4*n,2),l
     integer,intent(in)::N
     real(8),intent(inout)::x(4*n,3)
     integer::i
     real(8)::mol(3,3),cntr_mass(3),d_i(3,3)
     
     do i=1,N
         mol(1,:)=x(i,:)
         mol(2,:)=x(i+n,:)
         mol(3,:)=x(i+2*n,:)
         mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
         mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
         cntr_mass=delec(i,1)*mol(1,:)+delec(i+n,1)*mol(2,:)+delec(i+n*2,1)*mol(3,:)
         cntr_mass=cntr_mass/(delec(i,1)+delec(i+n,1)+delec(i+n*2,1))
         d_i(1,:)=mol(1,:)-cntr_mass
         d_i(2,:)=mol(2,:)-cntr_mass
         d_i(3,:)=mol(3,:)-cntr_mass
         cntr_mass=cntr_mass*(1+dt*mtk_nib)
         x(i,:)=cntr_mass+ d_i(1,:)
         x(i+n,:)=cntr_mass+ d_i(2,:)
         x(i+2*n,:)=cntr_mass+ d_i(3,:)
         x(i,:)=periodic_bound(x(i,:),l)
            x(i+N,:)=periodic_bound(x(i+N,:),l)
            x(i+2*N,:)=periodic_bound(x(i+2*N,:),l)
     end do
     end subroutine
     
    end module
        
        
        
        