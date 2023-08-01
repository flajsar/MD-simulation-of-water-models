
    
    module SPME_voda_3site
    use basic_fnc
    use omp_lib
    implicit none
        !vsebuje subroutines za izračun potenciala, sil in virial
        !lj_del_pot -> zračuna lj del potenciala,sile,virial med kisiki 
        ! ewald_sum-> zračuna elektrostatkse sile kot ewaldovo vsoto
    
    
    contains
    
    subroutine lj_del_pot(x,delec,a,N,l,a1,a2,potencial_lj,virial)

        implicit none
        integer, intent(in)::N
        real(8),intent(in)::x(3*N,3),delec(3*N,2)
        real(8),intent(in)::l,a1,a2
        real(8),intent(out)::a(3*N,3),potencial_lj,virial
        !temp spremenljivke:
        real(8),allocatable::f(:,:,:)
        real(8)::x1,x2,r,r2,temppot,v1(3),v2(3),mol(3,3)
        integer::i,j
    
        allocate(f(3*N,3*N,3))
        f=0
        
        
        !zračunamo LJ potencial in sledeče sile med kisiki
        do i=1,N
            do j=i+1,N
                 r2=len_sq(diff12(x(j,:),x(i,:),l))
                 r=dsqrt(r2)
                 x1=(a2/r)**3
                 x1=x1**2
                 temppot=4*a1*x1*(x1-1)
                 x2=-24*a1*x1*(2*x1-1)
                 potencial_lj=potencial_lj+temppot
                 f(i,j,:)=-x2/r2*diff12(x(j,:),x(i,:),l)
                 f(j,i,:)=-f(i,j,:)
                 virial=virial-dot_product(f(i,j,:),diff12(x(j,:),x(i,:),l))
                 
            end do
        end do
        do i=1,3*N
            a(i,1)=a(i,1)+sum(f(i,:,1)) !sešteje sile da dobimo vsoto sil na atom i 
            a(i,2)=a(i,2)+sum(f(i,:,2))
            a(i,3)=a(i,3)+sum(f(i,:,3))
        end do
        
        end subroutine
        
        
        
        
        subroutine ewald_sum(x,delec,a,N,l,G,eps0,virial,potencial_el,k,self_int)

        implicit none
        integer, intent(in)::N,k
        real(8),intent(in)::x(3*N,3),delec(3*N,2),eps0
        real(8),intent(in)::l,G,self_int 
        real(8),intent(out)::a(3*N,3),potencial_el,virial
        !temp spremenljivke:
        real(8),allocatable::f(:,:,:)
        real(8)::x1,x2,r,r2,pi,sqpi,fakt1,fakt2,temppot,potencial_el_t
        integer::i,j,k1,k2,k3,indx
        real(8)::vec_k(3),rec_v(3,3),cntr_mass(3),d_i(3,3),force(3),mol(3,3),vs_re(3*N),vs_im(3*N)
        complex(8)::vs1,vs2,im,y1,t,vs1t
        
        allocate(f(3*N,3*N,3))
        f=0
        pi=atan(1.)*4
        sqpi=dsqrt(pi)
        potencial_el_t=0
        !!zračunamo short-range del elektrostatskega potenciala za vse atome
        !$OMP PARALLEL DO default (none) &
        !$OMP& schedule(dynamic) &
        !$OMP& SHARED(N,x,l,G,pi,eps0,f,sqpi,delec) &
        !$OMP& reduction(+:potencial_el_t) &
        !$OMP& private(r2,r,fakt1,fakt2,x1,x2,temppot,j,i)
        do i=1,3*N
            do j=i+1,3*N
                r2=len_sq(diff12(x(j,:),x(i,:),l)) !kvadrat razdalje med atomoma
                r=sqrt(r2) 
                fakt1=nint(mod(real(j-i),real(N))/real(N) +0.5 -1.0/(N+1)) !fakt1 in fakt2 sta uporabljena zato da določita če sta atoma del iste molekule
                fakt2=((fakt1+1.0)/2.0)**(-1)-1  
                x1=fakt1*erfc(G*r)-fakt2*erf(G*r) !če sta del iste molekule: fakt1=0, fakt2=1 in obratno
                x2=delec(i,2)*delec(j,2) !naboj*naboj
                temppot=x2*x1/(r*4.0*pi*eps0)

                potencial_el_t=temppot+potencial_el_t
                

                f(i,j,:)=(diff12(x(j,:),x(i,:),l))*(temppot+(4.0*pi*eps0)**(-1)*  2 * G/sqpi * x2 * exp(-G**2 * r2))/r2
                f(j,i,:)=-f(i,j,:)
                
            end do
        end do
        !$omp end parallel do
        
        do i=1,3*N
            a(i,1)=a(i,1)+sum(f(i,:,1)) !sešteje sile da dobimo vsoto sil na atom i 
            a(i,2)=a(i,2)+sum(f(i,:,2))
            a(i,3)=a(i,3)+sum(f(i,:,3))
        end do
        virial=virial-potencial_el_t
        potencial_el=potencial_el_t+potencial_el
        
        !recipročni del potenciala
        temppot=0
        vs1=(0.0,0.0)
        vs2=vs1
        im=(0.0,1.0)
        rec_v(1,:)=(/1.d0,0.d0,0.d0/)
        rec_v(2,:)=(/0.d0,1.d0,0.d0/)
        rec_v(3,:)=(/0.d0,0.d0,1.d0/)
        rec_v(1,:)=rec_v(1,:)*2*pi/l
        rec_v(2,:)=rec_v(2,:)*2*pi/l
        rec_v(3,:)=rec_v(3,:)*2*pi/l
        do k1=-k,k
            do k2=-k,k
                do k3=-k,k
                    if ((k1 .eq. 0) .and. (k2 .eq. 0) .and. (k3 .eq. 0)) then
                    else
                        indx=(k1+k * (k2 + k *k3))+1
                        vec_k=k1*rec_v(1,:)+k2*rec_v(2,:)+k3*rec_v(3,:)
                        vs1=(0.0,0.0)
                        vs2=vs1
                        fakt1=exp(-len_sq(vec_k)/(4*g**2))/len_sq(vec_k)
                        vs_re=0
                        vs_im=0
                        !$OMP PARALLEL DO &
                        !$OMP& default (shared) &
                        !$OMP& private(vs1,i)  &
                        !$OMP& reduction(+:vs_re,vs_im)
                        do i=1,3*N
                            vs1t=delec(i,2)*cdexp(im*dot_product(vec_k,x(i,:)))
                            vs_re(i)=real(vs1t,8)
                            vs_im(i)=dimag(vs1t)
                        end do
                        !$omp end parallel do
                        vs1t=sum(vs_re)+sum(vs_im)*im
                        vs2=DCONJG(vs1t)
                        t=vs1t*vs2
                        temppot=temppot+fakt1*t%re
                        do j=1,3*N
                            y1=cdexp(im*dot_product(vec_k,x(j,:))) * vs2
                            force=vec_k(:)*(delec(j,2)/(eps0*l**3) * fakt1 * y1%IM)
                            a(j,:)=a(j,:)+force
                        end do
                        end if
                end do 
            end do
        end do

        
        potencial_el=temppot/(2*eps0*l**3)+potencial_el
        virial=virial-temppot/(2*eps0*l**3)
        
        !prištejemo self interaction potencial
        potencial_el=potencial_el-self_int
        virial=virial+self_int
        
        
        !atomic to molecular virial correction
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
            virial=virial+(dot_product(a(i,:),d_i(1,:))+dot_product(a(i+n,:),d_i(2,:))+dot_product(a(i+n*2,:),d_i(3,:)))
            
        end do
        
        end subroutine
    
    
end module
    
