module verlet
    use basic_fnc
    use omp_lib
    implicit none
    !rattle_položaji -> zračuna halfstep hiutrosti, unconstrained položaje, ki jih potem popravi
    !rattle_hitrosti -> ko z novimi položaji zračunamo sile, zračunamo hitrosti in jih popravimo
    !verlet1/2 -> osnovni velocity verlet algoritem

    contains 

        subroutine verlet1(N,x,a,v,delec,dt,l,ni_bar)
            implicit none
            integer,intent(in):: N
            real(8),intent(in):: dt,l,delec(4*N,2),ni_bar
            real(8),intent(inout)::x(4*N,3),a(3*N,3),v(3*N,3)
            integer::i
            !zračunamo halfstep hitrosti in uporabimo halfstep hitrosti in položaje da zračuna nove položaje
            do i=1,3*N
                v(i,:)=v(i,:)+a(i,:)*dt/(2*delec(i,1))           
                x(i,:)=x(i,:)*ni_bar +dt * v(i,:)
                x(i,:)=periodic_bound(x(i,:),l)
            end do

        end subroutine

        subroutine verlet2(N,x,a,v,delec,dt) 
                implicit none
                integer,intent(in):: N
                real(8),intent(in):: dt,delec(4*N,2)
                real(8),intent(inout)::x(4*N,3),a(3*N,3),v(3*N,3)
                integer::i
                do i=1,3*N 
                    v(i,:)=v(i,:) +dt * a(i,:) / (2.0 * delec(i,1) )
                end do
    
        end subroutine

        subroutine rattle_polozaji(N,x,a,v,delec,dt,l,dOH,dHH,virial)
            implicit none
            integer,intent(in):: N
            real(8),intent(in):: dt,l,dOH,dHH,delec(4*N,2),a(3*N,3)
            real(8),intent(inout)::x(4*N,3),v(3*N,3),virial

            real(8):: tempv(3,3),tempx(3,3),g,tempx_old(3,3),popravek(3,3),tempm(3)
            real(8):: skalarni,d_old(3),d_new(3),d_hoh(3,3)  
            real(8):: dtkv_inv,redm,odlocitev,g_vs(3)
            integer::i,j,k
        
        
        
            !zračunamo halfstep hitrosti
            do i=1,3*N
                v(i,:)=v(i,:)+a(i,:)*dt/(2*delec(i,1))
            end do
        
            dtkv_inv=dt**(-2) !inverz kvadrata dt
            !$OMP PARALLEL DO default (private) &
            !$OMP& schedule(dynamic) &
            !$OMP& SHARED(N,x,l,delec,dt,v,dhh,doh) 
            do i=1,N !uporabimo halfstep hitrosti in položaje da zračuna unconstrained položaje za eno molekulo
                odlocitev=1
                tempx_old(1,:)=x(i,:)
                tempx_old(2,:)=x(i+N,:)
                tempx_old(3,:)=x(i+2*N,:)!shranimo stare položaje, in zračunamo nove unconstrained položaje
                tempx(1,:)=x(i,:)+dt * v(i,:)
                tempx(1,:)=periodic_bound(tempx(1,:),l)
                tempx(2,:)=x(i+N,:) + dt * v(i+N,:)
                tempx(2,:)=periodic_bound(tempx(2,:),l)
                tempx(3,:)=x(i+2*N,:) + dt * v(i+2*N,:)
                tempx(3,:)=periodic_bound(tempx(3,:),l)
                
                tempm(1)=delec(i,1) !shranimo mase atomov molekule
                tempm(2)=delec(i+N,1)
                tempm(3)=delec(i+2*N,1)
                tempv(1,:)=v(i,:)
                tempv(2,:)=v(i+N,:)
                tempv(3,:)=v(i+2*N,:)
                d_hoh=doh
                d_hoh(2,3)=dhh    
                d_hoh(3,2)=dhh
            
                !zračunamo faktorje g za te tri atome

                !g_vs=0
                giteracije: do while (odlocitev > 0.00003)

                
                    do k=1,3
                        do j=k+1,3
                            redm=(1/tempm(k)+1/tempm(j))**(-1)
                            d_old=diff12(tempx_old(j,:),tempx_old(k,:),l)
                            d_new=diff12(tempx(j,:),tempx(k,:),l)
                            skalarni=dot_product(d_old,d_new)
                            g=redm * (d_hoh(k,j)**2-len_sq(d_new))/skalarni
                            
                            
                            
                            tempx_old(k,:)=tempx(k,:)
                            tempx_old(j,:)=tempx(j,:)
                            
                            
                            
                            
                            
                            tempx(k,:)=tempx(k,:)+ (1/(2*tempm(k))) * g * d_old
                            tempx(j,:)=tempx(j,:)- (1/(2*tempm(j))) * g * d_old
                        
                            tempv(k,:)=tempv(k,:)+ (1/(2*dt*tempm(k))) * g * d_old
                            tempv(j,:)=tempv(j,:)- (1/(2*dt*tempm(j))) * g * d_old
                            tempx(k,:)=periodic_bound(tempx(k,:),l)
                            tempx(j,:)=periodic_bound(tempx(j,:),l)
                        end do
                    end do
                    !zračunamo mamo vse g faktorje, lahko popravimmo pozicije
                    odlocitev=abs(len_sq(diff12(tempx(2,:),tempx(1,:),l)) -dOH**2)+ abs(len_sq(diff12(tempx(1,:),tempx(2,:),l))-dOH**2) +abs(len_sq(diff12(tempx(2,:),tempx(3,:),l))-dHH**2)  !zračunamo parameter odločitve                
                end do giteracije
                !sprejmemo nove položaje in hitrosti
                x(i,:)=tempx(1,:)
                x(i+N,:)=tempx(2,:)
                x(i+2*N,:)=tempx(3,:)
                v(i,:)=tempv(1,:)
                v(i+N,:)=tempv(2,:)
                v(i+2*N,:)=tempv(3,:)
                !preveri in upošteva periodne meje
                x(i,:)=periodic_bound(x(i,:),l)
                x(i+N,:)=periodic_bound(x(i+N,:),l)
                x(i+2*N,:)=periodic_bound(x(i+2*N,:),l)
            end do
            !$omp end parallel do
        end subroutine
        



        subroutine rattle_hit(N,x,a,v,delec,dt,l,dOH,dHH,virial) 
            implicit none
            integer,intent(in):: N
            real(8),intent(in):: dt,l,dOH,dHH,delec(4*N,2)
            real(8),intent(inout)::x(4*N,3),a(3*N,3),v(3*N,3),virial
    

            real(8):: tempv(3,3),h,tempx(3,3),popravek(3,3),tempm(3)
            real(8):: skalarni1,skalarni2,skalarni3
            real(8):: skalarni,d_old(3),v_new(3),d_hoh(3,3),h_vs(3)  
            real(8):: redm,d_kv,odlocitev
            integer::i,j,k
            
            
            
            !$OMP PARALLEL DO default (private) &
            !$OMP& schedule(dynamic) &
            !$OMP& SHARED(N,a,delec,dhh,doh,dt,l,x,v) 
            do i=1,N 
                odlocitev=1
                tempx(1,:)=x(i,:)
                tempv(1,:)=v(i,:) +dt * a(i,:) / (2.0 * delec(i,1) )
                tempx(2,:)=x(i+N,:)
                tempv(2,:)=v(i+N,:)+dt * a(i+N,:) /(2.0 * delec(i+N,1))
                tempx(3,:)=x(i+2*N,:)
                tempv(3,:)=v(i+2*N,:)+dt * a(i+2*N,:) /(2.0 * delec(i+2*N,1))
                d_hoh=doh
                d_hoh(2,3)=dhh
                tempm(1)=delec(i,1)
                tempm(2)=delec(i+N,1)
                tempm(3)=delec(i+2*N,1)
                
                !zračunamo faktorje h za te tri atome
                hiteracije: do while (odlocitev >0.001)
                    do k=1,3
                        do j=k+1,3
                            redm=(1/tempm(k)+1/tempm(j))**(-1)
                            d_old=diff12(tempx(j,:),tempx(k,:),l)
                            v_new=tempv(k,:)-tempv(j,:)
                            skalarni=dot_product(d_old,v_new)
                            d_kv=d_hoh(k,j)**2
                            h=-redm * skalarni / (d_kv)
                            tempv(k,:)=tempv(k,:)+(tempm(k))**(-1) * h * d_old
                            tempv(j,:)=tempv(j,:)-(tempm(j))**(-1) * h * d_old
                        end do
                    end do
                    !zračunamo mamo vse h faktorje, lahko popravimmo pozicije
                                    

                    skalarni1=dot_product((tempv(2,:)-tempv(1,:)),diff12(tempx(2,:),tempx(1,:),l))
                    skalarni2=dot_product((tempv(3,:)-tempv(1,:)),diff12(tempx(1,:),tempx(3,:),l))
                    skalarni3=dot_product((tempv(2,:)-tempv(3,:)),diff12(tempx(2,:),tempx(3,:),l))
                    odlocitev= abs(skalarni1) + abs(skalarni2)+abs(skalarni3) 
                                   
                end do hiteracije
            !sprejmemo nove hitrosti
            v(i,:)=tempv(1,:)
            v(i+N,:)=tempv(2,:)
            v(i+2*N,:)=tempv(3,:)
            end do
            !$omp end parallel do
        end subroutine


        
end module
