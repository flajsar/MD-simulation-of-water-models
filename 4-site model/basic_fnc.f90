module basic_fnc
implicit none
    contains
    !image- izracuna razliko med dvema koordinatama po mininum image konveciji
    !len_sq- izracuna kvadrat dolžine vektorja
    !diff12- izracuna razliko dveh vektorjev z upoštevanjem periodnih mej
    !periodic_bound- preveri če je atom prek meje škatle in to popravi da ga premakne na drugo stran
    !calc_pressure- z virial in hitrostmi zračuna tlak v sistemu po enačbi P=1/(3*V) * (2*K-Virial)
    !calc_rdf_pict- zračuna za posamezno sliko dodatek k  rdf 
    !i_j(n1,n2) vrne 0 če sta n1 in n2 enaka, drugače 0
    
    
    function image(x1,x2,l)
        implicit none
        real(8) :: x1,x2,l,x12,image
        x12=x2-x1
        x12=x12-l*int(2*x12/l) !če je x12 večji od l/2, se upošteva bližja slika
        image=x12
    end function image
    
    
    function len_sq(vektor)
        implicit none
        real(8):: vektor(3),len_sq
        len_sq=vektor(1)**2+vektor(2)**2+vektor(3)**2
    end function len_sq
    
    
    function diff12(v1,v2,l)
        implicit none
        real(8):: v1(3),v2(3),diff12(3),l
        diff12(1)=image(v1(1),v2(1),l)
        diff12(2)=image(v1(2),v2(2),l)
        diff12(3)=image(v1(3),v2(3),l)
    end function diff12
    
    
    function periodic_bound(v,l)
        implicit none
        real(8)::periodic_bound(3),v(3),l
        v(1)=v(1)-l*int(v(1)/l)-floor(v(1)/(1000*l))*l !del z int preveri če je večje od l in l odšteje
        v(2)=v(2)-l*int(v(2)/l)-floor(v(2)/(1000*l))*l !del z floor preveri če je negativno in l prišteje
        v(3)=v(3)-l*int(v(3)/l)-floor(v(3)/(1000*l))*l
        periodic_bound(:)=v(:)
    end function
    
    function i_j(n1,n2)
        integer::i_j,n1,n2
        
        if (n1==n2) then
            i_j=0
        else
            i_j=1
        end if 
    end function
    
    
    subroutine calc_pressure(virial,v,l,N,p,delec,x,a)
        implicit none
        integer,intent(in)::N
        real(8),intent(in)::l,virial,v(3*N,3),delec(4*N,2),x(4*N,3),a(3*N,3)
        real(8),intent(inout)::p
        
        integer::i
        real(8)::mass,v_cm(3),T,v_i(3)
        
        mass=0
        v_cm=0
        t=0
        do i=1,3*n
            v_cm=v_cm+v(i,:)*delec(i,1)
            mass=mass+delec(i,1)
            t=t+delec(i,1)*len_sq(v(i,:))
        end do
        v_cm=v_cm/mass
        t=t/(6*N-3)
        
        p=0
        do i=1,N
            v_i=(v(i,:)*delec(i,1)+v(i+n,:)*delec(i+n,1)+v(i+2*n,:)*delec(i+2*n,1))/(delec(i,1)+delec(i+n,1)+delec(i+2*n,1))
            p=p+len_sq(v_i(:)-v_cm(:))*(delec(i,1)+delec(i+n,1)+delec(i+2*n,1))
        end do
        p=p+3*t
        
        p=(p-virial)/(3*l**3)
    end subroutine
    
    subroutine calc_rdf_pict(x,N,rdf_n,rdf_oh,rdf_hh,del_r,rdf_max,n_points,l)
        implicit none
        integer,intent(in)::N,n_points
        real(8),intent(in)::x(4*N,3),del_r,rdf_max,l
        real(8),intent(inout)::rdf_n(n_points),rdf_oh(n_points),rdf_hh(n_points)
        
        integer::i,j,point
        real(8)::razd,rho,pi,c1,c2
        
        pi=atan(1.)*4
        rho=N/l**3
        do i=1,N
            do j=1,N
                razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
                if (razd<rdf_max) then 
                    point=1+floor(mod(razd,rdf_max)/del_r)
                    rdf_n(point)=rdf_n(point)+(1)*i_j(i,j)/(4*N*pi*rho*point**2*del_r**3)
                end if
                    
            end do
        end do 
        do i=2*N,3*N
            do j=2*N,3*N
                razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
                if ((razd<rdf_max) .and. ((i-j) .ne. N)) then 
                    point=1+floor(mod(razd,rdf_max)/del_r)
                    rdf_hh(point)=rdf_hh(point)+(1)*i_j(i,j)/(4*N*pi*rho*point**2*del_r**3)
                end if
                    
            end do
        end do 
        do i=1,N
            do j=2*N,3*N
                razd=sqrt(len_sq(diff12(x(i,:),x(j,:),l)))
                if ((razd<rdf_max) .and. ((i-j) .ne. N)) then  
                    point=1+floor(mod(razd,rdf_max)/del_r)
                    rdf_oh(point)=rdf_oh(point)+(1)*i_j(i,j)/(4*N*pi*rho*point**2*del_r**3)
                end if
                    
            end do
        end do 
        
    end subroutine    

    
    
    subroutine postavitev_m(x,N,r_om,l)
    implicit none
    integer,intent(in):: N
    real(8),intent(in)::r_om,l
    real(8),intent(inout)::x(4*N,3)
    real(8)::mol(3,3),vec(3),gama
    integer::i
    
    do i=1,N
        mol(1,:)=x(i,:)
        mol(2,:)=x(i+n,:)
        mol(3,:)=x(i+2*n,:)
        mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
        mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
        vec=0.5*(mol(2,:)+mol(3,:))-mol(1,:)
        gama=r_om/(dsqrt(len_sq(vec)))
        x(i+3*N,:)=mol(1,:)+gama*vec
        x(i+3*N,:)=periodic_bound(x(i+3*N,:),l)
    end do
    
    
    end subroutine
    
    
    subroutine project_realf(vf,a,N,r_om,x,l)
    implicit none
    integer,intent(in)::N
    real(8),intent(in)::vf(4*N,3),x(4*N,3),r_om,l
    real(8),intent(inout)::a(3*N,3)
    real(8)::gama,mol(4,3),f1(3),vec(3),vec1(3),fm(3)
    integer::i
    
    do i=1,3*n
        a(i,:)=a(i,:)+vf(i,:)
    end do
    
    do i=1,N
        mol(1,:)=x(i,:)
        mol(2,:)=x(i+n,:)
        mol(3,:)=x(i+2*n,:)
        mol(4,:)=x(i+3*n,:)
        mol(2,:)=mol(2,:)+l*int((mol(1,:)-mol(2,:))*2/l)
        mol(3,:)=mol(3,:)+l*int((mol(1,:)-mol(3,:))*2/l)
        mol(4,:)=mol(4,:)+l*int((mol(1,:)-mol(4,:))*2/l)
        vec=0.5*(mol(2,:)+mol(3,:))-mol(1,:)
        gama=r_om/(dsqrt(len_sq(vec)))
        vec1=mol(4,:)-mol(1,:)
        f1=vec1*(dot_product(vf(i+3*n,:),vec1)/dot_product(vec1,vec1))
        fm=vf(i+3*n,:)
        a(i,:)=a(i,:)+fm-gama*(fm-f1)
        a(i+n,:)=a(i+n,:)+gama*(fm-f1)*0.5
        a(i+2*n,:)=a(i+2*n,:)+gama*(fm-f1)*0.5
    end do
    
    
    
    end subroutine
    
end module