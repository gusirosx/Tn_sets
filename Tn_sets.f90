program TN_sets
    implicit none
    !----------------------------------------------------------
    !     Calculation of Thurgood quadrature sets for DOM
    ! G.S.Rodrigues (December 2020)
    !---------------------- Description -----------------------
    !  For more details see "The TN Quadrature Set for the 
    !     Discrete Ordinates Method" by Thurgood et al.(1995)     
    !          doi: https://doi.org/10.1115/1.2836285
    !----------------------------------------------------------
    ! <INPUT>
    ! N     - Discretization order
    ! <INTERNAL>
    ! nv    - Number of vertices
    ! nst   - Number of sub-triangles
    ! Lt    - Axis distance between two flat vertices
    ! ct    - sub-triangles counter variable
    ! ver   - Vertices of plane triangles
    ! verp  - Projected vertices in the unit sphere
    ! Cm    - Centers of gravity of plane triangles
    ! subT  - Sub-triangles vertices indicator
    ! Cm_N2 - Norm 2 of Cm
    ! ver_N2- Norm 2 of verp
    ! alpha - Solid angle verticies
    ! P     - Half perimeter
    ! T_flag- Upward/Downward triangle indicator
    ! r     - Auxiliary counter variable
    ! Tu    - Upward triangles counter variable
    ! Td    - Downward triangles counter variable
    ! Tr    - Downward triangles auxiliary counter variable
    ! aux   - Auxiliary variable
    ! <OUTPUT>
    ! mux   - x ordinate direction quadrature
    ! etay  - y ordinate direction quadrature
    ! xiz   - z ordinate direction quadrature
    ! Wq    - Sn quadrature weigth quadrature
    !----------------------------------------------------------
    integer :: i,j,k,N,nst,nv,ct,r,Tu,Td,Tr,T_flag
    double precision :: Lt,Cm_N2,ver_N2,P,aux,alpha(3)
    double precision,allocatable,dimension(:,:) :: ver,verp,Cm
    double precision, allocatable, dimension(:):: mux,etay,xiz,Wq
    integer,allocatable,dimension(:,:) :: subT
    character(len=4) :: fn
    N = 5
    nv = int((N+2)*(N+1)/2) !Number of vertices
    nst= N*N                !Number of sub-triangles
    Lt = 1.0d0/N !Axis distance between two flat vertices
    allocate(ver(3,nv),verp(3,nv),subT(3,nst),Cm(3,nst))
    allocate(mux(nst),etay(nst),xiz(nst),Wq(nst))
    !Determination of the basal triangle verticies
    k = 1
    do i = 1,N+1
        do j = 1,i
            ver(1,k) = 1 - (i-1)*Lt
            ver(2,k) = (j - 1)*Lt
            ver(3,k) = (i - j)*Lt
            k = k + 1
        end do
    end do
    !First Triangle Loop
    do ct = 1,3
        subT(ct,1) = ct
    end do
    !WARNING:For your own sake,DO NOT modify any of these values
    Tu = 1;Td = 0;Tr = 2;ct = 2;r = 0
    !Remaining Triangles Loop
    do i = 2,N
        T_flag = 0;Tu = Tu - 1
        do j = 1,i + r + 1
            if(T_flag == 0)then!Upward Triangles
                subT(1,ct) =   i + Tu 
                subT(2,ct) = 2*i + Tu 
                subT(3,ct) = subT(2,ct) + 1
                Tu = Tu + 1 
                T_flag = 1
            else if(T_flag == 1)then!Downward Triangles
                subT(1,ct) = i + Td
                subT(2,ct) = subT(1,ct) + 1
                subT(3,ct) = subT(2,ct) + Tr
                Td = Td + 1 
                T_flag = 0
            end if
            ct = ct + 1
        end do
        Tr = Tr + 1;r = r + 1
    end do
    !Determination of the centers of mass of the flat triangles
    do ct = 1,nst !Triangles Loop
        do i = 1,3!Coordinates Loop
            Cm(i,ct) = 0.0d0
            do k = 1,3 !Vertices Loop
                Cm(i,ct) = Cm(i,ct) + ver(i,subT(k,ct))
            end do
            Cm(i,ct) = Cm(i,ct)/3.0d0
        end do
    end do
    !Determination of the discrete directions cosines
    do ct = 1,nst !Triangles Loop
        Cm_N2 = 0.0d0
        do i = 1,3 !Coordinates Loop
            Cm_N2 = Cm_N2 + Cm(i,ct)**2.0d0
        end do
        Cm_N2 = sqrt(Cm_N2)
        mux(ct) = Cm(1,ct)/Cm_N2
        etay(ct)= Cm(2,ct)/Cm_N2
        xiz(ct) = Cm(3,ct)/Cm_N2
    end do
    !Projection of verticies in the unit sphere
    do k = 1,nv
        ver_N2 = 0.0d0
        do i = 1,3 !Coordinates Loop
            ver_N2 = ver_N2 + ver(i,k)**2.0d0
        end do
        ver_N2 = sqrt(ver_N2)
        do i = 1,3 !Coordinates Loop
            verp(i,k) = ver(i,k)/ver_N2
        end do
    end do
    !Determination of the discrete weights
    write (fn,'(I3.3)') N
    open(unit=5,file='T'//trim(fn)//'.dat')
    do ct = 1,nst !Triangles Loop
        alpha = 0.0d0
        do i = 1,3 !Coordinates Loop
            alpha(1) = alpha(1) + verp(i,subT(1,ct))*verp(i,subT(2,ct)) !alpha_12
            alpha(2) = alpha(2) + verp(i,subT(1,ct))*verp(i,subT(3,ct)) !alpha_13
            alpha(3) = alpha(3) + verp(i,subT(2,ct))*verp(i,subT(3,ct)) !alpha_23
        end do
        alpha = acos(alpha)
        P = 0.5d0*sum(alpha)
        aux = tan(0.5d0*P)
        do i = 1,3 !Angle Loop
            aux = aux*tan(0.5d0*(P-alpha(i)))
        end do
        Wq(ct) = 4.0d0*atan(sqrt(aux))
        write(5,'(I3,F14.8,F14.8,F14.8,F14.8)') ct,mux(ct),etay(ct),xiz(ct),Wq(ct)
    end do
    close(5)
end program TN_sets
