subroutine sor_rb(u,b,h,bcnd,res,n,n3,omega,eps,ksor,kmg,dig,rank,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve poisson equation using a weighted Gauss-
!                    Seidel method (Successive Over-Relaxation).
!     NOTE:          u(x,y,z) is defined as u(0:nx+2,0:ny+2,0:nz+2) 
!                    where 1 and n+1 are for the boundary conditions.
!                    Red-black ordering is used
!     --------------------------------------------------------------
  use mpi
  use mod_constants, only: eps0
  use mod_part_info, only: flag_pbc, flag_nmn
  implicit none
  integer:: i,is,j,k,ks,kr,kl,l,n(3),n3,shift(0:nproc-1), &
       length(0:nproc-1),ksor,kmg,dig,nr,m,u_0,u_n3p1,flag, &
       flag_c,flag_s,bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  parameter (u_0=1, u_n3p1=2)
  integer ierr,rank,nproc,status(MPI_STATUS_SIZE)
  real(kind=8):: h(3),u(0:n(1)+2,0:n(2)+2,-1:n(3)+2),mr, &
       b(0:n(1)+1,0:n(2)+1,0:n(3)+1),u_pbc(0:n(1)+2,0:n(2)+2,2)
  real(kind=8), allocatable:: u_tmp(:,:,:),b_tmp(:,:,:)
  integer, allocatable:: bcnd_tmp(:,:,:)
  ! Solver constants
  real(kind=8):: omega,eta,eps,rij,res_s,res,res_tmp ! residual
  include 'mg.h'

  !
  ! Initialization & multi-grid parameters
  !
  res_s=0.d0
  ksor=0
  nr=4
  if( kmg.eq.1 .and. dig.eq.0 ) nr=10**6
  eta=0.6d0

  !
  ! Compute size of local block
  !
  m= n3/nproc
  mr= real(n3)/real(nproc)

  flag=0
  if((real(m)-mr).eq.0.d0) flag=1

  flag_c=0 ! flag for array concatenation
  if( flag.eq.1 .and. (real(m/2)-mr/2.d0).lt.0 .and. &
       dig.lt.0 ) flag_c=1

  flag_s=0 ! flag for array splitting
  if( flag.eq.1 .and. (real(m/2)-mr/2.d0).lt.0 .and. &
       dig.gt.0 ) flag_s=1

  kl= 1
  if(flag.eq.1) then 
     if( (rank+1).eq.nproc ) then    
        kr= m+1
     else
        kr= m
     endif
  else ! Less than 1 node per proc
     m= n3
     kr= m+1
  endif

  !
  ! Split arrays u, bcnd, b when m=1 and move toward finer grids
  !
  if(flag_s.eq.1) then
     u(:,:,0:m+2)= u(:,:,rank*m:(rank+1)*m+2)
     b(:,:,0:m+1)= b(:,:,rank*m:(rank+1)*m+1) 
     bcnd(:,:,0:m+2)= bcnd(:,:,rank*m:(rank+1)*m+2)
  endif

  !
  ! SOR iteration
  !
  do l=1,nr

     !
     ! Red ordering based on old back values
     ! i=2 => j=1,3,5, etc.
     ! i=3 => j=2,4,6, etc.
     !
     
     !$OMP PARALLEL PRIVATE(is,ks) 
     !$OMP DO COLLAPSE(2)
     do k=kl,kr
        do j=1,n(2)+1
           
           ! Shift red and blacks for increasing k's
           if(mod(k+rank*m,2).eq.0) then
              ks=0
           else
              ks=1
           endif

           ! Starting points for "reds"
           if(mod(j,2).eq.0) then
              is=1
           else
              is=2
           endif
           is=is+ks
           if(is.gt.2) is=1

           do i=is,n(1)+1,2
              
              if(bcnd(i,j,k).le.0) then
                 
                 ! Neumann BC's (LHS only)
                 if( bcnd(i,j,k).eq.-2 ) then 
                    u(i,j,k)= (1.d0-omega)*u(i,j,k) + omega*(  b(i,j,k) - &
                         an*u(i,j+1,k) - as*u(i,j-1,k) - & 
                         2.d0*ae*u(i+1,j,k) - &
                         at*u(i,j,k+1) - ab*u(i,j,k-1) )/ac
                    goto 100
                 endif

                 ! Periodic BCs (calculate from 1->n)
                 if(bcnd(i,j,k).eq.0) then 
                    if(k.eq.m+1) goto 100 
                    if(j.eq.n(2)+1) goto 100 
                 endif
                 
                 ! Update u(i,j,k) [iter. k+1] inside domain     
                 u(i,j,k)= (1.d0-omega)*u(i,j,k) + omega*(  b(i,j,k) - &
                      aw*u(i-1,j,k) - an*u(i,j+1,k) - & 
                      as*u(i,j-1,k) - ae*u(i+1,j,k) - &
                      at*u(i,j,k+1) - ab*u(i,j,k-1) )/ac
                 
100              continue
              endif

           enddo
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL     

     !
     ! Send/receive ghost nodes to/from neighbour processes
     !
     if(flag.eq.1) then

        if( MOD(rank,2).eq.1 ) then ! rank is odd
           ! Send ghost value u(:,:,k=1) to rank-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, ierr)
           ! Received ghost value from rank-1 stored in u(:,:,-1:0)
           call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
           
           if(rank.lt.nproc-1) then
              ! Send ghost value u(:,:,m-1:m) to rank+1
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
              ! Received ghost value from rank+1 stored in u(:,:,k=m+1)
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
           endif
        else ! Rank is even 
           if(rank.gt.0) then
              call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, & 
                   0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
              call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                   0, MPI_COMM_WORLD, ierr)
           endif
           
           if(rank.lt.nproc-1) then
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           endif
        endif
        
     endif

     !
     ! Periodic boundary conditions
     !
     if(flag_pbc.eq.0) goto 110

     !$OMP PARALLEL
     !$OMP DO
     do k=kl-1,kr+1
        do i=0,n(1)+2
           if( bcnd(i,1,k).eq.0 .or. bcnd(i,1,k).eq.-2 ) u(i,0,k)= u(i,n(2),k)
           if( bcnd(i,n(2)+1,k).eq.0 .or. bcnd(i,n(2)+1,k).eq.-2 ) u(i,n(2)+1,k)= u(i,1,k)
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 
 
     if( nproc.gt.1 .and. flag.eq.1 ) then
        if(rank.eq.0) then
           ! Send periodic BC u(:,:,1) to nprocs-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                0, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.(nproc-1)) then
           ! Received periodic value from rank #0 : u(:,:,m+1)=u(:,:,1)
           call MPI_RECV(u_pbc(0,0,u_n3p1),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                0, MPI_COMM_WORLD, status, ierr)
           ! Send periodic BC u(:,:,m) to rank #0
           call MPI_SEND(u(0,0,m),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                1, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.0) then
           ! Received periodic value from nproc-1 : u(:,:,0)=u(:,:,m)
           call MPI_RECV(u_pbc(0,0,u_0),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                1, MPI_COMM_WORLD, status, ierr)
        endif
     endif
      
     !$OMP PARALLEL 
     !$OMP DO
     do j=0,n(2)+2
        do i=0,n(1)+2
           if( bcnd(i,j,1).eq.0 .or. bcnd(i,j,1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.0) u(i,j,0)= u_pbc(i,j,u_0)
              else
                 u(i,j,0)= u(i,j,m)
              endif
           endif
           if( bcnd(i,j,m+1).eq.0 .or. bcnd(i,j,m+1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.nproc-1) u(i,j,m+1)= u_pbc(i,j,u_n3p1)
              else
                 u(i,j,m+1)= u(i,j,1)
              endif
           endif
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 
        
110 continue

     !
     ! Black ordering based on new red values
     ! i=2 => j=2,4,6, etc.
     ! i=3 => j=1,3,5, etc.
     !     

     !$OMP PARALLEL PRIVATE(is,ks)
     !$OMP DO COLLAPSE(2)
     do k=kl,kr
        do j=1,n(2)+1
           
           ! Shift red and blacks for increasing k's
           if(mod(k+rank*m,2).eq.0) then
              ks=0
           else
              ks=1
           endif

           ! Starting points for "blacks"
           if(mod(j,2).eq.0) then
              is=2
           else
              is=1
           endif
           is=is+ks
           if(is.gt.2) is=1
           
           do i=is,n(1)+1,2
              
              if(bcnd(i,j,k).le.0) then
                 
                 ! Neumann BC's (LHS only)
                 if( bcnd(i,j,k).eq.-2 ) then 
                    u(i,j,k)= (1.d0-omega)*u(i,j,k) + omega*(  b(i,j,k) - &
                         an*u(i,j+1,k) - as*u(i,j-1,k) - & 
                         2.d0*ae*u(i+1,j,k) - &
                         at*u(i,j,k+1) - ab*u(i,j,k-1) )/ac
                    goto 120
                 endif

                 ! Periodic BCs (calculate from 1->n)
                 if(bcnd(i,j,k).eq.0) then 
                    if(k.eq.m+1) goto 120 
                    if(j.eq.n(2)+1) goto 120 
                 endif
                 
                 ! Update u(i,j,k) [iter. k+1] inside domain     
                 u(i,j,k)= (1.d0-omega)*u(i,j,k) + omega*(  b(i,j,k) - &
                      aw*u(i-1,j,k) - an*u(i,j+1,k) - & 
                      as*u(i,j-1,k) - ae*u(i+1,j,k) - &
                      at*u(i,j,k+1) - ab*u(i,j,k-1) )/ac
                 
120              continue
              endif
              
           enddo
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL

     !
     ! Send/receive ghost nodes to/from neighbour processes
     !
     if(flag.eq.1) then

        if(MOD(rank,2).eq.1) then ! rank is odd
           ! Send ghost value u(:,:,k=1) to rank-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, ierr)
        ! Received ghost value from rank-1 stored in u(:,:,-1:0)
           call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
           
           if(rank.lt.nproc-1) then
              ! Send ghost value u(:,:,m-1:m) to rank+1
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           ! Received ghost value from rank+1 stored in u(:,:,k=m+1)
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
           endif
        else ! Rank is even 
           if(rank.gt.0) then
              call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, & 
                   0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
              call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                   0, MPI_COMM_WORLD, ierr)
           endif
           
           if(rank.lt.nproc-1) then
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           endif
        endif

     endif

     !
     ! Periodic boundary conditions
     !     
     if(flag_pbc.eq.0) goto 130

     !$OMP PARALLEL
     !$OMP DO
     do k=kl-1,kr+1
        do i=0,n(1)+2
           if( bcnd(i,1,k).eq.0 .or. bcnd(i,1,k).eq.-2 ) u(i,0,k)= u(i,n(2),k)
           if( bcnd(i,n(2)+1,k).eq.0 .or. bcnd(i,n(2)+1,k).eq.-2 ) u(i,n(2)+1,k)= u(i,1,k)
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 

     if( nproc.gt.1 .and. flag.eq.1 ) then
        if(rank.eq.0) then
           ! Send periodic BC u(:,:,1) to nprocs-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                0, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.(nproc-1)) then
           ! Received periodic value from rank #0 : u(:,:,m+1)=u(:,:,1)
           call MPI_RECV(u_pbc(0,0,u_n3p1),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                0, MPI_COMM_WORLD, status, ierr)
           ! Send periodic BC u(:,:,m) to rank #0
           call MPI_SEND(u(0,0,m),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                1, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.0) then
           ! Received periodic value from nproc-1 : u(:,:,0)=u(:,:,m)
           call MPI_RECV(u_pbc(0,0,u_0),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                1, MPI_COMM_WORLD, status, ierr)
        endif
     endif
       
     !$OMP PARALLEL 
     !$OMP DO
     do j=0,n(2)+2
        do i=0,n(1)+2
           if( bcnd(i,j,1).eq.0 .or. bcnd(i,j,1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.0) u(i,j,0)= u_pbc(i,j,u_0)
              else
                 u(i,j,0)= u(i,j,m)
              endif
           endif
           if( bcnd(i,j,m+1).eq.0 .or. bcnd(i,j,m+1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.nproc-1) u(i,j,m+1)= u_pbc(i,j,u_n3p1)
              else
                 u(i,j,m+1)= u(i,j,1)
              endif
           endif
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 

130  continue

     !
     ! Calculate residual
     !
    
     res=0.d0
     !$OMP PARALLEL PRIVATE(rij) REDUCTION(+:res)
     !$OMP DO COLLAPSE(2)
     do k=kl,kr
        do j=1,n(2)+1
           do i=1,n(1)+1
              
              rij=0.d0              
              ! Interior points
              if(bcnd(i,j,k).eq.-1) then
                 rij= b(i,j,k) - ( aw*u(i-1,j,k) + an*u(i,j+1,k) +  &
                      as*u(i,j-1,k) + ae*u(i+1,j,k) + ac*u(i,j,k) + &
                      at*u(i,j,k+1) + ab*u(i,j,k-1) )
              endif
              
              res= res + ABS(rij)
              
           enddo
        enddo
     enddo

     !$OMP END DO NOWAIT
     !$OMP END PARALLEL

     !
     !  Sum residual between all processes
     !
     if(flag.eq.1) then
        call MPI_ALLREDUCE(res, res_tmp, 1, MPI_REAL8, MPI_SUM, &
             MPI_COMM_WORLD, ierr)
        res= res_tmp
     endif


     !
     ! Convergence test
     !
     ksor= ksor + 1 

     if(dig.eq.0) then
        if( res.le.eps ) exit
     else
        if (abs(res_s) .gt. max(tiny(1.d0),1.d-300)) then
           if ((res/res_s) .le. eta) exit
        endif
     endif
     
     res_s= res

  enddo

  !
  ! Concatenate arrays u, bcnd, b when m=1 and move toward coarser grids
  !
  if(flag_c.eq.1) then  

     allocate ( u_tmp(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
          b_tmp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
          bcnd_tmp(0:n(1)+2,0:n(2)+2,0:n(3)+2) )

     ! u()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+2
     enddo
     length(0)=length(0)+2
     length(nproc-1)=length(nproc-1)+2
     
     length= length*(n(1)+3)*(n(2)+3)
     shift= shift*(n(1)+3)*(n(2)+3)

     kl= 1
     if(rank.eq.0) kl= kl-2
     kr= m
     if(rank.eq.nproc-1) kr= kr+1
     
     call MPI_ALLGATHERV(u(0:n(1)+2,0:n(2)+2,kl:kr),length(rank),MPI_REAL8,u_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     u= u_tmp

     ! b()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+1
     enddo
     length(0)=length(0)+1
     length(nproc-1)=length(nproc-1)+1
     
     length= length*(n(1)+2)*(n(2)+2)
     shift= shift*(n(1)+2)*(n(2)+2)

     kl= 1
     if(rank.eq.0) kl= kl-1
     kr= m
     if(rank.eq.nproc-1) kr=kr+1

     call MPI_ALLGATHERV(b(0:n(1)+1,0:n(2)+1,kl:kr),length(rank),MPI_REAL8,b_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     b= b_tmp

     ! bcnd()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+1
     enddo
     length(0)=length(0)+1
     length(nproc-1)=length(nproc-1)+2
     
     length= length*(n(1)+3)*(n(2)+3)
     shift= shift*(n(1)+3)*(n(2)+3)

     kl= 1
     if(rank.eq.0) kl=kl-1
     kr= m
     if(rank.eq.nproc-1) kr=kr+1
     
     call MPI_ALLGATHERV(bcnd(0:n(1)+2,0:n(2)+2,kl:kr),length(rank),MPI_INTEGER,bcnd_tmp, &
          length,shift,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     bcnd= bcnd_tmp

     deallocate ( u_tmp, b_tmp, bcnd_tmp )

  endif

  return
end subroutine sor_rb


subroutine jacobi(u,b,h,bcnd,res,n,n3,eps,ksor,kmg,dig,rank,nproc)
!     ==============================================================
!     VERSION:         0.1
!     LAST MOD:      DEC/14
!     MOD AUTHOR:    G. Fubiani
!     COMMENTS:      Solve poisson equation using a Jacobi's method.
!     NOTE:          u(x,y,z) is defined as u(0:nx+2,0:ny+2,0:nz+2) 
!                    where 1 and n+1 are for the boundary conditions.
!     --------------------------------------------------------------
  use mpi
  use mod_constants, only: eps0
  use mod_part_info, only: flag_pbc, flag_nmn
  implicit none
  integer:: i,j,k,l,kl,kr,n(3),n3,shift(0:nproc-1),length(0:nproc-1), &
       ksor,kmg,dig,nr,m,u_0,u_n3p1,flag,flag_c,flag_s, &
       bcnd(0:n(1)+2,0:n(2)+2,0:n(3)+2)
  parameter (u_0=1, u_n3p1=2)
  integer ierr,rank,nproc,status(MPI_STATUS_SIZE)
  real(kind=8):: h(3),mr,u_pbc(0:n(1)+2,0:n(2)+2,2), &
       uold(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
       u(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
       b(0:n(1)+1,0:n(2)+1,0:n(3)+1)
  real(kind=8), allocatable:: u_tmp(:,:,:),b_tmp(:,:,:)
  integer, allocatable:: bcnd_tmp(:,:,:)
  ! Solver constants
  real(kind=8):: eta,eps,rij,res_s,res,res_tmp ! residual
  include 'mg.h'

  !
  ! Initialization & multi-grid parameters
  !
  res_s=0.d0
  ksor=0
  nr=4
  if( kmg.eq.1 .and. dig.eq.0 ) nr=10**6
  eta=0.6d0

  !
  ! Compute size of local block
  !
  m= n3/nproc
  mr= real(n3)/real(nproc)

  flag=0
  if((real(m)-mr).eq.0.d0) flag=1

  flag_c=0 ! flag for array concatenation
  if( flag.eq.1 .and. (real(m/2)-mr/2.d0).lt.0 .and. &
       dig.lt.0 ) flag_c=1

  flag_s=0 ! flag for array splitting
  if( flag.eq.1 .and. (real(m/2)-mr/2.d0).lt.0 .and. &
       dig.gt.0 ) flag_s=1
  
  kl= 1
  if(flag.eq.1) then 
     if( (rank+1).eq.nproc ) then    
        kr= m+1
     else
        kr= m
     endif
  else ! Less than 1 node per proc
     m= n3
     kr= m+1
  endif

  !
  ! Split arrays u, bcnd, b when m=1 and move toward finer grids
  !
  if(flag_s.eq.1) then 
     u(:,:,0:m+2)= u(:,:,rank*m:(rank+1)*m+2)
     b(:,:,0:m+1)= b(:,:,rank*m:(rank+1)*m+1) 
     bcnd(:,:,0:m+2)= bcnd(:,:,rank*m:(rank+1)*m+2)
  endif

  !
  ! SOR iteration
  !
  do l=1,nr

     !$OMP PARALLEL 
     !$OMP DO COLLAPSE(2)
     do k=kl-1,kr+1
        do j=0,n(2)+2           
           do i=0,n(1)+2                                  
                 uold(i,j,k)= u(i,j,k)
           enddo
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 

     res=0.d0
     !$OMP PARALLEL PRIVATE(rij) REDUCTION(+:res) 
     !$OMP DO COLLAPSE(2) 
     do k=kl,kr
        do j=1,n(2)+1           
           do i=1,n(1)+1
              
              if(bcnd(i,j,k).le.0) then
                 
                 ! Neumann BC's (LHS only)
                 if( bcnd(i,j,k).eq.-2 ) then 
                    rij= b(i,j,k) - &
                         an*uold(i,j+1,k) -  as*uold(i,j-1,k) - &
                         2.d0*ae*uold(i+1,j,k) -  &
                         at*uold(i,j,k+1) - ab*uold(i,j,k-1) - &
                         ac*uold(i,j,k)
                    goto 100
                 endif

                 ! Periodic BCs (calculate from 1->n)
                 if(bcnd(i,j,k).eq.0) then 
                    if(k.eq.m+1) goto 100 
                    if(j.eq.n(2)+1) goto 100 
                 endif
        
                 rij= b(i,j,k) - &
                      aw*uold(i-1,j,k) - an*uold(i,j+1,k) -  &
                      as*uold(i,j-1,k) - ae*uold(i+1,j,k) -  &
                      ac*uold(i,j,k) - at*uold(i,j,k+1) - &
                      ab*uold(i,j,k-1)
                 
                 res= res + ABS(rij)
                 
                 ! Update u(i,j,k) [iter. k+1] inside domain     
                 u(i,j,k)= uold(i,j,k) + rij/ac
                
100              continue

              endif

           enddo
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 

     !
     ! Send/receive ghost nodes to/from neighbour processes
     !
     if(flag.eq.1) then

        if(MOD(rank,2).eq.1) then ! rank is odd
           ! Send ghost value u(:,:,k=1) to rank-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, ierr)
           ! Received ghost value from rank-1 stored in u(:,:,-1:0)
           call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, &
                0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
           
           if(rank.lt.nproc-1) then
              ! Send ghost value u(:,:,m-1:m) to rank+1
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
              ! Received ghost value from rank+1 stored in u(:,:,k=m+1)
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
           endif
        else ! Rank is even 
           if(rank.gt.0) then
              call MPI_RECV(u(0,0,-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank-1, & 
                   0, MPI_COMM_WORLD, status, ierr) ! 2 LHS ghost nodes for MG
              call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank-1, &
                   0, MPI_COMM_WORLD, ierr)
           endif
           
           if(rank.lt.nproc-1) then
              call MPI_RECV(u(0,0,m+1),(n(1)+3)*(n(2)+3), MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, status, ierr)
              call MPI_SEND(u(0,0,m-1),(n(1)+3)*(n(2)+3)*2, MPI_REAL8, rank+1, &
                   0, MPI_COMM_WORLD, ierr) ! 2 LHS ghost nodes for MG
           endif
        endif

     endif

     !
     ! Periodic boundary conditions
     !     
     if(flag_pbc.eq.0) goto 110 
        
     !$OMP PARALLEL
     !$OMP DO
     do k=kl-1,kr+1
        do i=0,n(1)+2
           if( bcnd(i,1,k).eq.0 .or. bcnd(i,1,k).eq.-2 ) u(i,0,k)= u(i,n(2),k)
           if( bcnd(i,n(2)+1,k).eq.0 .or. bcnd(i,n(2)+1,k).eq.-2 ) u(i,n(2)+1,k)= u(i,1,k)
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 

     if( nproc.gt.1 .and. flag.eq.1 ) then
        if(rank.eq.0) then
           ! Send periodic BC u(:,:,1) to nprocs-1
           call MPI_SEND(u(0,0,1),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                0, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.(nproc-1)) then
           ! Received periodic value from rank #0 : u(:,:,m+1)=u(:,:,1)
           call MPI_RECV(u_pbc(0,0,u_n3p1),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                0, MPI_COMM_WORLD, status, ierr)
           ! Send periodic BC u(:,:,m) to rank #0
           call MPI_SEND(u(0,0,m),(n(1)+3)*(n(2)+3), MPI_REAL8, 0, &
                1, MPI_COMM_WORLD, ierr)
        endif
        if(rank.eq.0) then
           ! Received periodic value from nproc-1 : u(:,:,0)=u(:,:,m)
           call MPI_RECV(u_pbc(0,0,u_0),(n(1)+3)*(n(2)+3), MPI_REAL8, nproc-1, &
                1, MPI_COMM_WORLD, status, ierr)
        endif
     endif
     
     !$OMP PARALLEL 
     !$OMP DO
     do j=0,n(2)+2
        do i=0,n(1)+2
           if( bcnd(i,j,1).eq.0 .or. bcnd(i,j,1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.0) u(i,j,0)= u_pbc(i,j,u_0)
              else
                 u(i,j,0)= u(i,j,m)
              endif
           endif
           if( bcnd(i,j,m+1).eq.0 .or. bcnd(i,j,m+1).eq.-2 ) then
              if( nproc.gt.1 .and. flag.eq.1 ) then
                 if(rank.eq.nproc-1) u(i,j,m+1)= u_pbc(i,j,u_n3p1)
              else
                 u(i,j,m+1)= u(i,j,1)
              endif
           endif
        enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP END PARALLEL 
             
110  continue     

     !
     !  Sum residual between all processes
     !
     if(flag.eq.1) then
        call MPI_ALLREDUCE(res, res_tmp, 1, MPI_REAL8, MPI_SUM, &
             MPI_COMM_WORLD, ierr)
        res= res_tmp
     endif

     !
     ! Convergence test
     !
     ksor= ksor + 1 

     if(dig.eq.0) then
        if( res.le.eps ) exit
     else
      if (abs(res_s) .gt. max(tiny(1.d0),1.d-300)) then
         if ((res/res_s) .le. eta) exit
      endif
     endif
     
     res_s= res

  enddo

  !
  ! Concatenate arrays u, bcnd, b when m=1 and move toward coarser grids
  !
  if(flag_c.eq.1) then  

     allocate ( u_tmp(0:n(1)+2,0:n(2)+2,-1:n(3)+2), &
          b_tmp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
          bcnd_tmp(0:n(1)+2,0:n(2)+2,0:n(3)+2) )

     ! u()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+2
     enddo
     length(0)=length(0)+2
     length(nproc-1)=length(nproc-1)+2
     
     length= length*(n(1)+3)*(n(2)+3)
     shift= shift*(n(1)+3)*(n(2)+3)

     kl= 1
     if(rank.eq.0) kl= kl-2
     kr= m
     if(rank.eq.nproc-1) kr= kr+1
     
     call MPI_ALLGATHERV(u(0:n(1)+2,0:n(2)+2,kl:kr),length(rank),MPI_REAL8,u_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     u= u_tmp

     ! b()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+1
     enddo
     length(0)=length(0)+1
     length(nproc-1)=length(nproc-1)+1
     
     length= length*(n(1)+2)*(n(2)+2)
     shift= shift*(n(1)+2)*(n(2)+2)

     kl= 1
     if(rank.eq.0) kl= kl-1
     kr= m
     if(rank.eq.nproc-1) kr=kr+1

     call MPI_ALLGATHERV(b(0:n(1)+1,0:n(2)+1,kl:kr),length(rank),MPI_REAL8,b_tmp, &
          length,shift,MPI_REAL8,MPI_COMM_WORLD,ierr)
     b= b_tmp

     ! bcnd()
     shift(0)=0
     do i=0,nproc-1
        length(i)=m
        if(i.ge.1) shift(i)=i*m+1
     enddo
     length(0)=length(0)+1
     length(nproc-1)=length(nproc-1)+2
     
     length= length*(n(1)+3)*(n(2)+3)
     shift= shift*(n(1)+3)*(n(2)+3)

     kl= 1
     if(rank.eq.0) kl=kl-1
     kr= m
     if(rank.eq.nproc-1) kr=kr+1
     
     call MPI_ALLGATHERV(bcnd(0:n(1)+2,0:n(2)+2,kl:kr),length(rank),MPI_INTEGER,bcnd_tmp, &
          length,shift,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     bcnd= bcnd_tmp

     deallocate ( u_tmp, b_tmp, bcnd_tmp )

  endif

  return
end subroutine jacobi

