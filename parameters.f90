      module parameters
      implicit none
         !-- Input parameters

         real(kind=8),parameter::small_k=10**(-6), &
               G_c = 0.000027,   & !energy release rate
               l   = 0.02,       & !length scale
               E   = 1, & !elastic modulus
               miu   = 0.3,        & !poisson's ratio
               rho = 1,     & !mass density
               lamda = E*miu/(1+miu)/(1-2*miu), & !Lame constant
               velocity= sqrt(E*(1-miu)/(1+miu)/(1-2*miu)/rho)
         real(kind=8),dimension(3,3)::mat_mtx !material matrix
      contains
         function material_matrix(E,miu) result(mat_mtx)
            real(kind=8),intent(in)::E,miu
            real(kind=8),dimension(3,3)::mat_mtx
            integer::i,j
            mat_mtx(1,1) = 1-miu
            mat_mtx(1,2) = miu
            mat_mtx(1,3) = 0.0
            mat_mtx(2,1) = miu
            mat_mtx(2,2) = 1-miu
            mat_mtx(2,3) = 0
            mat_mtx(3,1) = 0
            mat_mtx(3,2) = 0
            mat_mtx(3,3) = 1-2*miu
            mat_mtx = E/(1+miu)/(1-2*miu)*mat_mtx
            print*, 'Material matrix is:'
            do i=1,3
               write(*,'(3F10.4)') (mat_mtx(i,j),j=1,3)
            enddo
         end function material_matrix
         !-- 
                       
      end module 
            

