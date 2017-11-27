      module connectivity                     
      implicit none
         type node_type
            integer,dimension(:),allocatable      :: node_id
            real(kind=8),dimension(:,:),allocatable :: coord
         endtype node_type

         type element_type
            integer,dimension(:),allocatable       :: element_id
            integer,dimension(:,:),allocatable     :: vertex
         endtype element_type


         character(len=*),parameter::node_file_name="node.txt"
         character(len=*),parameter::el_file_name="element.txt"
         character(len=*),parameter::BC_file_name="BC.txt"
         character(len=*),parameter::force_file_name="force.txt"
         character(len=*),parameter::crack_file_name="crack.txt"
         integer:: node_nums, element_nums, BC_nums, force_nums, &
                    crack_nums
         type(node_type)::node_data
         type(element_type)::element_data
         integer,dimension(:,:),allocatable::BC_data, force_data,&
                      crack_data
         integer,dimension(:),allocatable::new_BC_data,new_force_data
         integer,parameter :: node_file_id = 101, &
                              element_file_id   = 102, &
                              BC_file_id        = 103, &
                              force_file_id     = 104, &
                              crack_file_id     = 105, &
                              problem_dimension = 2 ,  &
                              element_dimension = 4,  &
                              BC_dimension      = 8,  &
                              force_dimension   = 8,  &              
                              crack_dimension   = 8
         real(kind=8), dimension(4,2)::node
         integer,dimension(1,4)::en 

      contains 

         subroutine element_info(el_id,element_data,node_data,en,node)
            integer,intent(in)::el_id
            type(node_type),intent(in)::node_data
            type(element_type),intent(in)::element_data
            integer,dimension(1,4),intent(out)::en
            real(kind=8),dimension(4,2),intent(out)::node

            en(1,1:4)  =element_data%vertex(el_id,1:4)
            node = node_data%coord(en(1,1:4),1:2)
         end subroutine element_info

         subroutine change_shape(BC_data,BC_nums,force_data,&
                     force_nums,new_BC_data,new_force_data)
            integer,intent(inout)::BC_nums,force_nums
            integer,dimension(BC_nums/8,8),intent(in)::BC_data
            integer,dimension(force_nums/8,8),intent(in)::force_data
            integer,dimension(:),allocatable::new_BC_data,new_force_data
            integer::i,j,bc_rows,force_rows,zero_count
            bc_rows = BC_nums/8
            force_rows = force_nums/8 
            zero_count = 0
            !do i=1,2
            !   print*, (BC_data(i,j),j=1,8)
            !enddo
            do j=0,7
               if (BC_data(bc_rows,8-j) /=0 ) then
                  zero_count = j
                  exit
               endif
            enddo
            !print*, 'number of zeros are:',zero_count
            allocate(new_BC_data(BC_nums-zero_count))
            do i=1,bc_rows-1
               new_BC_data(((i-1)*8+1):8*i) = BC_data(i,:)
            enddo
            new_BC_data((1+(bc_rows-1)*8):(8-zero_count+(bc_rows-1)*8))=&
            BC_data(bc_rows,1:(8-zero_count))
            BC_nums = BC_nums-zero_count
            
            
            zero_count = 0
            do j=0,7
               if (force_data(force_rows,8-j) /=0 ) then
                  zero_count = j
                  exit
               endif
            enddo

            allocate(new_force_data(force_nums-zero_count))
            do i=1,force_rows-1
               new_force_data(((i-1)*8+1):8*i) = force_data(i,:)
            enddo
            new_force_data((1+(bc_rows-1)*8):(8-zero_count+(bc_rows-1)*8))= &
            force_data(force_rows,1:(8-zero_count))
            force_nums = force_nums-zero_count
           ! print*, BC_nums,force_nums, zero_count
           ! print*, force_rows,BC_rows
!            do i=1,2
!               print*,(BC_data(i,j),j=1,8)
!            enddo
         end subroutine change_shape

         subroutine determine_size(node_nums,element_nums,BC_nums,&
                     force_nums, crack_nums)
            integer,intent(out) :: node_nums, element_nums, BC_nums,&
                              force_nums, crack_nums
            node_nums = get_file_size(node_file_id, &
                              node_file_name,problem_dimension)
            element_nums = get_file_size(element_file_id, &
                              el_file_name,element_dimension)
            BC_nums = get_BC_size(BC_file_id,BC_file_name)
            force_nums = get_BC_size(force_file_id,force_file_name)
            crack_nums = get_BC_size(crack_file_id,crack_file_name)
                              
            print*, 'number of the nodes are: ', node_nums         
            print*, 'number of the elements are: ', element_nums
            print*, 'number of the BC nodes are: ', BC_nums
            print*, 'number of the force nodes are: ', force_nums
            print*, 'number of the crack nodes are: ', crack_nums
         end subroutine determine_size 

         subroutine build_connectivity(node_nums, node_data, &
                        element_nums, element_data,BC_nums,  &
                        BC_data,force_nums, force_data,      &
                        crack_nums, crack_data)
            integer,intent(in) :: node_nums, element_nums, BC_nums,&
                  force_nums, crack_nums
            type(node_type),intent(out)         :: node_data
            type(element_type),intent(out)      :: element_data
            integer,dimension(:,:),allocatable, &
               intent(out):: BC_data, force_data, crack_data

            node_data=read_node_data(node_file_id,node_file_name, & 
                                 node_nums, problem_dimension)

            element_data=read_element_data(element_file_id, &
               el_file_name, element_nums, element_dimension)
            allocate(BC_data(BC_nums/8,BC_dimension))
            allocate(force_data(force_nums/8,force_dimension))
            allocate(crack_data(crack_nums/8,crack_dimension))
                       
            BC_data = read_BC_data(BC_file_id,BC_file_name,&
                              BC_nums,BC_dimension)

            force_data = read_BC_data(force_file_id,force_file_name,&
                              force_nums,force_dimension)

            crack_data = read_BC_data(crack_file_id,crack_file_name,&
                              crack_nums,crack_dimension)


          ! call print_node_data(node_nums,problem_dimension,node_data)
           ! call print_element_data(element_nums,element_dimension, &
           !                      element_data)
           ! call print_BC_data(BC_nums,8,BC_data)
           ! call print_BC_data(force_nums,8,force_data)
           ! call print_BC_data(crack_nums,crack_dimension,crack_data)

         end subroutine build_connectivity

         function get_BC_size(file_id,file_name) result(array_size)
            integer,intent(in):: file_id
            integer:: check,ierror,j
            character(len=*),intent(in)::file_name
            integer,parameter::col=8 !each line has 8 values
            integer:: line_num,array_size
            integer,dimension(1,col)::idata
           
            line_num = 0
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do 
               read(file_id,*,iostat=check) (idata(1,j),j=1,col)
               if (check > 0) then
                  print*, 'Error reading file. Please check'
                  exit
               else if (check < 0) then
                  !print*, 'End of the file'
                  exit
               else
                  line_num = line_num + 1
               endif

            enddo
            array_size = line_num*col
            close(file_id)
         end function get_BC_size

         function get_file_size(file_id,file_name,col) result(line_num)
            integer,intent(in) :: file_id,col
            character(len=*),intent(in)::file_name
            integer            :: check,line_num,ierror,i
            real,dimension(1,col)    :: idata
            line_num=0 
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do 
               read(file_id,*,iostat=check) (idata(1,i),i=1,col)
               if (check > 0) then
                  print*, 'Error reading file. Please check'
                  exit
               else if (check < 0) then
                  !print*, 'End of the file'
                  exit
               else
                  line_num = line_num + 1
               endif

            enddo
            close(file_id)
         end function get_file_size
         function read_node_data(file_id,file_name, array_size,col) &
         result(idata)
            integer,intent(in)   :: file_id, array_size, col
            character(len=*),intent(in)::file_name
            integer              :: i,j,check,ierror
            type(node_type)::idata

            allocate(idata%node_id(array_size))
            allocate(idata%coord(array_size,col))
                        
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do i=1,array_size
               read(file_id,*,iostat=check) (idata%coord(i,j), j=1,col)
               idata%node_id(i) = i
               if (check>0) then
                  print*, 'Error reading data at read_node_data'
               endif
            enddo 
            close(file_id)
         end function read_node_data
         
         function read_BC_data(file_id,file_name,array_size,idimension)&
                  result(idata)
            integer,intent(in)   :: file_id, array_size, idimension
            character(len=*),intent(in)::file_name
            integer              :: i,j,check,ierror
            integer,dimension(idimension,array_size/8):: idata_temp
            integer,dimension(array_size/8,idimension)::idata 
            integer::row
            integer,parameter::col = 8 !number of columns in txt file
            
            row = array_size/col
                        
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do i=1,row
               read(file_id,*,iostat=check) (idata_temp(j,i),j=1,col)
               !do j=1,col
               !   idata(j+8*(i-1),1) = j+8*(i-1)
               !enddo
               if (check>0) then
                  print*, 'Error reading data at read_node_data'
               endif
            enddo 
            idata = transpose(idata_temp)
           ! do i=1,row
           !    print*, (idata(i,j),j=1,8)
           ! enddo
            close(file_id)
         end function read_BC_data
        
         function read_element_data(file_id,file_name, array_size,col) &
         result(idata)
            integer,intent(in)   :: file_id, array_size, col
            character(len=*),intent(in)::file_name
            integer              :: i,j,check,ierror
            type(element_type)   :: idata

            allocate(idata%element_id(array_size))
            allocate(idata%vertex(array_size,col))
                        
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do i=1,array_size
               read(file_id,*,iostat=check) (idata%vertex(i,j), j=1,col)
               idata%element_id(i) = i
               if (check>0) then
                  print*, 'Error reading data at read_node_data'
               endif
            enddo 
            close(file_id)
         end function read_element_data
         subroutine print_element_data(array_size,col,idata)
            integer::i
            integer, intent(in)::array_size,col
            type(element_type),intent(in)::idata
            if (col==2) then
               
               do i=1,array_size
                  print*,idata%element_id(i),idata%vertex(i,1), &
                         idata%vertex(i,2)
               end do

            elseif (col==3) then
               
               do i=1,array_size
                  print*,idata%element_id(i),idata%vertex(i,1), &
                         idata%vertex(i,2),idata%vertex(i,3)
               end do
            elseif (col==4) then
               
               do i=1,array_size
                  print*,idata%element_id(i),idata%vertex(i,1), &
                         idata%vertex(i,2),idata%vertex(i,3), &
                         idata%vertex(i,4)
               end do

            endif
         end subroutine print_element_data

         subroutine print_node_data(array_size,col,idata)
            integer::i
            integer, intent(in)::array_size,col
            type(node_type),intent(in)::idata
            if (col==2) then
               
               do i=1,array_size
                  print*,idata%node_id(i),idata%coord(i,1), &
                         idata%coord(i,2)
               end do

            elseif (col==3) then
               
               do i=1,array_size
                  print*,idata%node_id(i),idata%coord(i,1), &
                         idata%coord(i,2),idata%coord(i,3)
               end do

            endif
         end subroutine print_node_data
         
         subroutine print_BC_data(array_size,col,idata)
            integer::i,j
            integer, intent(in)::array_size,col
            integer,dimension(array_size/8,col),intent(in)::idata
               
               do i=1,array_size/8
                  print*,(idata(i,j),j=1,8)
               end do

         end subroutine print_BC_data

 
      end module connectivity
