      module read_data
      implicit none

         type node_type
            integer,dimension(:),allocatable      :: node_id
            real(kind=8),dimension(:,:),allocatable :: coord
         endtype node_type

         type element_type
            integer,dimension(:),allocatable       :: element_id
            integer,dimension(:,:),allocatable     :: vertex
         endtype element_type

      contains
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
         
         function read_BC_data(file_id,file_name, array_size,col) &
         result(idata)
            integer,intent(in)   :: file_id, array_size, col
            character(len=*),intent(in)::file_name
            integer              :: i,j,check,ierror
            integer,dimension(:,:),allocatable   :: idata
            
            allocate(idata(array_size,col+1))
                        
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do i=1,array_size
               read(file_id,*,iostat=check) (idata(i,j), j=2,col)
               idata(i,1) = i
               if (check>0) then
                  print*, 'Error reading data at read_node_data'
               endif
            enddo 
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
         
         function get_BC_size(file_id,file_name) result(line_num)
            integer,intent(in):: file_id
            integer:: check,ierror
            character(len=*),intent(in)::file_name
            integer:: line_num, idata
           
            open(file_id,file=file_name,status="old",action="read", &
                  iostat=ierror)
            do 
               read(file_id,*,iostat=check) idata
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
            integer::i
            integer, intent(in)::array_size,col
            integer,dimension(array_size,col+1),intent(in)::idata
               
               do i=1,array_size
                  print*,idata(i,1),idata(i,2)
               end do

         end subroutine print_BC_data

      end module read_data

