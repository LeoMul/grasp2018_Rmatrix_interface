program read 
    !read_grasp2k.f90 
    !gfortran read_grasp2k.f90 -o read_grasp2k.x
    !author: Leo Patrick Mulholland, QUB PhD Student, ARC

    !just puts the grasp2k orbitals in the format 
    !radial grid(i), (large(i), small(i)) for each orbital 

    !e.g col 1 = radial grid
    !e.g col 2,3 = large,small of first orbital
    !e.g col 4,5 = large,small of second orbital
    !etc

    implicit none

    integer :: i,j,II
    character*6 :: header
    
    integer,parameter :: NUM_ORBITALS_FIXED   = 50
    integer,parameter :: NUM_RAD_POINTS_FIXED = 500

    integer :: NP(NUM_ORBITALS_FIXED),NAK(NUM_ORBITALS_FIXED),MFJ(NUM_ORBITALS_FIXED)
    real*8 :: E(NUM_ORBITALS_FIXED)
    real*8,allocatable :: large_component(:,:),small_component(:,:),radial_grid(:)
    real*8 ::PZ

    character*2 :: NH_ARRAY(NUM_ORBITALS_FIXED)
    character*2 :: angular_symbols_neg(5)
    character*2 :: angular_symbols_pos(4)
    integer ::num_orbitals, num_points 


    angular_symbols_neg = (/'S ','P ','D ','F ','G ' /)
    angular_symbols_pos = (/'P-','D-','F-','G-' /)

    write(*,100)
    print*,'Entering routine read_orbitals'
    call read_orbitals()    
    
    write(*,100)
    print*,'Entering routine write_out'
    call write_out()
    write(*,100)
    print*,'STOPPING NORMALLY'

    100 format('------------------------------------')

    contains

    subroutine write_out()

        !writes TARGET.INP
        !this file contains:
        !generalised occupation numbers
        !the largest radial grid used to store the orbitals
        !the relativistic orbitals, 
                !n, kappa
                !large component 
        !then 
                !n, kappa
                !small component 
        
        !where the components have been extended by zeros to fit into the complete grid described above.


        open(2,file='orbitalsgrasp2k.dat')
        print*,'Writing to orbitalsgrasp2k.dat'
        !write(2,3040) num_orbitals,num_points
        !WRITE (2,3050) (radial_grid(I),I=1,num_points)

        write(2,3030) '  #Radial Grid',(i*2,NP(i),nh_array(i),'large',i*2+1,NP(i),nh_array(i),'small',i=1,num_orbitals)

        do i = 1,num_points
            WRITE (2,3050) radial_grid(i), (large_component(j,i),small_component(j,i),j=1,num_orbitals)
        end do 

        close(2)

        3030 FORMAT(A16,1x,1000(I3,I3,A2,A5,3X))
        3050 FORMAT (1X,1P,1000E16.8)
    end subroutine

    subroutine read_orbitals()
        logical :: rwfn_existence

        inquire(file="rwfn.out",exist=rwfn_existence)


        if (.not. rwfn_existence) then 
            stop 'rwfn.out not found - did you save it to something else? are you in the correct directory?'
        end if 

        open(1,file = "rwfn.out",form='unformatted',status='old')

        !allocate the matrices with the fixed dimensions
        !there is probably a more elegant way to do this
        !however, it is unlikely for there to be more than 50 relativistic orbitals
        allocate(large_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(small_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(radial_grid(NUM_RAD_POINTS_FIXED))
    
        large_component = 0.0d0
        small_component = 0.0d0 

        READ (1) header
        write(*,1020)

        !I need to add an exception here for when the number of orbitals is larger
        !than NUM_ORBITALS_FIXED
        !either that or have a dynamic allocation should another be detected...
        

        !sketch of new algorithm:
        !do while iostat == True 
        !read NP(II),NAK(II),E(II),MFJ(ii)
        !extend allocation to allow for  another orbital
        !read orbital
        !repeat

        DO II = 1,NUM_ORBITALS_FIXED
            READ (1, iostat = j) NP(II),NAK(II),E(II),MFJ(ii)
            if (j .ne. -1) then 
                num_orbitals = num_orbitals + 1
                write(*,1010) NP(II),NAK(II),E(II),MFJ(ii)
                READ (1) PZ, (large_component(ii,I),I=1,MFJ(ii)), (small_component(ii,I),I=1,MFJ(ii))       !, (small_component(I),I=1,MFJ)
                READ (1) (radial_grid(I),I=1,MFJ(ii))

                if     (NAK(ii) .ge. 0) then 
                    NH_ARRAY(ii) = angular_symbols_pos(NAK(ii))
                elseif (NAK(ii) .le. 0) then
                    NH_ARRAY(ii) = angular_symbols_neg(abs(NAK(ii)))
                else   
                    print*,'invalid kappa in orbital',ii
                    stop 
                end if 

            else
                write(*,1030)num_orbitals
                exit
            end if
        end do 

        close(1)


        num_points = maxval(MFJ) 

        1010 format(I6,1X,I6,1X,F12.6,1X,I7)
        1020 format('PrincN',1X,' Kappa',1X,'  Eigenv(Ha)',1X,'NPoints')
        1030 format(I4,1X,'relativistic orbials found - exiting grasp2018 read')
    end subroutine



end program 