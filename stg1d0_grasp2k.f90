program read 
    !stg1d0_grasp2k.f90 
    !gfortran stg1d0_grasp2k.f90 -o stg1d0_grasp2k.x
    !author: Leo Patrick Mulholland, QUB PhD Student, ARC
    
    !reads in:
    !rwfn.out, rmcdhf.sum from grasp2018 runs
    !into the format employed by the DARC R-matrix codes (Norrington, Berrington, Ballance)
    
    !Is this a particularly elegant code? Not really. Does it work? Yes.
    !This code serves really as a proof of concept that it is quite easy
    !to employ the large atomic models used by other atomic physicists
    !to produce R-matrix data.
    !I hope to produce a similar code for fac, jac

    !In principle, their large models can be used to systematically produce good orbitals
    !With large amounts of CI in the eigenvectors.
    !If calculation size becomes a problem  (R-matrix is O(N^3) afterall), the important 
    !CSFS can be isolaed and used in conjunction with the fully optimised orbitals.

    !A point to consider is the implementation of pseudo-correlation-orbitals in the main
    !version of grasp2018. There are orbitals with large principal quantum numberincluded to 
    !help convergence of lower ones. The problem is these orbitals tend to have the incorrect
    !number of nodes, as this restriction is not enforced on them. Only the so-called spectroscopic
    !orbitals have this constraint. I have not tested this in R-matrix - but imagine results involving
    !excitation to those orbitals would probably be a little unphysical in terms of Collision strength
    !and so forth. 

    implicit none

    !Parameters to be onbtained.
    integer ::num_orbitals, num_points 

    !allocatable arrays.
    real*8 ,allocatable :: large_component(:,:),small_component(:,:),radial_grid(:)
    real*8 ,allocatable :: generalised_occupation(:)
    integer,allocatable :: princ_n(:),kappa(:),specific_num_pts(:)
    real*8 ,allocatable :: orbital_eigenvalue(:)
    character*2,allocatable:: angular_string(:)

    !subshell labels
    character*2 :: angular_symbols_neg(5) = (/'s ','p ','d ','f ','g ' /)
    character*2 :: angular_symbols_pos(4) = (/'p-','d-','f-','g-' /)

    !main function

    write(*,100)
    print*,'Entering routine read_orbitals'
    call read_orbitals()    
    
    write(*,100)
    print*,'Entering routine read_occupation_numbers'

    call read_occupation_numbers()
    write(*,100)
    print*,'Entering routine print_out_info'
    write(*,*)
    call print_out_info()
    write(*,*)
    write(*,100)

    print*,'Entering routine write_out_targetdotinp'
    call write_out_targetdotinp()
    write(*,100)
    print*,'STOPPING NORMALLY'

    100 format('----------------------------------------------------')

    contains



    subroutine get_num_orbitals_num_points()

        !This routine scans the rwfn.out file and calculates number of orbitals.
        !This results in a double read, but is porbably more efficient than 
        !reallocating large arrays.
        !It also sets the maximum number of radial points needed

        !This subroutine basically does a read through of rwfn.out and sets the dimensions.
        !This allows for the complete avoidance of preallocated arrays.

        integer :: line_counter
        integer :: iostat_store

        integer :: max_num_points

        !dummy arguments for reading
        integer :: a,b 
        real*8 :: c
        
        max_num_points = 0

        open(1,file = "rwfn.out",form='unformatted',status='old')

        !READ (1, iostat = iostat_store) 

        !so, first read the first two lines.
        
        line_counter = 2
        READ (1)
        READ (1, iostat = iostat_store) a,b,c,num_points
        
        !the do loop will read the next two lines, and then the header of the next orbital.
        !by having the header of the next orbital at the end of the looped routine,
        !it gaurantees an exit when the file is done.


        do while (iostat_store .ne. -1) 
            line_counter = line_counter + 3

            if (num_points .gt. max_num_points) then 
                max_num_points = num_points
            end if 

            READ (1)
            READ (1)
            READ (1, iostat = iostat_store) a,b,c,num_points

        end do 
        
        num_orbitals = (line_counter - 1)/3
        num_points = max_num_points

        write(*,1030) num_orbitals
        write(*,1040) max_num_points

        close(1)

        1030 format(I6,1X,'relativistic orbials found. ')
        1040 format(I6,' is the max radial grid length.')

    end subroutine

    subroutine read_orbitals()
        !reads the orbitals from the grasp2018 format
        !first checks that the rwfn file exists

        !then calls the subroutine get_num_orbitals_num_points() to set the dimensions.

        logical :: rwfn_existence
        integer :: iter_orbs
        integer :: iostat_store

        integer :: i
        real*8 ::PZ
        
        inquire(file="rwfn.out",exist=rwfn_existence)


        if (.not. rwfn_existence) then 
            stop 'rwfn.out not found - did you save it to something else? are you in the correct directory?'
        end if 

        !get dimensions
        call get_num_orbitals_num_points()

        !allocate everythng
        allocate(large_component(num_orbitals,num_points))
        allocate(small_component(num_orbitals,num_points))
        allocate(radial_grid(num_points))
        allocate(generalised_occupation(num_orbitals))
        allocate(princ_n(num_orbitals))
        allocate(kappa(num_orbitals))
        allocate(specific_num_pts(num_orbitals))
        allocate(orbital_eigenvalue(num_orbitals))
        allocate(angular_string(num_orbitals))

        !intiialize at zero
        large_component = 0.0d0
        small_component = 0.0d0 

        open(1,file = "rwfn.out",form='unformatted',status='old')

        READ (1) !read header

        !loop over rwfn.out, reading the format exactly as it is written in grasp2018.
        DO iter_orbs = 1,num_orbitals

            READ (1, iostat = iostat_store) & 
            princ_n(iter_orbs),kappa(iter_orbs),orbital_eigenvalue(iter_orbs),specific_num_pts(iter_orbs)

            READ (1) PZ, (large_component(iter_orbs,I),I=1,specific_num_pts(iter_orbs)),&
             (small_component(iter_orbs,I),I=1,specific_num_pts(iter_orbs))       

            READ (1) (radial_grid(I),I=1,specific_num_pts(iter_orbs))
            
            !this logic just generates the angular l(bar or no bar) label for each orbital.
            if ((kappa(iter_orbs).lt.(len(angular_symbols_pos)+1)) & 
                .OR. (kappa(iter_orbs) .ge. -(len(angular_symbols_pos)+1))) then
                !this logic just stops it running off the edge of either of the angular labels arrays

                if     (kappa(iter_orbs) .ge. 0) then 
                    angular_string(iter_orbs) = angular_symbols_pos(kappa(iter_orbs))
                elseif (kappa(iter_orbs) .le. 0) then
                    angular_string(iter_orbs) = angular_symbols_neg(abs(kappa(iter_orbs)))
                else   
                    !in principal this should never be reached. only if kappa==0
                    print*,'invalid kappa in orbital',iter_orbs
                    stop 
                end if 
            else
                !in the event I didn't put enough into the angular labels, come here and just **.
                angular_string(iter_orbs) = '**'
            end if 

        END DO 

        close(1)


    end subroutine


    subroutine read_occupation_numbers()
        !This routine loops through mcdhf.sum and gets the occupation numbers.
        
        !This file is in a formatted structure.
        !The file has the word 'Subshell' twice before the occupations are printed.

        !This routine this loops until Subshell is found twice, and then the occupations are 
        !found easily based on the table formatting from grasp2k.

        logical :: finished
        character*5 :: dummy_character
        real*8 :: dummy_real(5)
        character*20 :: string_detector
        integer :: counter ,ii,jj
        integer :: subshell_mention_counter
        
        !in principal this should be replaced with an iostat routine
        integer , parameter :: max_read = 100000

        logical :: sum_existence
            
        inquire(file="rmcdhf.sum",exist=sum_existence)

        if (.not. sum_existence) then 
            stop 'rmcdhf.sum not found - did you save it to something else? are you in the correct directory?'
        end if 
        open(3,file='rmcdhf.sum')


        finished=.false.
        counter = 0
        subshell_mention_counter = 0

        do while (finished .eqv. .false.)
            counter = counter + 1
            read(3, *) string_detector
            !generalised occ is after the second mention of Subshell
            !in the 2018 version of grasp (Fischer version)
            if (string_detector .eq. 'Subshell') then 
                subshell_mention_counter = subshell_mention_counter + 1
            end if 

            if (subshell_mention_counter .gt. 1) then 
                exit
            end if 

            if (counter .gt. max_read) then 
                stop 'max read reached, check rmcdhf.sum file'
                exit
            end if 

        end do 

        !reading in the generalised occupation numbers.
        !dummy character, dummy reals to just read in the occupation
        do ii = 1,num_orbitals
            read(3,*) dummy_character ,  (dummy_real(jj),jj=1,5) ,generalised_occupation(ii)
        end do 
        
   

        close(3)
        !sum of generalised occupation numbers should be #electrons.
        write(*,1010) sum(generalised_occupation)
        1010 format(' Total number of electrons is: ',F12.6)
    end subroutine

    subroutine print_out_info()
        !Writes out information about the debugging.
        integer :: ii 
        write(*,1020)
        do ii = 1,num_orbitals
            write(*,1010) princ_n(ii),angular_string(ii),kappa(ii),orbital_eigenvalue(ii) &
                ,specific_num_pts(ii),generalised_occupation(ii)
        end do 
        


        1020 format('   Subshell',1X,'Kappa',5X,'Eigenv(Ha)',2X,'NPoints',4X,'GenOcc')
        1010 format(3X,3X,I3,A2,3X,I3,3X,ES12.6,4X,I5,F10.4)

    end subroutine


    subroutine write_out_targetdotinp()
        
        integer :: iter_orbs
        integer :: iter_grid

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

        open(2,file='TARGET.INP')
        print*,'Writing to TARGET.INP'

        WRITE (2,3040) num_orbitals,num_points
        WRITE (2,3050) (generalised_occupation(iter_orbs),iter_orbs=1,num_orbitals)
        WRITE (2,3050) (radial_grid(iter_grid),iter_grid=1,num_points)

        do iter_orbs = 1,num_orbitals
            WRITE (2,3060) princ_n(iter_orbs),kappa(iter_orbs)
            WRITE (2,3050) (large_component(iter_orbs,iter_grid),iter_grid=1,num_points)
            WRITE (2,3060) princ_n(iter_orbs),kappa(iter_orbs)
            WRITE (2,3050) (small_component(iter_orbs,iter_grid),iter_grid=1,num_points)
        end do 

        close(2)

        3040 FORMAT (1X,2I7)
        3050 FORMAT (1X,1P,4E16.8)
        3060 FORMAT (1X,I4,2X,I4)

    end subroutine

end program 