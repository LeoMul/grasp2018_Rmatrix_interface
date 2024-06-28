program read 
    implicit none

    integer :: i,j,II,stat
    real*8 :: x
    character*6 :: header

    integer,parameter :: NUM_ORBITALS_FIXED   = 50
    integer,parameter :: NUM_RAD_POINTS_FIXED = 500

    integer :: NP(NUM_ORBITALS_FIXED),NAK(NUM_ORBITALS_FIXED),MFJ(NUM_ORBITALS_FIXED)
    real*8 :: E(NUM_ORBITALS_FIXED)
    real*8,allocatable :: large_component(:,:),small_component(:,:),radial_grid(:)
    real*8 ::PZ

    real*8,allocatable :: generalised_occupation(:)

    integer ::num_orbitals, num_points 

    logical :: end_of_file

    write(*,100)
    print*,'Entering routine read_orbitals'
    call read_orbitals()    
    
    write(*,100)
    print*,'Entering routine read_occupation_numbers'

    call read_occupation_numbers(num_orbitals)
    
    write(*,100)
    print*,'Entering routine write_out'
    call write_out()
    
    write(*,100)

    100 format('------------------------------------')

    contains

    subroutine write_out()

        open(2,file='TARGET.INP')

        write(2,3040) num_orbitals,num_points
        WRITE (2,3050) (generalised_occupation(i),i=1,num_orbitals)
        WRITE (2,3050) (radial_grid(I),I=1,num_points)

        do i = 1,num_orbitals
            WRITE (2,3060) NP(i),NAK(i)
            WRITE (2,3050) (large_component(I,j),j=1,num_points)
            WRITE (2,3060) NP(i),NAK(i)
            WRITE (2,3050) (small_component(I,j),j=1,num_points)
        end do 

        close(2)

        3040 FORMAT (1X,2I7)
        3050 FORMAT (1X,1P,4E16.8)
        3060 FORMAT (1X,I4,2X,I4)
    end subroutine

    subroutine read_orbitals()
        open(1,file = "rwfn.out",form='unformatted',status='old')

        allocate(large_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(small_component(NUM_ORBITALS_FIXED,NUM_RAD_POINTS_FIXED))
        allocate(radial_grid(NUM_RAD_POINTS_FIXED))
    
        large_component = 0.0d0
        small_component = 0.0d0 

        READ (1) header
        write(*,1020)

        DO II = 1,NUM_ORBITALS_FIXED
            READ (1, iostat = j) NP(II),NAK(II),E(II),MFJ(ii)
            if (j .ne. -1) then 
                num_orbitals = num_orbitals + 1
                write(*,1010) NP(II),NAK(II),E(II),MFJ(ii)
                READ (1) PZ, (large_component(ii,I),I=1,MFJ(ii)), (small_component(ii,I),I=1,MFJ(ii))       !, (small_component(I),I=1,MFJ)
                READ (1) (radial_grid(I),I=1,MFJ(ii))
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

    subroutine read_occupation_numbers(num_orbitals)
        
        integer,intent(in) :: num_orbitals

        logical :: finished

        character*5 :: dummy_character
        real*8 :: dummy_real(5)
        !real*8, allocatable :: generalised_occupation(:)

        character*20 :: string_detector
        integer :: counter ,ii,jj
        integer :: subshell_mention_counter
        integer , parameter :: max_read = 100000
        open(3,file='rmcdhf.sum')

        finished=.false.
        counter = 0
        subshell_mention_counter = 0
        allocate(generalised_occupation(num_orbitals))

        do while (finished .eqv. .false.)
            counter = counter + 1

            read(3, *) string_detector
            !print*,string_detector

            if (string_detector .eq. 'Subshell') then 
                subshell_mention_counter = subshell_mention_counter + 1
            end if 

            if (subshell_mention_counter .gt. 1) then 
                exit
            end if 

            if (counter .gt. max_read) then 
                exit
            end if 

        end do 

        do ii = 1,num_orbitals
            read(3,*) dummy_character ,  (dummy_real(jj),jj=1,5) ,generalised_occupation(ii)
        end do 
        
        close(3)
    end subroutine

end program 