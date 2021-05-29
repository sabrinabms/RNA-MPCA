PROGRAM MAIN_GENERALIZATION
    USE annGeneralization
    USE newTypes

    implicit none
    
    character(16) :: string
    integer :: nProcessors, nExperiments
    
    CALL getarg(1, string)
    READ (string, '(I10)') nExperiments
    CALL getarg(2, string)
    READ (string, '(I10)') nProcessors
    print*, '************************************************************'
    print*, 'Num. exps: ', nExperiments
    print*, 'Num. procs: ', nProcessors
    print*, '************************************************************'
    
    call annBackpropagation(nExperiments, nProcessors)
    
END PROGRAM MAIN_GENERALIZATION
