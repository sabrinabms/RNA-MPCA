!***************************************************************!
! MULTIPLE PARTICLE COLLISION ALGORITHM (MPCA)                  !
!***************************************************************!
! Developed by: Eduardo Favero Pacheco da Luz (CAP/INPE)        !
! Modified by:  Sabrina Sambati                                 !    
!               Juliana Anochi                                  ! 
!               Reynier Hernandez Torres (CAP/INPE)             !
! Based on PCA (Wagner F. Sacco)                                !
!***************************************************************!

PROGRAM MPCA

    USE annTraining
    USE annGeneralization
    USE mpcaFunctions
    USE newtypes

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !**********************
    ! VARIABLES DEFINITION
    !**********************
    integer :: contD, iError
    integer :: ndimensions, i, j
    real :: harvest
    character(len = 50) :: stringArg, str, str0, str1

    type (Particle) :: oldParticle
    type (Particle) :: newParticle
    type (Particle) :: bestParticle    
    type (annConfig) :: config
    type (optionsMPCA) :: op
    type (statusMPCA) :: st

    integer (kind = 8) :: nClasses
    integer (kind = 8) :: nClassesValidation
    integer (kind = 8) :: nClassesGeneralization
    integer (kind = 8) :: nClassesTest
    integer (kind = 8) :: nInputs
    integer (kind = 8) :: nOutputs
    real (kind = 8) :: targetError
    integer (kind = 8) :: nEpochs
    integer :: loadWeightsBias
    logical :: haveValidation
    logical :: tryInitialArchitecture
    integer (kind = 8) :: lower_Hidden_Layers
    integer (kind = 8) :: lower_First_Hidden_Layer
    integer (kind = 8) :: lower_Second_Hidden_Layer
    integer (kind = 8) :: lower_Activation_Function
    real (kind = 8) :: lower_Alpha
    real (kind = 8) :: lower_Eta
    integer (kind = 8) :: upper_Hidden_Layers
    integer (kind = 8) :: upper_First_Hidden_Layer
    integer (kind = 8) :: upper_Second_Hidden_Layer
    integer (kind = 8) :: upper_Activation_Function
    real (kind = 8) :: upper_Alpha
    real (kind = 8) :: upper_Eta
    integer (kind = 8) :: initial_Hidden_Layers
    integer (kind = 8) :: initial_First_Hidden_Layer
    integer (kind = 8) :: initial_Second_Hidden_Layer
    integer (kind = 8) :: initial_Activation_Function
    real (kind = 8) :: initial_Alpha
    real (kind = 8) :: initial_Eta
    logical :: doStopMPCA

    real (kind = 8) :: value_to_reach
    integer :: particles_processor
    integer (kind = 8) :: maximum_nfe_mpca
    integer :: cycle_blackboard_mpca
    integer (kind = 8) :: nfe_exploitation_mpca
    real (kind = 8) :: lower_exploitation_mpca
    real (kind = 8) :: upper_exploitation_mpca
    integer :: type_probability_mpca
    logical :: verbose

    NAMELIST /content/&
        nClasses, &
        nClassesValidation, &
        nClassesGeneralization, &
        nClassesTest, &
        nInputs, &
        nOutputs, &
        targetError, &
        nEpochs, &
        loadWeightsBias, &
        haveValidation, &
        tryInitialArchitecture

    NAMELIST /bounds/&
        lower_Hidden_Layers, &
        upper_Hidden_Layers, &
        lower_First_Hidden_Layer, &
        upper_First_Hidden_Layer, &
        lower_Second_Hidden_Layer, &
        upper_Second_Hidden_Layer, &
        lower_Activation_Function, &
        upper_Activation_Function, &
        lower_Alpha, &
        upper_Alpha, &
        lower_Eta, &
        upper_Eta

    NAMELIST /initial/ &
        initial_Hidden_Layers, &
        initial_First_Hidden_Layer, &
        initial_Second_Hidden_Layer, &
        initial_Activation_Function, &
        initial_Alpha, &
        initial_Eta

    NAMELIST /algorithm_configuration/&
        value_to_reach, &
        particles_processor, &
        maximum_nfe_mpca, &
        cycle_blackboard_mpca, &
        nfe_exploitation_mpca, &
        lower_exploitation_mpca, &
        upper_exploitation_mpca, &
        type_probability_mpca, &
        verbose

    !******************
    ! INITIALIZING MPI
    !******************
    CALL MPI_Init(iError)
    CALL MPI_Comm_Size(MPI_COMM_WORLD, op % nProcessors, iError)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, op % iProcessor, iError)

    !*********************************************
    ! SETTING PARAMETERS
    ! 1 - Number of the experiment
    !*********************************************
    CALL getarg(1, stringArg)
    read (stringArg, '(I10)') op % iExperiment
    
    op % ndimensions = 6

    OPEN(1, FILE='./config/configuration.ini', status='OLD', action='READ')
        read(1, algorithm_configuration)
    CLOSE(1)

    op % emin = value_to_reach
    op % nParticlesProcessor = particles_processor
    op % maxNFE = maximum_nfe_mpca
    op % iCycleBlackboard = cycle_blackboard_mpca
    op % iterPerturbation = nfe_exploitation_mpca
    op % lo_small = lower_exploitation_mpca
    op % up_small = upper_exploitation_mpca
    op % typeProbability = type_probability_mpca
    op % verbose = verbose

    !Output files
    if (op % iProcessor < 10) then
        write (str0, '(I1)') op % iProcessor
    else
        write (str0, '(I2)') op % iProcessor
    endif

    if (op % iExperiment < 10) then
        write (str1, '(I1)') op % iExperiment
    else
        write (str1, '(I2)') op % iExperiment
    endif

    if (op % iProcessor == 0) then
        if (op % iExperiment == 1) then
            OPEN(UNIT = 2, FILE = './output/final.out', ACCESS = 'APPEND')
                write(2,*)'F.OBJ, N.CAM, N.NEUR C1, N.NEUR C2, F.ATIV, ALFA, ETA'
            CLOSE(2)
        endif
    endif

    allocate(oldParticle % solution(op % ndimensions))
    allocate(newParticle % solution(op % ndimensions))
    allocate(bestParticle % solution(op % ndimensions))

    OPEN(1, file='./config/configuration.ini', status='old', action='read')
        read(1, content)
        read(1, bounds)
        read(1, initial)
    CLOSE(1)

    config % nClasses = nClasses
    config % nClassesValidation = nClassesValidation
    config % nClassesGeneralization = nClassesGeneralization
    config % nClassesTest = nClassesTest
    config % nInputs = nInputs
    config % nOutputs = nOutputs
    config % targetError = targetError
    config % nEpochs = nEpochs
    config % loadWeightsBias = loadWeightsBias
    config % haveValidation = haveValidation
    config % tryInitialArchitecture = tryInitialArchitecture 

    op % lowerBound(1) = lower_Hidden_Layers
    op % lowerBound(2) = lower_First_Hidden_Layer
    op % lowerBound(3) = lower_Second_Hidden_Layer
    op % lowerBound(4) = lower_Activation_Function
    op % lowerBound(5) = lower_Alpha
    op % lowerBound(6) = lower_Eta
    op % upperBound(1) = upper_Hidden_Layers
    op % upperBound(2) = upper_First_Hidden_Layer
    op % upperBound(3) = upper_Second_Hidden_Layer
    op % upperBound(4) = upper_Activation_Function
    op % upperBound(5) = upper_Alpha
    op % upperBound(6) = upper_Eta

    !------------------------------------------------------------!
    !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
    !------------------------------------------------------------!
    allocate(config % x(config % nInputs, config % nClasses))
    allocate(config % y(config % nOutputs, config % nClasses))

    OPEN (2, file = './data/x.txt')
        do i = 1, config % nInputs
            read(2, *) (config % x(i, j), j = 1, config % nClasses)
        enddo
    CLOSE (2)

    OPEN (2, file = './data/y.txt')
        do I = 1, config % nOutputs
            read(2, *) (config % y(i, j), j = 1, config % nClasses)
        enddo
    CLOSE (2)

    if (config % haveValidation .eqv. .true.) then
       
        allocate(config % x_valid(config % nInputs, config % nClassesValidation))
        allocate(config % y_valid(config % nOutputs, config % nClassesValidation))

        OPEN (1, file = './data/y_valid.txt')
            do i = 1, config % nOutputs
                read(1, *) (config % y_valid(i, j), j = 1, config % nClassesValidation)
            enddo
        CLOSE (1)

        OPEN (2, file = './data/x_valid.txt')
            do i = 1, config % nInputs
                read(2, *) (config % x_valid(i, j), j = 1, config % nClassesValidation)
            enddo
        CLOSE (2)

    endif

    ! RANDOM NUMBER CONFIGURATION
    CALL init_random_seed(op)

    st % NFE = 0
    st % higherNFE = 0
    st % lastUpdate = 0
    st % totalNFE = 0
    st % iBest = 0
    st % flag = .false.
    st % bestObjectiveFunction = huge(0.D0)
    st % doStop = .false.
    doStopMPCA = .false.

    bestParticle % fitness = huge(0.D0)


    !*****************************
    ! CREATING INITIAL POPULATION
    !*****************************
    if (config % tryInitialArchitecture .eqv. .true.) then
            oldParticle % solution(1) = initial_Hidden_Layers
            oldParticle % solution(2) = initial_First_Hidden_Layer
            oldParticle % solution(3) = initial_Second_Hidden_Layer
            oldParticle % solution(4) = initial_Activation_Function
            oldParticle % solution(5) = initial_Alpha
            oldParticle % solution(6) = initial_Eta
    else
        do contD = 1, op % ndimensions
            CALL RANDOM_NUMBER(harvest)
            oldParticle % solution(contD) = (harvest &
            & * (op % upperBound(contD) - op % lowerBound(contD))) &
            & + op % lowerBound(contD)
        enddo
    endif

    oldParticle % fitness = neuralNetworkTraining(oldParticle % solution, op, st, config)
    st % NFE = st % NFE + 1

    if (oldParticle % fitness < bestParticle % fitness) then
        bestParticle = oldParticle
    endif

    !***************************************************************************
    ! PRINCIPAL LOOP
    !***************************************************************************
    ! do while ((st % higherNFE .LE. op % maxNFE / op % nProcessors) &
    !     .and. (st % NFE .LE. op % maxNFE / op % nProcessors) &
    !     .and. (st % totalNFE .LE. op % maxNFE) &
    !     .and. (.not. doStopMPCA))
    do while (st % NFE .LE. op % maxNFE / op % nProcessors)

        CALL Perturbation(oldParticle, newParticle, bestParticle, op, st, config)

        if (newParticle % fitness < oldParticle % fitness) then

            if (newParticle % fitness < bestParticle % fitness) then
                bestParticle = newParticle
            endif

            oldParticle = newParticle
            CALL Exploration(oldParticle, newParticle, bestParticle, op, st, config)
        
        else
        
            CALL Scattering(oldParticle, newParticle, bestParticle, op, st, config)
        
        endif

        if(st % doStop) then
            exit
        endif

        ! BLACKBOARD UPDATE
        if (((st % NFE - st % lastUpdate) > op % iCycleBlackboard) &
            & .AND. (st % higherNFE < op % maxNFE / op % nProcessors) &
            & .AND. ((st % higherNFE + op % iCycleBlackboard) < op % maxNFE / op % nProcessors)) &
            & then

            CALL blackboard(bestParticle, st % NFE, st % higherNFE, st % totalNFE, st % doStop, &
            & doStopMPCA, op, st)

            st % lastUpdate = st % NFE

        endif

    enddo

    !*************************
    ! FINAL BLACKBOARD UPDATE
    !*************************

    CALL blackboard(bestParticle, st % NFE, st % higherNFE, st % totalNFE, st % doStop, &
            & doStopMPCA, op, st)

    if (op % iProcessor == 0) then
        CALL copyFileBest(op, st)
    endif

    if (op % iProcessor == 0) then
        OPEN(UNIT = 20, FILE = './output/final.out', ACCESS = 'APPEND')
            write(20, '(ES14.6E2)',ADVANCE="NO") bestParticle % fitness
            write(20, '(I2)',ADVANCE="NO") ceiling(bestParticle % solution(1))
            write(20, '(I3)',ADVANCE="NO") ceiling(bestParticle % solution(2))

            if (ceiling(bestParticle % solution(1)) .GT. 1) then
                write(20, '(I3)',ADVANCE="NO") ceiling(bestParticle % solution(3))
            else
                write(20, '(I3)',ADVANCE="NO") 0
            endif

            write(20, '(I2)',ADVANCE="NO") ceiling(bestParticle % solution(4))
            write(20, '(ES14.6E2)',ADVANCE="NO") bestParticle % solution(5)
            write(20, '(ES14.6E2)',ADVANCE="NO") bestParticle % solution(6)
        CLOSE(20)

        if (op % verbose .eqv. .true.) then
            write(*,FMT="(A1,A,t25,I10,A,I10)",ADVANCE="NO") achar(13), &
                & "NFE (total): ", &
                & st % totalNFE, &
                & " of ", &
                & op % maxNFE
            write(*,*)

            CALL write_formatted('Best objective function value: ', 'bright normal', &
                real_to_string_scientific(bestParticle % fitness, 1, 4, 3), 'normal')
            CALL write_formatted('Number of hidden layers: ', 'bright normal', &
                integer_to_string(ceiling(bestParticle % solution(1)), 2), 'normal')
            CALL write_formatted('Neurons in hidden layer 1: ', 'bright normal', &
                integer_to_string(ceiling(bestParticle % solution(2)), 2), 'normal')

            if (ceiling(bestParticle % solution(1)) == 2) then
                CALL write_formatted('Neurons in hidden layer 2: ', 'bright normal', &
                    integer_to_string(ceiling(bestParticle % solution(3)), 2), 'normal')
            endif

            CALL write_formatted('Activation function: ', 'bright normal', &
                integer_to_string(ceiling(bestParticle % solution(4)), 2), 'normal')
            CALL write_formatted('Alpha: ', 'bright normal', &
                real_to_string(bestParticle % solution(5), 1, 4), 'normal')
            CALL write_formatted('Eta: ', 'bright normal', &
                real_to_string(bestParticle % solution(6), 1, 4), 'normal')
            CALL end_section('', 'normal')

        endif

    endif

    ! FINALIZING

    deallocate(oldParticle % solution)
    deallocate(newParticle % solution)
    deallocate(bestParticle % solution)
    deallocate(config % x)
    deallocate(config % y)

    if (config % haveValidation .eqv. .true.) then
        deallocate(config % x_valid)
        deallocate(config % y_valid)
    endif

    CALL MPI_Finalize(iError)

END PROGRAM MPCA