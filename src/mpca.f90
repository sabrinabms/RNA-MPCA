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
    USE newTypes

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    !**********************
    ! VARIABLES DEFINITION
    !**********************
    integer :: iSeed, contD, contP, iCycleBlackboard, iError
    integer :: nSeed, un, istat, it, nDimensions, i, j
    integer, allocatable, DIMENSION(:) :: vSeed, vSeed2
	REAL (kind = 8), allocatable, DIMENSION(:) :: minB, maxB
    REAL (kind = 8), allocatable, DIMENSION(:) :: realSolution
    real :: harvest, rRandom
    CHARACTER(len = 50) :: string, str0, str1, fString, str
    character(len = 100) :: filename
    character(len = 100) :: format_string
    integer :: tmp_unit

    TYPE (Particle), allocatable, DIMENSION(:) :: oldParticle
    TYPE (Particle), allocatable, DIMENSION(:) :: newParticle
    TYPE (Particle), allocatable, DIMENSION(:) :: bestParticle
    TYPE (Particle) :: bestParticleProcessor
    TYPE (OptionsMPCA) :: op
    TYPE (annConfig) :: config
    TYPE (StatusMPCA) :: st

    INTEGER (kind = 8) :: nClasses
    INTEGER (kind = 8) :: nClassesValidation
    INTEGER (kind = 8) :: nClassesGeneralization
    INTEGER (kind = 8) :: nClassesTest
    INTEGER (kind = 8) :: nInputs
    INTEGER (kind = 8) :: nOutputs
    REAL (kind = 8) :: targetError
    INTEGER (kind = 8) :: nEpochs
    integer :: loadWeightsBias
    LOGICAL :: haveValidation
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
    logical :: doStopMPCA = .false.


    NAMELIST /content/ nClasses, &
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

    NAMELIST /bounds/ lower_Hidden_Layers, &
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

    NAMELIST /initial/ initial_Hidden_Layers, &
        initial_First_Hidden_Layer, &
        initial_Second_Hidden_Layer, &
        initial_Activation_Function, &
        initial_Alpha, &
        initial_Eta

    REAL (kind = 8) :: value_to_reach
    INTEGER :: particles_processor
    INTEGER (kind = 8) :: maximum_nfe_mpca
    INTEGER :: cycle_blackboard_mpca
    INTEGER (kind = 8) :: nfe_exploitation_mpca
    REAL (kind = 8) :: lower_exploitation_mpca
    REAL (kind = 8) :: upper_exploitation_mpca
    INTEGER :: type_probability_mpca
    LOGICAL :: verbose

    NAMELIST /algorithm_configuration/ value_to_reach, &
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
    CALL getarg(1, string)
    READ (string, '(I10)') op % iExperiment

    op % nDimensions = 6

    OPEN(1, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
    READ(1, algorithm_configuration)
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
    IF (op % iProcessor < 10) THEN
        WRITE (str0, '(I1)') op % iProcessor
    ELSE
        WRITE (str0, '(I2)') op % iProcessor
    END IF

    IF (op % iExperiment < 10) THEN
        WRITE (str1, '(I1)') op % iExperiment
    ELSE
        WRITE (str1, '(I2)') op % iExperiment
    END IF

    if (op % iProcessor == 0) then
        if (op % verbose .eqv. .true.) then
            str = integer_to_string(op % iExperiment, 3)
            CALL start_section2('Experiment ' // trim(str), 'normal')
        end if
        if (op % iExperiment == 1) then
            OPEN(UNIT = 20, FILE = './output/final.out', ACCESS = 'APPEND')
            WRITE(20,*)'F.OBJ, N.CAM, N.NEUR C1, N.NEUR C2, F.ATIV, ALFA, ETA'
            CLOSE(20)
        end if
    end if

    !Allocating space for dynamic variables
	allocate(minB(op % nDimensions))
	allocate(maxB(op % nDimensions))
    allocate(oldParticle(op % nParticlesProcessor))
    allocate(newParticle(op % nParticlesProcessor))
    allocate(bestParticle(op % nParticlesProcessor))
	allocate(st % minB(op % nDimensions))
	allocate(st % maxB(op % nDimensions))

    DO contP = 1, op % nParticlesProcessor
        allocate(oldParticle(contP) % solution(op % nDimensions))
        allocate(newParticle(contP) % solution(op % nDimensions))
        allocate(bestParticle(contP) % solution(op % nDimensions))
    END DO

	allocate(realSolution(op % nDimensions))

    OPEN(1, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
    read(1, content)
    read(1, bounds)
    read(1, initial)
    close(1)

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

    fString = '(F8.5)'
    write(fString(2:7), '(I6)') config % nClasses

    OPEN (2, file = './data/x.txt')
    DO i = 1, config % nInputs
        READ(2, *) (config % x(i, j), j = 1, config % nClasses)
    END DO
    CLOSE (2)

    OPEN (1, file = './data/y.txt')
    DO I = 1, config % nOutputs
        READ(1, *) (config % y(i, j), j = 1, config % nClasses)
    END DO
    CLOSE (1)

    if (config % haveValidation .eqv. .true.) then
        allocate(config % x_valid(config % nInputs, config % nClassesValidation))
        allocate(config % y_valid(config % nOutputs, config % nClassesValidation))

        write(fString(2:7), '(I6)') config % nClassesValidation

        OPEN (1, file = './data/y_valid.txt')
        DO i = 1, config % nOutputs
            READ(1, *) (config % y_valid(i, j), j = 1, config % nClassesValidation)
        END DO
        CLOSE (1)

        OPEN (2, file = './data/x_valid.txt')
        DO i = 1, config % nInputs
            READ(2, *) (config % x_valid(i, j), j = 1, config % nClassesValidation)
        END DO
        CLOSE (2)

    end if

    ! RANDOM NUMBER CONFIGURATION
    CALL init_random_seed(op)

    st % NFE = 0
    
    st % higherNFE = 0
    !st % it !estah sendo usado?
    st % lastUpdate = 0
    st % totalNFE = 0
    st % iBest = 0
    st % flag = .false.
    st % bestObjectiveFunction = huge(0.D0)
    ! st % fileUpdated ! estah sendo usado?
    st % doStop = .false.
    doStopMPCA = .false.

    DO contP = 1, op % nParticlesProcessor
        bestParticle(contP) % fitness = huge(0.D0)
    end do
    bestParticleProcessor % fitness = huge(0.D0)

    !*****************************
    ! CREATING INITIAL POPULATION
    !*****************************

    DO contP = 1, op % nParticlesProcessor
        ! print*, 'Entrou no loop : populacao inicial'

        if (contP == 1 .and. config % tryInitialArchitecture .eqv. .true.) then
            oldParticle(contP) % solution(1) = initial_Hidden_Layers
            oldParticle(contP) % solution(2) = initial_First_Hidden_Layer
            oldParticle(contP) % solution(3) = initial_Second_Hidden_Layer
            oldParticle(contP) % solution(4) = initial_Activation_Function
            oldParticle(contP) % solution(5) = initial_Alpha
            oldParticle(contP) % solution(6) = initial_Eta
        else
            do contD = 1, op % nDimensions
                CALL RANDOM_NUMBER(harvest)
                oldParticle(contP) % solution(contD) = (harvest &
                    & * (op % upperBound(contD) - op % lowerBound(contD))) &
                    & + op % lowerBound(contD)
            end do
        end if

        oldParticle(contP) % fitness = neuralNetworkTraining(oldParticle(contP) % solution, op, st, config)
        st % NFE = st % NFE + 1

        IF (oldParticle(contP) % fitness < bestParticle(contP) % fitness) THEN
            bestParticle(contP) = oldParticle(contP)
        END IF

    END DO

    DO contP = 1, op % nParticlesProcessor
        if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
            bestParticleProcessor = bestParticle(contP)
        end if
    END DO

    !***************************************************************************
    ! PRINCIPAL LOOP
    !***************************************************************************
    DO WHILE ((st % higherNFE .LE. op % maxNFE / op % nProcessors) &
        .and. (st % NFE .LE. op % maxNFE / op % nProcessors) &
        .and. (st % totalNFE .LE. op % maxNFE) &
        .and. (.not. doStopMPCA) &
        )
        !PARTICLES LOOP (OpenMP)
        DO contP = 1, op % nParticlesProcessor
            CALL Perturbation(oldParticle(contP), newParticle(contP), bestParticle(contP), op, st, config)
            IF (newParticle(contP) % fitness < oldParticle(contP) % fitness) THEN
                IF (newParticle(contP) % fitness < bestParticle(contP) % fitness) THEN
                    bestParticle(contP) = newParticle(contP)
                END IF
                oldParticle(contP) = newParticle(contP)

                CALL Exploration(oldParticle(contP), newParticle(contP), bestParticle(contP), op, st, config)
            ELSE
                CALL Scattering(oldParticle(contP), newParticle(contP), bestParticle(contP), op, st, config)
            END IF

            if(st % doStop) then
                exit
            end if
        END DO

        ! BLACKBOARD UPDATE
        IF (((st % NFE - st % lastUpdate) > iCycleBlackboard) &
            & .AND. (st % higherNFE < op % maxNFE / op % nProcessors) &
            & .AND. ((st % higherNFE + iCycleBlackboard) < op % maxNFE / op % nProcessors)) &
            & THEN
            DO contP = 1, op % nParticlesProcessor
                if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
                    bestParticleProcessor = bestParticle(contP)
                end if
            END DO

            call blackboard(bestParticleProcessor, st % NFE, st % higherNFE, st % totalNFE, st % doStop, &
            & doStopMPCA, op, st)

            DO contP = 1, op % nParticlesProcessor
                bestParticle(contP) = bestParticleProcessor
            END DO
            st % lastUpdate = st % NFE

        END IF
    END DO

    !*************************
    ! FINAL BLACKBOARD UPDATE
    !*************************
    DO contP = 1, op % nParticlesProcessor
        if (bestParticle(contP) % fitness < bestParticleProcessor % fitness) then
            bestParticleProcessor = bestParticle(contP)
        end if
    END DO

    call blackboard(bestParticleProcessor, st % NFE, st % higherNFE, st % totalNFE, st % doStop, &
            & doStopMPCA, op, st)

    DO contP = 1, op % nParticlesProcessor
        bestParticle(contP) = bestParticleProcessor
    END DO

    IF (op % iProcessor == 0) THEN
        call copyFileBest(op, st)
    END IF

    IF (op % iProcessor == 0) THEN
        OPEN(UNIT = 20, FILE = './output/final.out', ACCESS = 'APPEND')
        write(20, '(ES14.6E2)',ADVANCE="NO") bestParticleProcessor % fitness
        write(20, '(I2)',ADVANCE="NO") ceiling(bestParticleProcessor % solution(1))
        write(20, '(I3)',ADVANCE="NO") ceiling(bestParticleProcessor % solution(2))
        IF (ceiling(bestParticleProcessor % solution(1)) .GT. 1) THEN
            write(20, '(I3)',ADVANCE="NO") ceiling(bestParticleProcessor % solution(3))
        ELSE
            write(20, '(I3)',ADVANCE="NO") 0
        END IF
        write(20, '(I2)',ADVANCE="NO") ceiling(bestParticleProcessor % solution(4))
        write(20, '(ES14.6E2)',ADVANCE="NO") bestParticleProcessor % solution(5)
        write(20, '(ES14.6E2)',ADVANCE="NO") bestParticleProcessor % solution(6)
        CLOSE(20)

        IF (op % verbose .eqv. .true.) THEN
            write(*,FMT="(A1,A,t25,I10,A,I10)",ADVANCE="NO") achar(13), &
                & "NFE (total): ", &
                & st % totalNFE, &
                & " of ", &
                & op % maxNFE
            write(*,*)

            call write_formatted('Best objective function value: ', 'bright normal', &
                real_to_string_scientific(bestParticleProcessor % fitness, 1, 4, 3), 'normal')
            call write_formatted('Number of hidden layers: ', 'bright normal', &
                integer_to_string(ceiling(bestParticleProcessor % solution(1)), 2), 'normal')
            call write_formatted('Neurons in hidden layer 1: ', 'bright normal', &
                integer_to_string(ceiling(bestParticleProcessor % solution(2)), 2), 'normal')
            IF (ceiling(bestParticleProcessor % solution(1)) == 2) THEN
                call write_formatted('Neurons in hidden layer 2: ', 'bright normal', &
                    integer_to_string(ceiling(bestParticleProcessor % solution(3)), 2), 'normal')
            END IF

            call write_formatted('Activation function: ', 'bright normal', &
                integer_to_string(ceiling(bestParticleProcessor % solution(4)), 2), 'normal')

            call write_formatted('Alpha: ', 'bright normal', &
                real_to_string(bestParticleProcessor % solution(5), 1, 4), 'normal')
            call write_formatted('Eta: ', 'bright normal', &
                real_to_string(bestParticleProcessor % solution(6), 1, 4), 'normal')

            call end_section('', 'normal')
        END IF
    END IF

    ! FINALIZING

    deallocate(oldParticle)
    deallocate(newParticle)
	deallocate(minB)
	deallocate(maxB)
    deallocate(bestParticle)
	deallocate(realSolution)
    deallocate(config % x)
    deallocate(config % y)
    if (config % haveValidation .eqv. .true.) then
        deallocate(config % x_valid)
        deallocate(config % y_valid)
    end if

    CALL MPI_Finalize(iError)

END PROGRAM MPCA

SUBROUTINE init_random_seed(op)
    use newTypes

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    TYPE (OptionsMPCA), intent(in) :: op

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    seed = seed * (op % iProcessor + 1)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
END SUBROUTINE