!***************************************************************!
! Optimization of ANN architecture by metaheuristic MPCA        !
!***************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)  !
! Modified by: Reynier Hernandez Torres (CAP/INPE)              !
! Updated: 24-Mar-2016                                          !
!***************************************************************!
MODULE annGeneralization
    USE newTypes

CONTAINS

    SUBROUTINE annBackpropagation(nExperiments, nProcessors)
        IMPLICIT NONE

        double precision :: bfitness, eta, alpha, fitness, efitness
        real(8) :: mse
        integer :: best, iExperiment, bExperiment, iProcessor, bProcessor, i, j, k
        integer, intent(in) :: nExperiments, nProcessors
        character(8) :: dummy
        character(10) :: str0, str1
        character(32) :: path, fString
        integer :: hiddenLayers, neuronsLayer(2), activationFunction
        INTEGER (kind = 8) :: nClasses
        INTEGER (kind = 8) :: nClassesValidation
        INTEGER (kind = 8) :: nClassesGeneralization
        INTEGER (kind = 8) :: nClassesActivation
        INTEGER (kind = 8) :: nInputs
        INTEGER (kind = 8) :: nOutputs
        REAL (kind = 8) :: targetError
        INTEGER (kind = 8) :: nEpochs
        integer :: loadWeightsBias
        LOGICAL :: haveValidation
        logical :: tryInitialArchitecture

        TYPE(annConfig) :: config

        NAMELIST /content/ nClasses, nClassesValidation,&
                nClassesGeneralization, &
                nClassesActivation, &
                nInputs, nOutputs, &
                targetError, nEpochs, &
                loadWeightsBias, &
                haveValidation, &
                tryInitialArchitecture
        !tryFixedConfiguration

        OPEN(1, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
            read(1, content)
        close(1)

        config % nClasses = nClassesGeneralization
        config % nInputs = nInputs
        config % nOutputs = nOutputs

        open (unit = 20, file = './output/final.out', status = 'old', action = 'read')
        read(20, '(A)') dummy

        bfitness = 1e+10
        do iExperiment = 1, nExperiments
            read(20, '(ES14.6E2,I2,I3,I3,I2,ES14.6E2,ES14.6E2)') &
                    efitness, &
                    hiddenLayers, &
                    neuronsLayer(1), &
                    neuronsLayer(2), &
                    activationFunction, &
                    eta, &
                    alpha

            if (efitness < bfitness) then
                bfitness = efitness
                bExperiment = iExperiment
                config % hiddenLayers = hiddenLayers
                config % neuronsLayer(1) = neuronsLayer(1)
                config % neuronsLayer(2) = neuronsLayer(2)
                config % activationFunction = activationFunction
            end if
        end do

        PRINT*, 'MELHOR RESULTADO:'
        PRINT*, 'Função objetivo MPCA:', bfitness
		PRINT*, 'Numero camadas', config % hiddenLayers
		PRINT*, 'N.Neuronios C.1:', config % neuronsLayer(1)
        if (config % hiddenLayers < 1) then
            PRINT*, 'N.Neuronios C.2:',  config % neuronsLayer(2)
        end if
		PRINT*, 'Funcao de ativacao',  config % activationFunction

        close(20)

        IF (bExperiment < 10) THEN
            WRITE (str1, '(I1)') bExperiment
        ELSE IF (bExperiment < 100) THEN
            WRITE (str1, '(I2)') bExperiment
        ELSE
            WRITE (str1, '(I3)') bExperiment
        END IF

        PRINT*, 'MELHOR RESULTADO: ./output/ann' // trim(str1) // '.best'

        open(12, FILE = './output/ann' // trim(str1) // '.best')
        open(13, file = './output/nn.best')

        read(12, '(A)') dummy
!        print*, dummy
        read(12, *) fitness
        write(13, *) fitness
!        PRINT*,fitness

        read(12, '(A)') dummy
!        PRINT*, dummy
        read(12, *) activationFunction
        write(13, '(I3)') activationFunction
!        PRINT*, activationFunction

        read(12, '(A)') dummy
!        print*, dummy
        read(12, *) hiddenLayers
        write(13, '(I3)') hiddenLayers
!        PRINT*, hiddenLayers

        read(12, '(A)') dummy
!        print*, dummy
        read(12, *) neuronsLayer(1)
        write(13, '(I3)') neuronsLayer(1)
!        PRINT*, neuronsLayer(1)

        if (config % hiddenLayers > 1) then
            read(12, '(A)') dummy
!            print*, dummy
            read(12, *) neuronsLayer(2)
            write(13, '(I3)') neuronsLayer(2)
!            PRINT*,  config % neuronsLayer(2)
        end if

        read(12, *) dummy
!        print*, dummy
        write(13, '(A)') dummy


        allocate(config % bh1(config % neuronsLayer(1)))
        allocate(config % bs(config % nOutputs))
        allocate(config % wh1(config % nInputs, config % neuronsLayer(1)))

        fString = '(   F11.5)'
        write(fString(2:4), '(I3)') config % neuronsLayer(1)

        DO i = 1, config % nInputs
            read(12, *) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
            write(13, fString) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
        ENDDO

        read(12, *) dummy
        write(13, '(A)') dummy

        read(12, *) (config % bh1(k), k = 1, config % neuronsLayer(1))
        write(13, fString) (config % bh1(k), k = 1, config % neuronsLayer(1))

        if (config % hiddenLayers == 2) then
            allocate(config % bh2(config % neuronsLayer(2)))
            allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
            allocate(config % ws(config % neuronsLayer(2), config % nOutputs))

            read(12, *) dummy
            write(13, '(A)') dummy
            fString = '(   F11.5)'
            write(fString(2:4), '(I3)') config % neuronsLayer(2)
            DO i = 1, config % neuronsLayer(1)
                read(12, *) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
                write(13, fString) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
            ENDDO

            read(12, *) dummy
            write(13, '(A)') dummy
            read(12, *) (config % bh2(k), k = 1, config % neuronsLayer(2))
            write(13, fString) (config % bh2(k), k = 1, config % neuronsLayer(2))

            read(12, *) dummy
            write(13, '(A)') dummy
            fString = '(   F11.5)'
            write(fString(2:4), '(I3)') config % nOutputs
            DO i = 1, config % neuronsLayer(2)
                read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
                write(13, fString) (config % ws(i, k), k = 1, config % nOutputs)
            ENDDO
        else
            read(12, *) dummy
            write(13, '(A)') dummy
            fString = '(   F11.5)'
            write(fString(2:4), '(I3)') config % nOutputs
            allocate(config % ws(config % neuronsLayer(1), config % nOutputs))

            DO i = 1, config % neuronsLayer(1)
                read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
                write(13, fString) (config % ws(i, k), k = 1, config % nOutputs)
            ENDDO
        end if

        read(12, *) dummy
        write(13, '(A)') dummy
        read(12, *) (config % bs(k), k = 1, config % nOutputs)
        write(13, fString) (config % bs(k), k = 1, config % nOutputs)

        write(*, *) 'Activating ANN'
        mse = neuralNetwork(config)

        write(13, '(A)') 'mse'
        write(13, '(ES14.6E2)') mse

        deallocate(config % bh1)
        deallocate(config % bs)
        deallocate(config % wh1)
        deallocate(config % ws)
        
        if (config % hiddenLayers == 2) then
            deallocate(config % bh2)
            deallocate(config % wh2)
        end if

        close(12)
        close(13)

END SUBROUTINE annBackpropagation


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REAL(8) FUNCTION neuralNetwork(config)

    IMPLICIT NONE

    TYPE(annConfig), intent(in) :: config

    !Variaveis usadas para o calculo da funcao objetivo
    double precision :: bestF
    double precision, allocatable, dimension(:,:) :: x_gen
    double precision, allocatable, dimension(:,:) :: y_gen
    double precision, allocatable, dimension(:) :: vs
    double precision, allocatable, dimension(:,:) :: ys
    double precision, allocatable, dimension(:,:) :: error
    double precision :: eqm
    double precision, allocatable, dimension(:) :: vh1, vh2
    double precision, allocatable, dimension(:) :: yh1, yh2

    double precision, parameter :: a = 1

    integer :: i, j, k ! variavel para controle do programa
    character(32) :: fString

    !Bias 1. camada oculta
    allocate(error(config % nOutputs, config % nClasses))
    allocate(x_gen(config % nInputs, config % nClasses))
    allocate(y_gen(config % nOutputs, config % nClasses))
    allocate(vh1(config % neuronsLayer(1)))
    allocate(yh1(config % neuronsLayer(1)))

    if (config % hiddenLayers == 2) then
        allocate(vh2(config % neuronsLayer(2)))
        allocate(yh2(config % neuronsLayer(2)))
    end if

    allocate(vs(config % nOutputs))
    allocate(ys(config % nOutputs, config % nClasses))

    !------------------------------------------------------------!
    !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
    !------------------------------------------------------------!
    OPEN (1, file = './data/y_gen.txt')
    DO I = 1, config % nOutputs
        READ(1, *) (y_gen(I, J), J = 1, config % nClasses)
    END DO
    CLOSE (1)

    OPEN (2, file = './data/x_gen.txt')
    DO I = 1, config % nInputs
        READ(2, *) (x_gen(I, J), J = 1, config % nClasses)
    END DO
    CLOSE (2)

    !----------------------------------------------------------------------!
    ! INICIO DA REDE: FEEDFORWARD
    !----------------------------------------------------------------------!

    DO i = 1, config % nClasses
        ! ATIVACAO CAMADA OCULTA 1
        vh1 = 0.d0
        vh1 = matmul(x_gen(:, i), config % wh1)
        vh1 = vh1 - config % bh1;
        select case(config % activationFunction)
        case (1) !LOGISTICA
            yh1 = 1.d0/(1.d0 + DEXP(-a * vh1))
        case (2) !TANGENTE
            yh1 = (1.d0 - DEXP(-vh1))/(1.d0 + DEXP(-vh1))
        case (3) !GAUSS
            yh1 = DEXP(-vh1)
        end select

        if (config % hiddenLayers == 2) then
            ! ATIVACAO CAMADA OCULTA 2
            vh2 = 0.d0
            vh2 = matmul(yh1, config % wh2)
            vh2 = vh2 - config % bh2;
            select case(config % activationFunction)
            case (1) !LOGISTICA
                yh2 = 1.d0/(1.d0 + DEXP(-a * vh2))
            case (2) !TANGENTE
                yh2 = (1.d0 - DEXP(-vh2))/(1.d0 + DEXP(-vh2))
            case (3) !GAUSS
                yh2 = DEXP(-vh2)
            end select
        end if

        ! ATIVACAO: CAMADA DE SAIDA
        vs = 0.d0
        if (config % hiddenLayers == 1) then
            vs = matmul(yh1, config % ws)
        else if (config % hiddenLayers == 2) then
            vs = matmul(yh2, config % ws)
        end if
        vs = vs - config % bs

        select case(config % activationFunction)
        case (1) !LOGISTICA
            ys(:, i) = 1.d0/(1.d0 + DEXP(-a * vs(:)))
        case (2) !TANGENTE
            ys(:, i) = (1.d0 - DEXP(-vs(:)))/(1.d0 + DEXP(-vs(:)))
        case (3) !GAUSS
            ys(:, i) = DEXP(-vs(:))
        end select

        ! CALCULO ERRO TREINAMENTO
        error(:, i) = y_gen(:, i) - ys(:, i)
!	write(*,*) y_gen(:, i), ys(:, i)

    ENDDO
    eqm = sum(error)
    eqm = (1.d0/(config % nClasses)) * eqm

    write(*,*) "Passou por aqui - result_ys.out"
    open(12, file = './output/result_ys.out')
    fString = '(   F11.5)'
    write(fString(2:4), '(I3)') config % nOutputs
    do i = 1, config % nOutputs
        write(12, fString) (ys(i, j), j = 1, config % nClasses)
    end do
    close(12)

    neuralNetwork = eqm

    deallocate(error)
    deallocate(x_gen)
    deallocate(y_gen)
    deallocate(vh1)
    deallocate(yh1)

    if (config % hiddenLayers == 2) then
        deallocate(vh2)
        deallocate(yh2)
    end if

    deallocate(vs)
    deallocate(ys)

END FUNCTION neuralNetwork

END MODULE annGeneralization
