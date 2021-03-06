!***************************************************************!
! Optimization of ANN architecture by metaheuristic MPCA        !
!***************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)  !
! Modified by: Reynier Hernandez Torres (CAP/INPE)              !
!***************************************************************!
MODULE annGeneralization
    USE newTypes
    USE rnaFunctions

CONTAINS

    REAL(8) FUNCTION neuralNetwork(config)

        IMPLICIT NONE

        TYPE(annConfig), intent(in) :: config

        !Variaveis usadas para o calculo da funcao objetivo
        double precision :: bestF
        double precision, allocatable, dimension(:,:) :: x_gen
        double precision, allocatable, dimension(:,:) :: y_gen
        double precision, allocatable, dimension(:) :: vs
        double precision, allocatable, dimension(:) :: ys
        real(8) :: eqm
        double precision, allocatable, dimension(:,:) :: error
        real (8), allocatable, dimension(:) :: errorClass
        double precision, allocatable, dimension(:) :: vh1, vh2
        double precision, allocatable, dimension(:) :: yh1, yh2

        double precision, parameter :: a = 1

        integer :: i, j, k ! variavel para controle do programa
        character(32) :: fString

        !Bias 1. camada oculta
        allocate(x_gen(config % nInputs, config % nClassesGeneralization))
        allocate(y_gen(config % nOutputs, config % nClassesGeneralization))
        allocate(vh1(config % neuronsLayer(1)))
        allocate(yh1(config % neuronsLayer(1)))
        allocate(error(config % nOutputs, config % nClassesGeneralization))
        allocate(errorClass(config % nClassesGeneralization))

        if (config % hiddenLayers == 2) then
            allocate(vh2(config % neuronsLayer(2)))
            allocate(yh2(config % neuronsLayer(2)))
        end if

        allocate(vs(config % nOutputs))
        allocate(ys(config % nOutputs))
        

        !------------------------------------------------------------!
        !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
        !------------------------------------------------------------!
        OPEN (1, file = './datain/y_gen.txt')
        DO I = 1, config % nOutputs
            READ(1, *) (y_gen(I, J), J = 1, config % nClassesGeneralization)
        END DO
        CLOSE (1)

        OPEN (2, file = './datain/x_gen.txt')
        DO I = 1, config % nInputs
            READ(2, *) (x_gen(I, J), J = 1, config % nClassesGeneralization)
        END DO
        CLOSE (2)

        !----------------------------------------------------------------------!
        ! INICIO DA REDE: FEEDFORWARD
        !----------------------------------------------------------------------!

        DO i = 1, config % nClassesGeneralization
        ! ACTIVATING HIDDEN LAYER 1
            vh1 = matmul(x_gen(:, i), config % wh1) - config % bh1
            yh1 = activation(vh1, config % activationFunction)

            if (config % hiddenLayers == 1) then
                vs = matmul(yh1, config % ws) - config % bs
            else
                vh2 = matmul(yh1, config % wh2) - config % bh2
                yh2 = activation(vh2, config % activationFunction)
                vs = matmul(yh2, config % ws) - config % bs
            end if

            ! ACTIVATING OUTPUT
            ys = activation(vs, config % activationFunction)

            ! CALCULO ERRO GENERALIZACAO
            error(:, i) = y_gen(:, i) - ys
            errorClass(i) = sum(error(:, i), dim = 1)
            errorClass(i) = 0.5d0 * (errorClass(i)**2.d0)

        ENDDO

        eqm = sum(errorClass) / dfloat(config % nClassesGeneralization)

    ! TRECHO COMENTADO PARA SER RESOLVIDO DEPOIS, NAO PODE SER APAGADO!!!!!
    !    open(12, file = './dataout/result_generalization.out')
    !    fString = '(   F11.5)'
    !    write(fString(2:4), '(I3)') config % nOutputs
    !    do i = 1, config % nOutputs
    !        write(12, fString) (ys(i, j), j = 1, config % nClassesGeneralization)
    !    end do
    !    close(12)

        neuralNetwork = eqm

        deallocate(error)
        deallocate(errorClass)
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
