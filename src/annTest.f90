!*******************************************************************!
! Optimization of the architecture of Artificial Neural Networks    !
! using Multi-Particle Collision Algorithm                          !
!*******************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)      !
! Modified by: Reynier Hernandez Torres (CAP/INPE)                  !
!*******************************************************************!

MODULE annTest
    use newTypes
    ! use foul
    use rnaFunctions

CONTAINS

    subroutine neuralNetworkTest(config)

        IMPLICIT NONE

        TYPE(annConfig) :: config
        ! TYPE(annConfig), intent(in) :: config

        double precision, allocatable, dimension(:,:) :: x_test
        double precision, allocatable, dimension(:,:) :: y_test
        double precision, allocatable, dimension(:) :: vs
        double precision, allocatable, dimension(:) :: ys
        double precision, allocatable, dimension(:) :: vh1, vh2
        double precision, allocatable, dimension(:) :: yh1, yh2

        double precision, allocatable, dimension(:,:) :: error
        real (8), allocatable, dimension(:) :: errorClass
        real(8) :: eqm


        integer :: i, j, k ! variavel para controle do programa
        character(32) :: fString
        character (100) :: DUMMY, fileNameBestofBest

        allocate(x_test(config % nInputs, config % nClassesTest))
        allocate(y_test(config % nOutputs, config % nClassesTest))

        allocate(vs(config % nOutputs))
        allocate(ys(config % nOutputs))
        allocate(config % bs(config % nOutputs))

        allocate(error(config % nOutputs, config % nClassesTest))
        allocate(errorClass(config % nClassesTest))

        
        !------------------------------------------------------------!
        !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
        !------------------------------------------------------------!
        OPEN (1, file = './datain/y_test.txt')
        DO I = 1, config % nOutputs
            READ(1, *) (y_test(I, J), J = 1, config % nClassesTest)
        END DO
        CLOSE (1)

        OPEN (2, file = './datain/x_test.txt')
        DO I = 1, config % nInputs
            READ(2, *) (x_test(I, J), J = 1, config % nClassesTest)
        END DO
        CLOSE (2)

        print*, " Number of test classes: ", config % nClassesTest

        open(12, file =  "./dataout/nn.best", STATUS = "old")

        read(12, *) dummy !Funcao objetivo do MPCA
        read(12, *) dummy
        read(12, *) dummy !Funcao ativação
        read(12, *) config % activationFunction
        print*, " Activation Function: ", config % activationFunction

        read(12, *) dummy !Número de camadas
        read(12, *) config % hiddenLayers
        read(12, *) dummy !Numero neuronios camada 1
        read(12, *) config % neuronsLayer(1)
        print*,' Number of hidden layers: ',  config %hiddenLayers
        print*,' Neurons in hidden layer 1: ',  config % neuronsLayer(1)
        
        if (config % hiddenLayers > 1) then
            read(12, *) dummy !Numero neuronios camada 1
            read(12, *) config %neuronsLayer(2)
            print*,' Neurons in hidden layer 2: ', config % neuronsLayer(2)
        endif

        allocate(vh1(config % neuronsLayer(1)))
        allocate(yh1(config % neuronsLayer(1)))
        allocate(config % bh1(config % neuronsLayer(1)))
        allocate(config % wh1(config % nInputs, config % neuronsLayer(1)))

        if (config % hiddenLayers > 1) then
            allocate(vh2(config % neuronsLayer(2)))
            allocate(yh2(config % neuronsLayer(2)))
            allocate(config % bh2(config % neuronsLayer(2)))
            allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
            allocate(config % ws(config % neuronsLayer(2), config % nOutputs))
        else
            allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
        endif

        read(12, *) dummy ! wh1
        do i = 1, config % nInputs
            read(12,*)(config % wh1(i,k), k = 1, config % neuronsLayer(1))
        enddo 

        read(12, *) dummy ! bh1
        read(12,*)(config % bh1(k), k = 1, config % neuronsLayer(1))

        if (config % hiddenLayers > 1) then

            read(12, *) dummy ! wh2
            do i = 1, config % neuronsLayer(1)
                read(12,*)(config % wh2(i,k), k = 1, config % neuronsLayer(2))
            enddo 

            read(12, *) dummy ! bh2
            read(12,*)(config % bh2(k), k = 1, config % neuronsLayer(2))

            read(12, *) dummy ! ws
            do i = 1, config % neuronsLayer(2)
                read(12,*)(config % ws(i,k), k = 1, config % nOutputs)
            enddo
        else
            read(12, *) dummy ! ws
            do i = 1, config % neuronsLayer(1)
                read(12,*)(config % ws(i,k), k = 1, config % nOutputs)
            enddo
        endif

        read(12, *) dummy ! bs
        read(12,*)(config % bs(k), k = 1, config % nOutputs) 
        close(12)
    !----------------------------------------------------------------------!
    ! INICIO DA REDE: FEEDFORWARD
    !----------------------------------------------------------------------!
    open(11, file = './dataout/result_test.out', access = 'append')
    open(12, file = './dataout/errors_test.out', access = 'append')
!    fString = '(   F11.5)'
!    write(fString(2:4), '(I3)') config % nOutputs

    DO i = 1, config % nClassesTest
	
        !ACTIVATING HIDDEN LAYER 1

	    vh1 = matmul(x_test(:, i), config % wh1) - config % bh1
        yh1 = activation(vh1, config % activationFunction)

        if (config % hiddenLayers == 1) then
        	vs = matmul(yh1, config % ws) - config % bs
        else
        	vh2 = matmul(yh1, config % wh2) - config % bh2
        	yh2 = activation(vh2, config % activationFunction)
        	vs = matmul(yh2, config % ws) - config % bs
        end if

        !ACTIVATING OUTPUT
        ys = activation(vs, config % activationFunction)
        write(11, *) ys !(ys(i, j), j = 1, config % nOutputs)

        !CALCULO ERRO TESTE
        error(:, i) = y_test(:, i) - ys
        errorClass(i) = sum(error(:, i), dim = 1)
        errorClass(i) = 0.5d0 * (errorClass(i)**2.d0)
        write(12, *)error

    ENDDO
    close(11)
    close(12)
    open(12, file = './dataout/eqm_test.out')
    eqm = sum(errorClass) / dfloat(config % nClassesTest)
    write(12,*)eqm
    close(12)


    END subroutine neuralNetworkTest

END MODULE annTest