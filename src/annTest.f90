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

        double precision, allocatable, dimension(:) :: vs
        double precision, allocatable, dimension(:) :: ys

        real(8) :: eqm
        double precision, allocatable, dimension(:,:) :: error
        real (8), allocatable, dimension(:) :: errorClass
        double precision, allocatable, dimension(:) :: vh1, vh2
        double precision, allocatable, dimension(:) :: yh1, yh2

        integer :: i, j, k ! variavel para controle do programa
        character(32) :: fString
        character (100) :: DUMMY, fileNameBestofBest

        allocate(config % x_test(config % nInputs, config % nClassesTest))
        allocate(config % y_test(config % nOutputs, config % nClassesTest))

        allocate(vs(config % nOutputs))
        allocate(ys(config % nOutputs))
        allocate(config % bs(config % nOutputs))

        allocate(error(config % nOutputs, config % nClassesTest))
        allocate(errorClass(config % nClassesTest))

        
        !------------------------------------------------------------!
        !LENDO OS PARAMETROS DO ARQUIVO DE ENTRADA
        !------------------------------------------------------------!
        OPEN (1, file = './data/y_test.txt')
        DO I = 1, config % nOutputs
            READ(1, *) (config % y_test(I, J), J = 1, config % nClassesTest)
        END DO
        CLOSE (1)

        OPEN (2, file = './data/x_test.txt')
        DO I = 1, config % nInputs
            READ(2, *) (config % x_test(I, J), J = 1, config % nClassesTest)
        END DO
        CLOSE (2)

        print*, " Number of test classes: ", config % nClassesTest

        open(12, file =  "./output/nn.best", STATUS = "old")

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
    open(11, file = './output/result_test.out')
    open(12, file = './output/errors_test.out')
!    fString = '(   F11.5)'
!    write(fString(2:4), '(I3)') config % nOutputs

    DO i = 1, config % nClassesTest
	
        !ACTIVATING HIDDEN LAYER 1

	    vh1 = matmul(config % x_test(:, i), config % wh1) - config % bh1
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
        error(:, i) = config % y_test(:, i) - ys
        errorClass(i) = sum(error(:, i), dim = 1)
        errorClass(i) = 0.5d0 * (errorClass(i)**2.d0)
        write(12, *)error

    ENDDO
    close(12)
    open(12, file = './output/eqm_test.out')
    eqm = sum(errorClass) / dfloat(config % nClassesTest)
    write(12,*)eqm
    close(12)


    END subroutine neuralNetworkTest

END MODULE annTest



!     neuralNetwork = eqm

!     deallocate(error)
!     deallocate(x_test)
!     deallocate(y_test)
!     deallocate(vh1)
!     deallocate(yh1)

!     if (config % hiddenLayers == 2) then
!         deallocate(vh2)
!         deallocate(yh2)
!     end if

!     deallocate(vs)
!     deallocate(ys)


!     !     real (8), allocatable, dimension(:) :: vs
!     !     real (8), allocatable, dimension(:,:) :: error
!     !     real (8), allocatable, dimension(:) :: yh1
!     !     real (8), allocatable, dimension(:) :: vh1
!     !     real (8), allocatable, dimension(:) :: errorClass
!     !     real (8), allocatable, dimension(:) :: yh2
!     !     real (8), allocatable, dimension(:) :: vh2
!     !     real (8), allocatable, dimension(:,:) :: y


!     !     ! real (8), parameter :: alphaObj = 0.1D0
!     !     ! real (8), parameter :: betaObj = 1.0D0
!     !     ! real (8), parameter :: p1 = 5.0e-8
!     !     ! real (8), parameter :: p2 = 5.0e-5

!     !     character (100) :: str0
!     !     character (100) :: str1
!     !     character (100) :: fString
!     !     character (100) :: dummy

!     !     integer :: i
!     !     integer :: j
!     !     integer :: k

!     !     !Allocating space for error variables
!     !     allocate(error(config % nOutputs, config % nClassesTest))
!     !     allocate(errorClass(config % nClassesTest))

!     !     ! Allocating space for v and y variables
!     !     allocate(vh1(config % neuronsLayer(1)))
!     !     allocate(yh1(config % neuronsLayer(1)))
!     !     if (config % hiddenLayers > 1) then
!     !         allocate(vh2(config % neuronsLayer(2)))
!     !         allocate(yh2(config % neuronsLayer(2)))
!     !     end if
!     !     allocate(vs(config % nOutputs))

!     !     neuralNetworkActivation = 0 ! inicializando variavel

!     !     !----------------------------------------------------------------------!
!     !     ! ACTIVATION
!     !     !----------------------------------------------------------------------!
!     !     allocate(y(config % nOutputs, config % nClassesTest))

!     !     fString = '(      F8.5)'
!     !     write(fString(2:7), '(I6)') config % nClassesTest

!     !     do i = 1, config % nClassesTest
!     !         ! ACTIVATING HIDDEN LAYER 1
!     !         vh1 = matmul(config % x(:, i), config % wh1) - config % bh1
!     !         fString = '(F8.5)'
!     !         write(fString(2:7), '(I6)') config % neuronsLayer(1)

!     !         yh1 = activation(vh1, config % activationFunction)

!     !         if (config % hiddenLayers == 1) then
!     !             vs = matmul(yh1, config % ws) - config % bs
!     !         end if

!     !         ! ACTIVATING HIDDEN LAYER 2
!     !         if (config % hiddenLayers == 2) then
!     !             vh2 = matmul(yh1, config % wh2) - config % bh2
!     !             yh2 = activation(vh2, config % activationFunction)
!     !             vs = matmul(yh2, config % ws) - config % bs
!     !         end if

!     !         ! ACTIVATING OUTPUT
!     !         y(:, i) = activation(vs, config % activationFunction)
!     !         error(:, i) = config % y(:, i) - y(:, i)
!     !     end do

!     !     neuralNetworkActivation = sum(error) / dfloat(config % nClassesTest)

!     !     fString = '(      F8.5)'
!     !     OPEN (2, file = './output/y_activation.txt')
!     !     DO i = 1, config % nOutputs
!     !         write(2, *) (y(i, j) , j = 1, config % nClassesTest)
!     !     END DO
!     !     CLOSE (2)

!     !     deallocate(error)
!     !     deallocate(errorClass)
!     !     deallocate(vh1)
!     !     deallocate(yh1)
!     !     if (config % hiddenLayers == 2) then
!     !         deallocate(vh2)
!     !         deallocate(yh2)
!     !     end if
!     !     deallocate(vs)
!     !     deallocate(y)


!     ! real(kind = 8) elemental function activation(v, afunction) result(y)
!     !     implicit none

!     !     real (kind = 8), intent(in) :: v
!     !     integer, intent(in) :: afunction

!     !     real (kind = 8), parameter :: a = 1

!     !     select case(afunction)
!     !     case (1) !LOGISTIC
!     !         y = 1.d0/(1.d0 + exp(-a * v))
!     !     case (2) !TANGENT
!     !         y = (1.d0 - exp(-v))/(1.d0 + exp(-v))
!     !     case (3) !GAUSS
!     !         y = exp(-v)
!     !     end select
!     ! end function activation

!     ! real(kind = 8) elemental function derivate(v, y, afunction) result(d)
!     !     implicit none
!     !     real (kind = 8), intent(in) :: v
!     !     real (kind = 8), intent(in) :: y
!     !     integer, intent(in) :: afunction
!     !     real (kind = 8), parameter :: a = 1

!     !     select case(afunction)
!     !     case (1)!LOGISTICA: e^-x / (1+e^-x)^2
!     !         d = ((a * exp(-a * v))/((1.d0 + exp(-a * v))**2.d0))
!     !     case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
!     !         d = (2 * exp(-v)) / ((1 + exp(-v))**2)
!     !     case (3)!GAUSS
!     !         d = -y/a
!     !     end select
!     ! end function derivate