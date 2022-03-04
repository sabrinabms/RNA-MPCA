!*******************************************************************!
! Optimization of the architecture of Artificial Neural Networks    !
! using Multi-Particle Collision Algorithm                          !
!*******************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)      !
! Modified by: Reynier Hernandez Torres (CAP/INPE)                  !
!*******************************************************************!

MODULE annActivation
    use newTypes
    ! use foul

CONTAINS

REAL(8) FUNCTION neuralNetworkActivation(fileNameBestofBest, config)

    IMPLICIT NONE

    TYPE(annConfig), intent(in) :: config

    double precision, allocatable, dimension(:,:) :: x_test
    double precision, allocatable, dimension(:,:) :: y_test
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
    allocate(x_test(config % nInputs, config % nClassesGeneralization))
    allocate(y_test(config % nOutputs, config % nClassesGeneralization))
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
    OPEN (1, file = './data/y_test.txt')
    DO I = 1, config % nOutputs
        READ(1, *) (y_test(I, J), J = 1, config % nClassesTest)
    END DO
    CLOSE (1)

    OPEN (2, file = './data/x_test.txt')
    DO I = 1, config % nInputs
        READ(2, *) (x_test(I, J), J = 1, config % nClassesTest)
    END DO
    CLOSE (2)

    !
    ! Lendo os melhores parametros achados pelo MPCA
    !
    open(12, fileNameBestofBest, STATUS = "old")
    read(12, *) dummy !Valor funcao objetivo
    read(12, *) activationFunction
    read(12, *) hiddenLayers
    read(12, *) neuronsLayer(1)
    print*,' Number of hidden layers: ',  hiddenLayers
    print*,' Neurons in hidden layer 1: ',  neuronsLayer(1)
    
    if (hiddenLayers > 1) then
        read(12, *) neuronsLayer(2)
        print*,' Neurons in hidden layer 2: ', neuronsLayer(2)
    endif

    ! PARAMOS AQUI!!!
    ! ESTAMOS LENDO O ARQUIVO QUE CONTEM OS MELHORES PARAMETROS E ALTERANDO O CODIGO DE ACORDO COM ROTINA
    ! DE GENERALIZACAO, TEMOS QUE TOMAR CUIDADO PORQUE NESSE CODIGO ESTAMOS USANDO ERRONEAMENTO O CONFIG
    !----------------------------------------------------------------------!
    ! INICIO DA REDE: FEEDFORWARD
    !----------------------------------------------------------------------!

    DO i = 1, config % nClassesTest
	! ACTIVATING HIDDEN LAYER 1
	vh1 = matmul(x_test(:, i), config % wh1) - config % bh1
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
        error(:, i) = y_test(:, i) - ys
        errorClass(i) = sum(error(:, i), dim = 1)
        errorClass(i) = 0.5d0 * (errorClass(i)**2.d0)

    ENDDO

    eqm = sum(errorClass) / dfloat(config % nClassesGeneralization)

! TRECHO COMENTADO PARA SER RESOLVIDO DEPOIS, NAO PODE SER APAGADO!!!!!
!    open(12, file = './output/result_generalization.out')
!    fString = '(   F11.5)'
!    write(fString(2:4), '(I3)') config % nOutputs
!    do i = 1, config % nOutputs
!        write(12, fString) (ys(i, j), j = 1, config % nClassesGeneralization)
!    end do
!    close(12)

    neuralNetwork = eqm

    deallocate(error)
    deallocate(x_test)
    deallocate(y_test)
    deallocate(vh1)
    deallocate(yh1)

    if (config % hiddenLayers == 2) then
        deallocate(vh2)
        deallocate(yh2)
    end if

    deallocate(vs)
    deallocate(ys)


    !     real (8), allocatable, dimension(:) :: vs
    !     real (8), allocatable, dimension(:,:) :: error
    !     real (8), allocatable, dimension(:) :: yh1
    !     real (8), allocatable, dimension(:) :: vh1
    !     real (8), allocatable, dimension(:) :: errorClass
    !     real (8), allocatable, dimension(:) :: yh2
    !     real (8), allocatable, dimension(:) :: vh2
    !     real (8), allocatable, dimension(:,:) :: y


    !     ! real (8), parameter :: alphaObj = 0.1D0
    !     ! real (8), parameter :: betaObj = 1.0D0
    !     ! real (8), parameter :: p1 = 5.0e-8
    !     ! real (8), parameter :: p2 = 5.0e-5

    !     character (100) :: str0
    !     character (100) :: str1
    !     character (100) :: fString
    !     character (100) :: dummy

    !     integer :: i
    !     integer :: j
    !     integer :: k

    !     !Allocating space for error variables
    !     allocate(error(config % nOutputs, config % nClassesTest))
    !     allocate(errorClass(config % nClassesTest))

    !     ! Allocating space for v and y variables
    !     allocate(vh1(config % neuronsLayer(1)))
    !     allocate(yh1(config % neuronsLayer(1)))
    !     if (config % hiddenLayers > 1) then
    !         allocate(vh2(config % neuronsLayer(2)))
    !         allocate(yh2(config % neuronsLayer(2)))
    !     end if
    !     allocate(vs(config % nOutputs))

    !     neuralNetworkActivation = 0 ! inicializando variavel

    !     !----------------------------------------------------------------------!
    !     ! ACTIVATION
    !     !----------------------------------------------------------------------!
    !     allocate(y(config % nOutputs, config % nClassesTest))

    !     fString = '(      F8.5)'
    !     write(fString(2:7), '(I6)') config % nClassesTest

    !     do i = 1, config % nClassesTest
    !         ! ACTIVATING HIDDEN LAYER 1
    !         vh1 = matmul(config % x(:, i), config % wh1) - config % bh1
    !         fString = '(F8.5)'
    !         write(fString(2:7), '(I6)') config % neuronsLayer(1)

    !         yh1 = activation(vh1, config % activationFunction)

    !         if (config % hiddenLayers == 1) then
    !             vs = matmul(yh1, config % ws) - config % bs
    !         end if

    !         ! ACTIVATING HIDDEN LAYER 2
    !         if (config % hiddenLayers == 2) then
    !             vh2 = matmul(yh1, config % wh2) - config % bh2
    !             yh2 = activation(vh2, config % activationFunction)
    !             vs = matmul(yh2, config % ws) - config % bs
    !         end if

    !         ! ACTIVATING OUTPUT
    !         y(:, i) = activation(vs, config % activationFunction)
    !         error(:, i) = config % y(:, i) - y(:, i)
    !     end do

    !     neuralNetworkActivation = sum(error) / dfloat(config % nClassesTest)

    !     fString = '(      F8.5)'
    !     OPEN (2, file = './output/y_activation.txt')
    !     DO i = 1, config % nOutputs
    !         write(2, *) (y(i, j) , j = 1, config % nClassesTest)
    !     END DO
    !     CLOSE (2)

    !     deallocate(error)
    !     deallocate(errorClass)
    !     deallocate(vh1)
    !     deallocate(yh1)
    !     if (config % hiddenLayers == 2) then
    !         deallocate(vh2)
    !         deallocate(yh2)
    !     end if
    !     deallocate(vs)
    !     deallocate(y)

    ! END FUNCTION neuralNetworkActivation


    ! real(kind = 8) elemental function activation(v, afunction) result(y)
    !     implicit none

    !     real (kind = 8), intent(in) :: v
    !     integer, intent(in) :: afunction

    !     real (kind = 8), parameter :: a = 1

    !     select case(afunction)
    !     case (1) !LOGISTIC
    !         y = 1.d0/(1.d0 + exp(-a * v))
    !     case (2) !TANGENT
    !         y = (1.d0 - exp(-v))/(1.d0 + exp(-v))
    !     case (3) !GAUSS
    !         y = exp(-v)
    !     end select
    ! end function activation

    ! real(kind = 8) elemental function derivate(v, y, afunction) result(d)
    !     implicit none
    !     real (kind = 8), intent(in) :: v
    !     real (kind = 8), intent(in) :: y
    !     integer, intent(in) :: afunction
    !     real (kind = 8), parameter :: a = 1

    !     select case(afunction)
    !     case (1)!LOGISTICA: e^-x / (1+e^-x)^2
    !         d = ((a * exp(-a * v))/((1.d0 + exp(-a * v))**2.d0))
    !     case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
    !         d = (2 * exp(-v)) / ((1 + exp(-v))**2)
    !     case (3)!GAUSS
    !         d = -y/a
    !     end select
    ! end function derivate

END MODULE annActivation