PROGRAM MAIN_ACTIVATION

    USE annActivation
    USE newTypes

    implicit none

    TYPE(annConfig) :: config
    REAL(kind = 8) :: MSE
    CHARACTER(len = 50) :: string, str1
    character (100) :: fString, dummy, CMD
    character (30) :: fileNameBest, fileNameBestofBest
    integer :: i
    integer :: j
    integer :: k
    integer :: nExperiments
    integer :: best, iExperiment, bExperiment
    double precision :: fitness, bfitness, efitness, eta, alpha

    integer :: hiddenLayers, neuronsLayer(2), activationFunction
    INTEGER (kind = 8) :: nClasses
    INTEGER (kind = 8) :: nClassesValidation
    INTEGER (kind = 8) :: nClassesGeneralization
    INTEGER (kind = 8) :: nClassesTest
    INTEGER (kind = 8) :: nInputs
    INTEGER (kind = 8) :: nOutputs
    INTEGER (kind = 8) :: nEpochs
    REAL (kind = 8) :: targetError
    integer :: loadWeightsBias
    LOGICAL :: haveValidation
    logical :: tryInitialArchitecture

    NAMELIST /content/ nClasses, nClassesValidation,&
            nClassesGeneralization, &
            nClassesValidation, &
            nClassesTest, &
            nInputs, nOutputs, &
            targetError, nEpochs, &
            loadWeightsBias, &
            haveValidation, &
            tryInitialArchitecture

    !*********************************************
    ! SETTING PARAMETERS
    ! 1 - Number of the experiment
    CALL getarg(1, string)
    READ (string, '(I10)')nExperiments      
    print*,'verificando -->', nExperiments

    print*,'****************************************************************'
    print*, '                   RNA Activation'
    print*,'****************************************************************'
    OPEN(10, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
    READ(10, content)
    CLOSE(10)

    config % nClassesTest = nClassesTest
    config % nInputs = nInputs
    config % nOutputs = nOutputs

    allocate(config % x(config % nInputs, config % nClassesTest))
    allocate(config % y(config % nOutputs, config % nClassesTest))

    fString = '(      F8.5)'
    write(fString(2:7), '(I6)') config % nClassesTest
    print*,' Number of activation classes: ', config % nClassesTest

    open(unit=20, file='./output/final.out', STATUS='OLD',ACTION='READ')
    READ(20,'(A)') dummy

    bfitness = 1e+10
    do iExperiment = 1, nExperiments
        read(20, '(ES14.6E2,I2,I3,I3,I2,ES14.6E2,ES14.6E2)')&
                efitness,&
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
        endif
    enddo

    print*,'Best objective function value:', bfitness
    print*,'Number of hidden layers:', config % hiddenLayers
    print*,'Neuron in hidden layer 1:', config % neuronsLayer(1)
    if(config % hiddenLayers < 1) then
        print*, 'Neurons in hidden layer 2:', config % neuronsLayer(2)
    end if
    print*,'Activation function:', config % activationFunction
    print*,'Alpha: ', alpha
    print*,'Eta: ', eta
    close(20)

    if (bExperiment < 10) then
        write(str1, '(i1)') bExperiment
    else if (bExperiment < 100) then
        write(str1, '(i2)') bExperiment
    else
        write(str1, '(i3)') bExperiment
    endif
    ! ann_e2.best

    fileNameBest = './output/ann_e' // trim(str1) // '.best'
    fileNameBestofBest = './output/nn.best'

    print*,'Best result:', fileNameBest
    
    CMD = 'cp ' //trim(fileNameBest)// ' ' // trim(fileNameBestofBest)

    print*, CMD
    CALL system(CMD)

    call annActivation(fileNameBestofBest, config)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ! OPEN (2, file = './data/x_test.txt')
    ! DO i = 1, config % nInputs
    !     READ(2, *) (x_test(i, j), j = 1, config % nClassesTest)
    ! END DO
    ! CLOSE (2)
    
    ! OPEN (2, file = './data/y_test.txt')
    ! DO i = 1, config % nOutputs
    !     READ(2, *) (y_test(i, j), j = 1, config % nClassesTest)
    ! END DO
    ! CLOSE (2)

    ! !--------------------------------------------------------------------!
    ! !WEIGHTS AND BIASES
    ! !--------------------------------------------------------------------!

    ! open(12, fileNameBestofBest, STATUS = "old")
    ! read(12, *) dummy !Valor funcao objetivo
    ! read(12, *) config % activationFunction
    ! read(12, *) config % hiddenLayers
    ! read(12, *) config % neuronsLayer(1)
    ! print*,' Number of hidden layers: ', config % hiddenLayers
    ! print*,' Neurons in hidden layer 1: ',  config % neuronsLayer(1)
    
    ! if (config % hiddenLayers > 1) then
    !     read(12, *) config % neuronsLayer(2)
    !     print*,' Neurons in hidden layer 2: ', config % neuronsLayer(2)
    ! endif
    
    ! ! Allocating space for config
    ! allocate(config % wh1(config % nInputs, config % neuronsLayer(1)))
    ! allocate(config % bh1(config % neuronsLayer(1)))

    ! if (config % hiddenLayers > 1) then
    !     allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
    !     allocate(config % bh2(config % neuronsLayer(2)))
    !     allocate(config % ws(config % neuronsLayer(2), config % nOutputs))
    ! else
    !     allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
    ! end if

    ! allocate(config % bs(config % nOutputs))

    ! fString = '(   F11.5)'
    ! write(fString(2:4), '(I3)') config % neuronsLayer(1)

    ! read(12, '(A)') dummy
    ! DO i = 1, config % nInputs
    !     read(12, *) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
    ! ENDDO

    ! read(12, *) dummy
    ! read(12, *) (config % bh1(k), k = 1, config % neuronsLayer(1))

    ! if (config % hiddenLayers == 2) then
    !     read(12, *) dummy
    !     fString = '(   F11.5)'
    !     write(fString(2:4), '(I3)') config % neuronsLayer(2)
    !     DO i = 1, config % neuronsLayer(1)
    !         read(12, *) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
    !     ENDDO

    !     read(12, *) dummy
    !     read(12, *) (config % bh2(k), k = 1, config % neuronsLayer(2))

    !     read(12, *) dummy
    !     fString = '(   F11.5)'
    !     write(fString(2:4), '(I3)') config % nOutputs
    !     DO i = 1, config % neuronsLayer(2)
    !         read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
    !     ENDDO
    ! else
    !     read(12, *) dummy
    !     fString = '(   F11.5)'
    !     write(fString(2:4), '(I3)') config % nOutputs
    !     allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
    !     DO i = 1, config % neuronsLayer(1)
    !         read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
    !     ENDDO
    ! end if

    ! read(12, *) dummy
    ! read(12, *) (config % bs(k), k = 1, config % nOutputs)
    ! close(12)

    ! MSE = neuralNetworkActivation(config)
    ! print*,'Activation - Mean Square Error: ', MSE

    ! deallocate(config % bh1)
    ! deallocate(config % bs)
    ! deallocate(config % wh1)
    ! deallocate(config % ws)
    ! if (config % hiddenLayers == 2) then
    !     deallocate(config % bh2)
    !     deallocate(config % wh2)
    ! end if
    
END PROGRAM MAIN_ACTIVATION