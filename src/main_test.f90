PROGRAM MAIN_TEST

    USE annTest
    USE newTypes

    implicit none

    type(annConfig) :: config
    character(len = 50) :: str2, str1
    character (len =100) :: fstr, DUMMY, CMD
    character (len =30) :: fileNameBest, fileNameBestofBest
    integer :: i
    integer :: j
    integer :: k
    integer :: nExperiments
    integer :: best, iExperiment, bExperiment
    double precision :: fitness, bFitness, eFitness, eta, alpha

    integer :: hiddenLayers, neuronsLayer(2), activationFunction
    integer :: nClasses
    integer :: nClassesValidation
    integer :: nClassesGeneralization
    integer :: nClassesTest
    integer :: nInputs
    integer :: nOutputs
    integer :: nEpochs
    real :: targetError
    integer :: loadWeightsBias
    logical :: haveValidation
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

    CALL getarg(1, str2)
    
    read (str2, '(I10)')nExperiments      

    print*,'****************************************************************'
    print*, '                   RNA Test'
    print*,'****************************************************************'

    open(10, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
    read(10, content)
    close(10)

    config % nClassesTest = nClassesTest
    config % nInputs = nInputs
    config % nOutputs = nOutputs

    fstr = '(      F8.5)'
    write(fstr(2:7), '(I6)') config % nClassesTest
    open(unit=20, file='./output/final.out', STATUS='OLD',ACTION='READ')
    read(20,'(A)') DUMMY

    bFitness = 1e+10
    bExperiment = 10000
    do iExperiment = 1, nExperiments
        read(20, '(ES14.6E2,I2,I3,I3,I2,ES14.6E2,ES14.6E2)')&
                eFitness,&
                hiddenLayers, &
                neuronsLayer(1), &
                neuronsLayer(2), &
                activationFunction, &
                eta, &
                alpha
        if (eFitness < bFitness) then
            bFitness = eFitness
            bExperiment = iExperiment
            config % hiddenLayers = hiddenLayers
            config % neuronsLayer(1) = neuronsLayer(1)
            config % neuronsLayer(2) = neuronsLayer(2)
            config % activationFunction = activationFunction
            config % eta = eta
            config % alpha = alpha
        endif
    enddo

    ! print*,'Best objective function value:', bFitness
    ! print*,'Number of hidden layers:', config % hiddenLayers
    ! print*,'Neuron in hidden layer 1:', config % neuronsLayer(1)
    ! if(config % hiddenLayers < 1) then
    !     print*, 'Neurons in hidden layer 2:', config % neuronsLayer(2)
    ! end if
    ! print*,'Activation function:', config % activationFunction
    ! print*,'Alpha: ', alpha
    ! print*,'Eta: ', eta
    close(20)

    if (bExperiment < 10) then
        write(str1, '(i1)') bExperiment
    else if (bExperiment < 100) then
        write(str1, '(i2)') bExperiment
    else
        write(str1, '(i3)') bExperiment
    endif

    fileNameBest = './output/ann_e' // trim(str1) // '.best'
    fileNameBestofBest = './output/nn.best'

    
    CMD = 'cp ' //trim(fileNameBest)// ' ' // trim(fileNameBestofBest)

    ! print*, CMD
    CALL system(CMD)

    CALL neuralNetworkTest(config)
    
END PROGRAM MAIN_TEST