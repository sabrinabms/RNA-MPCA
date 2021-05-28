PROGRAM MAIN_ACTIVATION
    USE annActivation
    USE newTypes

    implicit none

    TYPE(annConfig) :: config
    REAL(kind = 8) :: MSE
    character (100) :: fString, dummy
    integer :: i
    integer :: j
    integer :: k

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

    NAMELIST /content/ nClasses, nClassesValidation,&
            nClassesGeneralization, &
            nClassesValidation, &
            nClassesActivation, &
            nInputs, nOutputs, &
            targetError, nEpochs, &
            loadWeightsBias, &
            haveValidation, &
            tryInitialArchitecture
            !tryFixedConfiguration

    print*, 'ENTROU NA ATIVACAO DA RNA'
    OPEN(1, FILE='./config/configuration.ini', STATUS='OLD', ACTION='READ')
    READ(1, content)
    CLOSE(1)

    config % nClassesActivation = nClassesActivation
    config % nInputs = nInputs
    config % nOutputs = nOutputs

    print*, 'N. Classes', config % nClassesActivation
    print*, 'N. Entradas', config % nInputs
    print*, 'N. Saidas', config % nOutputs

    allocate(config % x(config % nInputs, config % nClassesActivation))
    allocate(config % y(config % nOutputs, config % nClassesActivation))

    fString = '(      F8.5)'
    write(fString(2:7), '(I6)') config % nClassesActivation

    OPEN (2, file = './data/x_gen.txt')
    DO i = 1, config % nInputs
        READ(2, *) (config % x(i, j), j = 1, config % nClassesActivation)
    END DO
    CLOSE (2)
    
    OPEN (2, file = './data/y_gen.txt')
    DO i = 1, config % nOutputs
        READ(2, *) (config % y(i, j), j = 1, config % nClassesActivation)
    END DO
    CLOSE (2)

    !--------------------------------------------------------------------!
    !WEIGHTS AND BIASES
    !--------------------------------------------------------------------!

    open(12, file = trim('./output/nn.best'), STATUS = "old")
    read(12, *) dummy
    read(12, *) config % activationFunction
    read(12, *) config % hiddenLayers
    read(12, *) config % neuronsLayer(1)
    if (config % hiddenLayers > 1) then
        read(12, *) config % neuronsLayer(2)
    end if

    ! Allocating space for config
    allocate(config % wh1(config % nInputs, config % neuronsLayer(1)))
    allocate(config % bh1(config % neuronsLayer(1)))

    if (config % hiddenLayers > 1) then
        allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
        allocate(config % bh2(config % neuronsLayer(2)))
        allocate(config % ws(config % neuronsLayer(2), config % nOutputs))
    else
        allocate(config % wh2(1, 1))
        allocate(config % bh2(1))
        allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
    end if
    allocate(config % bs(config % nOutputs))

    fString = '(   F11.5)'
    write(fString(2:4), '(I3)') config % neuronsLayer(1)

    read(12, '(A)') dummy
    DO i = 1, config % nInputs
        read(12, *) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
    ENDDO

    read(12, *) dummy
    read(12, *) (config % bh1(k), k = 1, config % neuronsLayer(1))

    if (config % hiddenLayers == 2) then
        read(12, *) dummy
        fString = '(   F11.5)'
        write(fString(2:4), '(I3)') config % neuronsLayer(2)
        DO i = 1, config % neuronsLayer(1)
            read(12, *) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
        ENDDO

        read(12, *) dummy
        read(12, *) (config % bh2(k), k = 1, config % neuronsLayer(2))

        read(12, *) dummy
        fString = '(   F11.5)'
        write(fString(2:4), '(I3)') config % nOutputs
        DO i = 1, config % neuronsLayer(2)
            read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
        ENDDO
    else
        read(12, *) dummy
        fString = '(   F11.5)'
        write(fString(2:4), '(I3)') config % nOutputs
        allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
        DO i = 1, config % neuronsLayer(1)
            read(12, *) (config % ws(i, k), k = 1, config % nOutputs)
        ENDDO
    end if

    read(12, *) dummy
    read(12, *) (config % bs(k), k = 1, config % nOutputs)
    close(12)

    config % activationFunction = 2
    MSE = neuralNetworkActivation(config)

    deallocate(config % bh1)
    deallocate(config % bs)
    deallocate(config % wh1)
    deallocate(config % ws)
    if (config % hiddenLayers == 2) then
        deallocate(config % bh2)
        deallocate(config % wh2)
    end if
    
END PROGRAM MAIN_ACTIVATION