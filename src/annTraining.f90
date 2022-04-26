!***************************************************************!
! Optimization of ANN architecture by metaheuristic MPCA        !
!***************************************************************!
! Developed by: Juliana Anochi and Sabrina Sambatti (CAP/INPE)  !
! Modified by: Reynier Hernandez Torres (CAP/INPE)              !
!***************************************************************!
MODULE annTraining
    USE newTypes
    USE foul
    USE rnaFunctions
    USE annGeneralization !alteracao para a funcao objetivo do Haroldo


CONTAINS

    REAL(kind = 8) FUNCTION neuralNetworkTraining(solution, op, st, config)

        IMPLICIT NONE

        real (8), intent(in) :: solution(:)
        type(OptionsMPCA), intent(in) :: op
        type(StatusMPCA), intent(inout) :: st
        type(annConfig) :: config

        real (8) :: rNumber
        real (8) :: dOutput
        real (8) :: penaltyObj
        real (8) :: aux
        real (8), allocatable, dimension(:) :: vs
        real (8), allocatable, dimension(:) :: vh1
        real (8), allocatable, dimension(:) :: vh2
        real (8), allocatable, dimension(:) :: yh1
        real (8), allocatable, dimension(:) :: yh2
        real (8), allocatable, dimension(:) :: ys
        real (8), allocatable, dimension(:) :: gradientOutput
        real (8), allocatable, dimension(:) :: gradientHiddenLayer1
        real (8), allocatable, dimension(:) :: gradientHiddenLayer2

        real (8), allocatable, dimension(:) :: errorClassTrain
        real (8), allocatable, dimension(:) :: errorClassValid
        real (8), allocatable, dimension(:,:) :: errorTrain
        real (8), allocatable, dimension(:,:) :: errorValid

        real (8), allocatable, dimension(:,:) :: deltaWeightOutput
        real (8), allocatable, dimension(:,:) :: deltaWeightHiddenLayer1
        real (8), allocatable, dimension(:,:) :: deltaWeightHiddenLayer2

        real (8), allocatable, dimension(:) :: deltaBiasOutput
        real (8), allocatable, dimension(:) :: deltaBiasHiddenLayer1
        real (8), allocatable, dimension(:) :: deltaBiasHiddenLayer2

        real (8), allocatable, dimension(:,:) :: deltaWeightOutputLast
        real (8), allocatable, dimension(:,:) :: deltaWeightHiddenLayer1Last
        real (8), allocatable, dimension(:,:) :: deltaWeightHiddenLayer2Last

        real (8), allocatable, dimension(:) :: deltaBiasOutputLast
        real (8), allocatable, dimension(:) :: deltaBiasHiddenLayer1Last
        real (8), allocatable, dimension(:) :: deltaBiasHiddenLayer2Last
        
        real (8), parameter :: alphaObj = 0.1D0
        real (8), parameter :: betaObj = 1.0D0
        real (8), parameter :: p1 = 5.0e-8
        real (8), parameter :: p2 = 5.0e-5

        character (100) :: str0
        character (100) :: str1
        character (100) :: fString
        character (100) :: dummy

        integer :: i
        integer :: j
        integer :: jj
        integer :: k
        integer :: epoch

        config % hiddenLayers = ceiling(solution(1))
        config % neuronsLayer(1) = ceiling(solution(2))
        config % neuronsLayer(2) = ceiling(solution(3))
        config % activationFunction = ceiling(solution(4))
        config % alpha = solution(5)
        config % eta = solution(6)

        ! Allocating space for config
        !PARA PRIMEIRA CAMADA: ALOCANDO VARIAVEL
        allocate(config % wh1(config % nInputs, config % neuronsLayer(1)))
        allocate(config % bh1(config % neuronsLayer(1)))
        allocate(vh1(config % neuronsLayer(1)))
        allocate(yh1(config % neuronsLayer(1)))

        !SE EXISTIR SEGUNDA CAMADA: ALOCANDO VARIAVEL
        if (config % hiddenLayers == 2) then
            allocate(config % wh2(config % neuronsLayer(1), config % neuronsLayer(2)))
            allocate(config % bh2(config % neuronsLayer(2)))
            allocate(vh2(config % neuronsLayer(2)))
            allocate(yh2(config % neuronsLayer(2)))
            allocate(config % ws(config % neuronsLayer(2), config % nOutputs))
        else
            allocate(config % ws(config % neuronsLayer(1), config % nOutputs))
        endif

        !PARA CAMADA DE SAIDA: ALOCANDO VARIAVEL
        allocate(config % bs(config % nOutputs))
        allocate(vs(config % nOutputs))
        allocate(ys(config % nOutputs))

        !Allocating space for error variables
        allocate(errorTrain(config % nOutputs, config % nClasses))
        allocate(errorClassTrain(config % nClasses))
	    if (config % haveValidation .eqv. .true.) then
		    allocate(errorValid(config % nOutputs, config % nClassesValidation))
	        allocate(errorClassValid(config % nClassesValidation))
	    endif
        
        !Allocating space for delta variables
        !PARA PRIMEIRA CAMADA: ALOCANDO VARIAVEL
        allocate(deltaBiasHiddenLayer1(config % neuronsLayer(1)))
        allocate(deltaBiasHiddenLayer1Last(config % neuronsLayer(1)))
        allocate(deltaWeightHiddenLayer1(config % nInputs, config % neuronsLayer(1)))
        allocate(deltaWeightHiddenLayer1Last(config % nInputs, config % neuronsLayer(1)))

        !SE EXISTIR SEGUNDA CAMADA: ALOCANDO VARIAVEL
        if (config % hiddenLayers == 2) then
            allocate(deltaBiasHiddenLayer2(config % neuronsLayer(2)))
            allocate(deltaBiasHiddenLayer2Last(config % neuronsLayer(2)))
            allocate(deltaWeightHiddenLayer2(config % neuronsLayer(1), config % neuronsLayer(2)))
            allocate(deltaWeightHiddenLayer2Last(config % neuronsLayer(1), config % neuronsLayer(2)))
        endif

        allocate(deltaBiasOutput(config % nOutputs))
        allocate(deltaBiasOutputLast(config % nOutputs))
        allocate(deltaWeightOutput(config % neuronsLayer(config % hiddenLayers), config % nOutputs))
        allocate(deltaWeightOutputLast(config % neuronsLayer(config % hiddenLayers), config % nOutputs))

        deltaWeightOutput = 0
        deltaWeightHiddenLayer1 = 0

        if (config % hiddenLayers == 2) then
	        deltaWeightHiddenLayer2 = 0
	        deltaBiasHiddenLayer2 = 0
	    endif

        deltaBiasOutput = 0
        deltaBiasHiddenLayer1 = 0

        !Allocating space for gradient variables
        allocate(gradientHiddenLayer1(config % neuronsLayer(1)))

        if (config % hiddenLayers == 2) then
            allocate(gradientHiddenLayer2(config % neuronsLayer(2)))
        endif

        allocate(gradientOutput(config % nOutputs))

        gradientOutput = 0
        gradientHiddenLayer1 = 0

	    if (config % hiddenLayers == 2) then
        	gradientHiddenLayer2 = 0
	    endif

        neuralNetworkTraining = 0
        
        !--------------------------------------------------------------------!
        !INITIAL WEIGHTS AND BIASES
        !--------------------------------------------------------------------!

        SELECT CASE (config % loadWeightsBias)
        CASE (0) ! All parameters initialized with 0.5
            do i = 1, config % nInputs
                do k = 1, config % neuronsLayer(1)
                    config % wh1(i, k) = 0.5
                enddo
            enddo

            do k = 1, config % neuronsLayer(1)
                config % bh1(k) = 0.5
            enddo

            if (config % hiddenLayers == 2) then
                do i = 1, config % neuronsLayer(1)
                    do k = 1, config % neuronsLayer(2)
                        config % wh2(i, k) = 0.5
                    enddo
                enddo
                do k = 1, config % neuronsLayer(2)
                    config % bh2(k) = 0.5
                enddo
            endif

            do i = 1, config % neuronsLayer(config % hiddenLayers)
                do k = 1, config % nOutputs
                    config % ws(i, k) = 0.5
                enddo
            enddo
            do k = 1, config % nOutputs
                config % bs(k) = 0.5
            enddo

        CASE (1) ! All parameters randomly initialized
            do i = 1, config % nInputs
                do k = 1, config % neuronsLayer(1)
                    call random_number(rNumber)
                    config % wh1(i, k) = rNumber
                enddo
            enddo
            do k = 1, config % neuronsLayer(1)
                call random_number(rNumber)
                config % bh1(k) = rNumber
            enddo

            if (config % hiddenLayers == 2) then
                do i = 1, config % neuronsLayer(1)
                    do k = 1, config % neuronsLayer(2)
                        call random_number(rNumber)
                        config % wh2(i, k) = rNumber
                    enddo
                enddo
                do k = 1, config % neuronsLayer(2)
                    call random_number(rNumber)
                    config % bh2(k) = rNumber
                enddo
            endif

            do i = 1, config % neuronsLayer(config % hiddenLayers)
                do k = 1, config % nOutputs
                    call random_number(rNumber)
                    config % ws(i, k) = rNumber
                enddo
            enddo
            do k = 1, config % nOutputs
                call random_number(rNumber)
                config % bs(k) = rNumber
            enddo

        CASE (2) !Empirical parameters
            open(12, file = trim('./datain/nn.initial'), STATUS = "old")
            read(12, *) config % hiddenLayers
            read(12, *) config % neuronsLayer(1)
            if(config % hiddenLayers == 2) then
                read(12, *) config % neuronsLayer(2)
            endif

            read(12, '(A)') dummy
            fString = '(   F11.5)'
            write(fString(2:4), '(I3)') config % neuronsLayer(1)
            do i = 1, config % nInputs
                read(12, fString) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
            enddo

            read(12, '(A)') dummy
            read(12, fString) (config % bh1(k), k = 1, config % neuronsLayer(1))

            if (config % hiddenLayers == 2) then
                read(12, '(A)') dummy
                write(fString(2:4), '(I3)') config % neuronsLayer(2)
                do i = 1, config % neuronsLayer(1)
                    read(12, fString) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
                enddo

                read(12, '(A)') dummy
                read(12, fString) (config % bh2(k), k = 1, config % neuronsLayer(2))
                
                write(fString(2:4), '(I3)') config % nOutputs
                read(12, '(A)') dummy
                do i = 1, config % neuronsLayer(2)
                    read(12, fString) (config % ws(i, k), k = 1, config % nOutputs)
                enddo
            else
                write(fString(2:4), '(I3)') config % nOutputs
                read(12, '(A)') dummy
                do i = 1, config % neuronsLayer(1)
                    read(12, fString) (config % ws(i, k), k = 1, config % nOutputs)
                enddo            
            endif

            read(12, '(A)') dummy
            read(12, fString) (config % bs(k), k = 1, config % nOutputs)

            close(12)
        END SELECT

        !----------------------------------------------------------------------!
        ! FEEDFORWARD: BEGIN OF THE ANN
        !----------------------------------------------------------------------!
        epoch = 0
        
        if (op % iProcessor == 0) then
            ! achar(): returns the character located at position I in the ASCII collating sequence
            write(*,FMT="(A1,A,t25,I10,A,I10)",ADVANCE="NO") achar(13), &
                & "NFE (in processor 0): ", &
                & st % NFE, &
                & " of ", &
                & op % maxNFE / op % nProcessors
        endif

        do while (epoch .LT. config % nEpochs)
        
            deltaWeightOutputLast = deltaWeightOutput
            deltaWeightHiddenLayer1Last = deltaWeightHiddenLayer1
	        
            if (config % hiddenLayers == 2) then
            	deltaWeightHiddenLayer2Last = deltaWeightHiddenLayer2
                deltaBiasHiddenLayer2Last = deltaBiasHiddenLayer2
            endif
	        
            deltaBiasOutputLast = deltaBiasOutput
            deltaBiasHiddenLayer1Last = deltaBiasHiddenLayer1
            epoch = epoch + 1

            do i = 1, config % nClasses
                ! ACTIVATING HIDDEN LAYER 1
                vh1 = matmul(config % x(:, i), config % wh1) - config % bh1
                yh1 = activation(vh1, config % activationFunction)

                if (config % hiddenLayers == 1) then
                    vs = matmul(yh1, config % ws) - config % bs
		        else
                    vh2 = matmul(yh1, config % wh2) - config % bh2
                    yh2 = activation(vh2, config % activationFunction)
                    vs = matmul(yh2, config % ws) - config % bs
                endif

                ! ACTIVATING OUTPUT
                ys = activation(vs, config % activationFunction)
                errorTrain(:, i) = config % y(:, i) - ys

                !CALCULO PADRAO do ERRO
                errorClassTrain(i) = sum(errorTrain(:, i), dim = 1)
                errorClassTrain(i) = 0.5d0 * (errorClassTrain(i)**2.d0)

                !-------------------------------------------------------------------------!
                !                        BACKPROPAGATION
                !-------------------------------------------------------------------------!
                !TRAINING OUTPUT LAYER
                do j = 1, config % nOutputs
                    dOutput = derivate(vs(j), ys(j), config % activationFunction)
                    gradientOutput(j) = errorTrain(j, i) * dOutput

                    if (config % hiddenLayers == 1) then
                        deltaWeightOutput(:, j) = config % eta * gradientOutput(j) * yh1
                    else
                        deltaWeightOutput(:, j) = config % eta * gradientOutput(j) * yh2
                    endif

                    deltaBiasOutput(j) = dfloat(-1) * config % eta * gradientOutput(j)

                    config % ws(:, j) = config % ws(:, j) &
                        + deltaWeightOutput(:, j) &
                        + config % alpha * deltaWeightOutputLast(:, j)

                    config % bs(j) = config % bs(j) &
                        + deltaBiasOutput(j) &
                        + config % alpha * deltaBiasOutputLast(j)
                enddo

                !TRAINING HIDDEN LAYER 2
                if (config % hiddenLayers == 2) then
                    do j = 1, config % neuronsLayer(2)
                        dOutput = derivate(vh2(j), yh2(j), config % activationFunction)
			            aux = 0.d0               

			            do k = 1, config % nOutputs
                            aux = aux + (gradientOutput(k) * config % ws(j, k))
                        enddo

                        gradientHiddenLayer2(j) = aux * dOutput
                        deltaWeightHiddenLayer2(:, j) = config % eta * gradientHiddenLayer2(j) * yh1
                        deltaBiasHiddenLayer2(j) = dfloat(-1) * config % eta * gradientHiddenLayer2(j)

                        config % wh2(:, j) = config % wh2(:, j) &
                                & + deltaWeightHiddenLayer2(:, j) &
                                & + config % alpha * deltaWeightHiddenLayer2Last(:, j)

                        config % bh2(j) = config % bh2(j) &
                                & + deltaBiasHiddenLayer2(j) &
                                & + config % alpha * deltaBiasHiddenLayer2Last(j)
                    end do
                endif

                ! TRAINING HIDDEN LAYER 1
                do j = 1, config % neuronsLayer(1)
                    dOutput = derivate(vh1(j), yh1(j), config % activationFunction)
	                aux = 0.d0

                    if (config % hiddenLayers == 1) then
                        do k = 1, config % nOutputs
                            aux = aux + (gradientOutput(k) * config % ws(j, k))
                        enddo
                    else
                        do jj = 1, config % neuronsLayer(2)
                            do k = 1, config % nOutputs
                                aux = aux + (gradientHiddenLayer2(jj) * config % wh2(jj, k))
                            enddo
                        enddo
                    endif

                    gradientHiddenLayer1(j) =  aux * dOutput
                    deltaWeightHiddenLayer1(:, j) = config % eta * gradientHiddenLayer1(j) * config % x(:, i)
                    deltaBiasHiddenLayer1(j) = dfloat(-1) * config % eta * gradientHiddenLayer1(j)

                    config % wh1(:, j) = config % wh1(:, j) &
                        & + deltaWeightHiddenLayer1(:, j) &
                        & + config % alpha * deltaWeightHiddenLayer1Last(:, j)

                    config % bh1(j) = config % bh1(j) &
                        & + deltaBiasHiddenLayer1(j) &
                        & + config % alpha * deltaBiasHiddenLayer1Last(j)
                enddo
            enddo

            config % MeanSquaredErrorTrain = sum(errorClassTrain) / dfloat(config % nClasses)

            !***********************************************************************
            ! CROSS VALIDATION
            !***********************************************************************
            if (config % haveValidation .eqv. .true.) then
                do i = 1, config % nClassesValidation
                    ! ACTIVATING HIDDEN LAYER 1
                    vh1 = matmul(config % x_valid(:, i), config % wh1) - config % bh1
                    yh1 = activation(vh1, config % activationFunction)

                    if (config % hiddenLayers == 1) then
                        vs = matmul(yh1, config % ws) - config % bs
                    else
                        vh2 = matmul(yh1, config % wh2) - config % bh2
                        yh2 = activation(vh2, config % activationFunction)
                        vs = matmul(yh2, config % ws) - config % bs
                    endif

                    ! ACTIVATING OUTPUT
                    ys = activation(vs, config % activationFunction)

                    errorValid(:, i) = config % y_valid(:, i) - ys

                    errorClassValid(i) = sum(errorValid(:, i), dim = 1)
                    errorClassValid(i) = 0.5d0 * (errorClassValid(i)**2.d0)
                end do

                config % MeanSquaredErrorValidation = sum(errorClassValid) / dfloat(config % nClassesValidation)
            endif
        enddo

        ! Objective Function (Carvalho,2011)
        penaltyObj = 1
            !& + 10000 * p1 * exp(dfloat(config % neuronsLayer(1))) &
            !& + 1e+6 * p1 * exp(dfloat(config % neuronsLayer(2))) &
            !& + p2 * dfloat(config % nEpochs) &

! 	 ABAIXO EXISTEM 3 FUNCOES OBJETIVOS QUE DEVEM SER COMPARADAS
!        PRIMEIRA: USANDO ERRO DA VALIDACAO CRUZADA
!        SEGUNDA: USANDO APENAS O ERRO QUADRATICO DO TREINAMENTO
!        TERCEIRA: USANDO O ERRO DE GENERALIZACAO

!        if (config % haveValidation .eqv. .true.) THEN   

! PRIMEIRA FUNCAO OBJETIVO
!            neuralNetworkTraining = penaltyObj &
!                & * ((alphaObj * config % MeanSquaredErrorValidation + betaObj * config % MeanSquaredErrorValidation) &
!                & / (alphaObj + betaObj))
!        else

! SEGUNDA FUNCAO OBJETIVO
!            neuralNetworkTraining = penaltyObj * config % MeanSquaredError
!        ENDif

! TERCEIRA FUNCAO OBJETIVO

        config % meanSquaredErrorGen = neuralNetwork(config) ! funcao objetivo do Haroldo, aqui é chamado: annGeneralization.f90
        neuralNetworkTraining = penaltyObj &
                & * ((alphaObj * config % MeanSquaredErrorTrain + betaObj * config % meanSquaredErrorGen) &
                & / (alphaObj + betaObj))

        ! Store configuration if objFunction is best
        if (neuralNetworkTraining < st % bestObjectiveFunction) then
            
            st % bestObjectiveFunction = neuralNetworkTraining

            if (op % iProcessor < 9) THEN
                WRITE (str0, '(I1)') op % iProcessor + 1
            else if (op % iProcessor < 99) THEN
                WRITE (str0, '(I2)') op % iProcessor + 1
            else
                WRITE (str0, '(I3)') op % iProcessor + 1
            endif

            if (op % iExperiment < 10) THEN
                WRITE (str1, '(I1)') op % iExperiment
            else if (op % iExperiment < 100) THEN
                WRITE (str1, '(I2)') op % iExperiment
            else
                WRITE (str1, '(I3)') op % iExperiment
            endif

            OPEN(12, FILE = './dataout/ann_e' // trim(str1) // '_p' // trim(str0) // '.out')

            write(12, '(A)') 'Objective Function (MPCA): '
            write(12, '(ES14.6E2)') st % bestObjectiveFunction
            write(12, '(A)') 'Activation Function: '
     	    write(12, '(I3)') config % activationFunction
     	    write(12, '(A)') 'Hidden Layers: '
            write(12, '(I3)') config % hiddenLayers
            write(12, '(A)') 'Neurons First Layer: '
            write(12, '(I3)') config % neuronsLayer(1)

            if (config % hiddenLayers == 2) then
                write(12, '(A)') 'Neurons Second Layer: '
                write(12, '(I3)') config % neuronsLayer(2)
            endif

            write(12, '(A)') 'wh1'
            fString = '(   F11.5)'
            write(fString(2:4), '(I3)') config % neuronsLayer(1)

            do i = 1, config % nInputs
                write(12, fString) (config % wh1(i, k), k = 1, config % neuronsLayer(1))
            enddo

            write(12, '(A)') 'bh1'
            write(12, fString) (config % bh1(k), k = 1, config % neuronsLayer(1))

            if (config % hiddenLayers == 2) then
                write(12, '(A)') 'wh2'
                write(fString(2:4), '(I3)') config % neuronsLayer(2)
                do i = 1, config % neuronsLayer(1)
                    write(12, fString) (config % wh2(i, k), k = 1, config % neuronsLayer(2))
                enddo

                write(12, '(A)') 'bh2'
                write(12, fString) (config % bh2(k), k = 1, config % neuronsLayer(2))

                write(fString(2:4), '(I3)') config % nOutputs
                write(12, '(A)') 'ws'
                do i = 1, config % neuronsLayer(2)
                    write(12, fString) (config % ws(i, k), k = 1, config % nOutputs)
                enddo
            else
                write(fString(2:4), '(I3)') config % nOutputs
                write(12, '(A)') 'ws'
                do i = 1, config % neuronsLayer(1)
                    write(12, fString) (config % ws(i, k), k = 1, config % nOutputs)
                enddo
            endif

            write(12, '(A)') 'bs'
            write(12, fString) (config % bs(k), k = 1, config % nOutputs)

            write(12, '(A)', advance = 'no') 'Alpha: '
            write(12, '(F13.4)') config % alpha
            write(12, '(A)', advance = 'no') 'Eta: '
            write(12, '(F13.4)') config % eta
            write(12, '(A)', advance = 'no') 'Activation Function: '

            select case(config % activationFunction)
            case (1)
                write(12, '(A)') 'LOGISTIC'
            case (2)
                write(12, '(A)') 'TANGENT'
            case (3)
                write(12, '(A)') 'GAUSS'
            end select

            write(12, '(A)', advance = 'no') 'MeanSquaredError: '
            write(12, '(ES14.6E2)') config % MeanSquaredErrorTrain

            if (config % haveValidation .eqv. .true.) then
                WRITE(12, '(A)', advance = 'no') 'MeanSquaredError - validation: '
                WRITE(12, '(ES14.6E2)') config % MeanSquaredErrorValidation
            endif

            WRITE(12, '(A)', advance = 'no') 'Epoch :'
            WRITE(12, '(I10)') epoch

            WRITE(12, '(A)', advance = 'no') 'NFE: '
            WRITE(12, '(I10)') st % NFE

            CLOSE(12)
        endif

        neuralNetworkTraining = st % bestObjectiveFunction

        deallocate(errorTrain)
        deallocate(errorClassTrain)

        if (config % haveValidation .eqv. .true.) then
                deallocate(errorValid)
                deallocate(errorClassValid)
        endif

        deallocate(vh1)
        deallocate(yh1)

        if (config % hiddenLayers == 2) then
            deallocate(vh2)
            deallocate(yh2)
        endif

        deallocate(vs)
        deallocate(ys)

        deallocate(config % bh1)
        deallocate(config % bs)
        deallocate(config % wh1)
        deallocate(config % ws)

        if (config % hiddenLayers == 2) then
            deallocate(config % bh2)
            deallocate(config % wh2)
        endif

        deallocate(deltaBiasHiddenLayer1)
        deallocate(deltaBiasHiddenLayer1Last)

        if (config % hiddenLayers == 2) then
	        deallocate(deltaBiasHiddenLayer2)
            deallocate(deltaBiasHiddenLayer2Last)
	        deallocate(deltaWeightHiddenLayer2)
            deallocate(deltaWeightHiddenLayer2Last)
            deallocate(gradientHiddenLayer2)
	    endif

	    deallocate(deltaBiasOutput)
        deallocate(deltaBiasOutputLast)
        deallocate(deltaWeightHiddenLayer1)
        deallocate(deltaWeightHiddenLayer1Last)
        deallocate(deltaWeightOutput)
        deallocate(deltaWeightOutputLast)
        deallocate(gradientHiddenLayer1)
        deallocate(gradientOutput)

    END FUNCTION neuralNetworkTraining

 

END MODULE annTraining
