MODULE newTypes

IMPLICIT NONE

TYPE :: Particle
    REAL (kind = 8) :: fitness
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: solution
END TYPE Particle

TYPE :: OptionsMPCA
    INTEGER :: nDimensions
    INTEGER :: typeProbability
    INTEGER :: iCycleBlackboard
    INTEGER :: nProcessors
    INTEGER :: iProcessor
    INTEGER :: iExperiment
    INTEGER :: nParticlesProcessor
    INTEGER (kind = 8) :: iterPerturbation
    INTEGER (kind = 8) :: maxNFE
    REAL (kind = 8), DIMENSION(6) :: lowerBound
    REAL (kind = 8), DIMENSION(6) :: upperBound
    REAL (kind = 8) :: lo_small
    REAL (kind = 8) :: up_small
    REAL (kind = 8) :: fmin
    REAL (kind = 8) :: emin
    REAL (kind = 8) :: rho
    LOGICAL :: verbose

END TYPE OptionsMPCA
    
TYPE :: StatusMPCA
    INTEGER (kind = 8) :: NFE
    INTEGER (kind = 8) :: it
    INTEGER (kind = 8) :: higherNFE
    INTEGER (kind = 8) :: lastUpdate
    INTEGER (kind = 8) :: totalNFE
    INTEGER (kind = 8) :: iBest
    LOGICAL :: flag
	REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: minB
	REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: maxB
    REAL (kind = 8) :: bestObjectiveFunction
    LOGICAL :: fileUpdated
    logical :: doStop
END TYPE StatusMPCA

TYPE :: annConfig
    integer :: nClasses
    integer :: nClassesValidation
    integer :: nClassesGeneralization
    integer :: nClassesTest
    integer :: nInputs
    integer :: nOutputs
    integer :: hiddenLayers
    integer :: neuronsLayer(2)
    integer :: activationFunction
    real (kind = 8) :: alpha
    real (kind = 8) :: eta
    real (kind = 8) :: MeanSquaredErrorTrain
    real (kind = 8) :: MeanSquaredErrorValidation
    REAL (kind = 8) :: targetError
    INTEGER (kind = 8) :: nEpochs
    integer :: loadWeightsBias
    LOGICAL :: haveValidation
    LOGICAL :: tryInitialArchitecture
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: x, x_valid, x_gen, x_test
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: y, y_valid, y_gen, y_test
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: wh1
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: wh2
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:,:) :: ws
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bh1
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bh2
    REAL (kind = 8), ALLOCATABLE, DIMENSION(:) :: bs
END TYPE annConfig

END MODULE newTypes
