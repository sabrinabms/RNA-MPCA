MODULE rnaFunctions

    USE newTypes

CONTAINS

   real(kind = 8) elemental function activation(v, afunction) result(y)
        implicit none

        real (kind = 8), intent(in) :: v
        integer, intent(in) :: afunction

        real (kind = 8), parameter :: a = 1

        select case(afunction)
        case (1) !LOGISTIC
            y = 1.d0/(1.d0 + exp(-a * v))
        case (2) !TANGENT
            y = (1.d0 - exp(-v))/(1.d0 + exp(-v))
        case (3) !GAUSS
            y = exp(-v)
        end select
    end function activation

    real(kind = 8) elemental function derivate(v, y, afunction) result(d)
        implicit none
        real (kind = 8), intent(in) :: v
        real (kind = 8), intent(in) :: y
        integer, intent(in) :: afunction
        real (kind = 8), parameter :: a = 1

        select case(afunction)
        case (1)!LOGISTICA: e^-x / (1+e^-x)^2
            d = ((a * exp(-a * v))/((1.d0 + exp(-a * v))**2.d0))
        case (2)!TANGENTE: 2e^-x / (1+e^-x)^2
            d = (2 * exp(-v)) / ((1 + exp(-v))**2)
        case (3)!GAUSS
            d = -y/a
        end select
    end function derivate

END MODULE rnaFunctions
