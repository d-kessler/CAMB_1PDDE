    module DarkEnergyInterface
    use precision
    use interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkEnergyModel
        logical :: is_cosmological_constant = .true.
        integer :: num_perturb_equations = 0
    contains
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: diff_rhopi_Add_Term
    procedure :: PerturbationInitial
    procedure :: PerturbationEvolve
    procedure :: PrintFeedback
    ! do not have to implement w_de or grho_de if BackgroundDensityAndPressure is inherited directly
    procedure :: w_de
    procedure :: grho_de
    procedure :: Effective_w_wa !Used as approximate values for non-linear corrections
    end type TDarkEnergyModel

    type, extends(TDarkEnergyModel) :: TDarkEnergyEqnOfState
        !Type supporting w, wa or general w(z) table
        real(dl) :: w_lam = -1_dl !p/rho for the dark energy (an effective value, used e.g. for halofit)
        real(dl) :: wa = 0._dl !may not be used, just for compatibility with e.g. halofit
        real(dl) :: cs2_lam = 1_dl !rest-frame sound speed, though may not be used
        logical :: use_tabulated_w = .false.  !Use interpolated table; note this is quite slow.
        logical :: no_perturbations = .false. !Don't change this, no perturbations is unphysical
        integer :: dde_model = 0 !dynamical dark energy model to use
        integer :: pm_in = 1 !determines the behavior of dark energy equation of state
        real(dl) :: dde_freq = 1_dl !frequency of oscillations for some dynamical models
        !Interpolations if use_tabulated_w=.true.
        Type(TCubicSpline) :: equation_of_state, logdensity
    contains
    procedure :: ReadParams => TDarkEnergyEqnOfState_ReadParams
    procedure :: Init => TDarkEnergyEqnOfState_Init
    procedure :: SetwTable => TDarkEnergyEqnOfState_SetwTable
    procedure :: PrintFeedback => TDarkEnergyEqnOfState_PrintFeedback
    procedure :: w_de => TDarkEnergyEqnOfState_w_de
    procedure :: grho_de => TDarkEnergyEqnOfState_grho_de
    procedure :: Effective_w_wa => TDarkEnergyEqnOfState_Effective_w_wa
    end type TDarkEnergyEqnOfState

    public TDarkEnergyModel, TDarkEnergyEqnOfState
    contains

    function w_de(this, a)
    class(TDarkEnergyModel) :: this
    real(dl) :: w_de, al
    real(dl), intent(IN) :: a

    w_de = -1._dl

    end function w_de  ! equation of state of the PPF DE

    function grho_de(this, a)  !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyModel) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    grho_de =0._dl

    end function grho_de

    subroutine PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyModel) :: this
    integer, intent(in) :: FeedbackLevel

    end subroutine PrintFeedback


    subroutine Init(this, State)
    use classes
    class(TDarkEnergyModel), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    end subroutine Init

    subroutine BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w

    ! if (this%is_cosmological_constant) then
    !    grhov_t = grhov * a * a
    !    if (present(w)) w = -1_dl
    ! else
    !    ! Ensure a valid result
    !    if (a > 1e-10) then
    !        grhov_t = grhov * this%grho_de(a) / (a * a)
    !    else
    !        grhov_t = 0._dl
    !    end if
    !    if (present(w)) w = this%w_de(a)
    ! end if

    ! Ensure a valid result
    if (a > 1e-10) then
        grhov_t = grhov * this%grho_de(a) / (a * a)
    else
        grhov_t = 0._dl
    end if
    if (present(w)) w = this%w_de(a)

    end subroutine BackgroundDensityAndPressure

    subroutine Effective_w_wa(this, w, wa)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = -1
    wa = 0

    end subroutine Effective_w_wa


    subroutine PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyModel), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe=0
    dgqe=0

    end subroutine PerturbedStressEnergy


    function diff_rhopi_Add_Term(this, dgrhoe, dgqe,grho, gpres, w, grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(in) :: dgrhoe, dgqe, grho, gpres, grhok, w, adotoa, &
        k, grhov_t, z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot

    ! Ensure, that the result is set, when the function is not implemented by
    ! subclasses
    ppiedot = 0._dl

    end function diff_rhopi_Add_Term

    subroutine PerturbationEvolve(this, ayprime, w, w_ix, a, adotoa, k, z, y)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, y(:), w
    integer, intent(in) :: w_ix
    end subroutine PerturbationEvolve

    subroutine PerturbationInitial(this, y, a, tau, k)
    class(TDarkEnergyModel), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    !Get intinitial values for perturbations at a (or tau)
    !For standard adiabatic perturbations can usually just set to zero to good accuracy

    y = 0

    end subroutine PerturbationInitial


    subroutine TDarkEnergyEqnOfState_SetwTable(this, a, w, n)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), w(n)
    real(dl), allocatable :: integral(:)

    if (abs(a(size(a)) -1) > 1e-5) error stop 'w table must end at a=1'

    this%use_tabulated_w = .true.
    call this%equation_of_state%Init(log(a), w)

    allocate(integral(this%equation_of_state%n))
    ! log (rho) =  -3 int dlna (1+w)
    call this%equation_of_state%IntegralArray(integral)
    integral  = -3*( (this%equation_of_state%X-this%equation_of_state%X(1)) + integral) + 4*this%equation_of_state%X
    integral = integral - integral(this%equation_of_state%n) !log(a^4 rho_de)) normalized to 0 at a=1
    call this%logdensity%Init(this%equation_of_state%X, integral)
    !Set w and wa to values today (e.g. as too simple first guess for approx fittings etc).
    this%w_lam = w(size(a))
    this%wa = -this%equation_of_state%Derivative(0._dl)

    end subroutine TDarkEnergyEqnOfState_SetwTable

    ! EQUATIONS OF STATE

    function w_lcdm(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0
    end function w_lcdm

    function w_sine(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (1.0_dl + pm * sin(1.0_dl - a))
    end function w_sine
    
    function w_ba(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (1.0_dl + pm * (1.0_dl - a) / (a**2 + (1.0_dl - a)**2))
    end function w_ba
    
    function w_sine_osc0(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (1.0_dl + pm * ( - a * sin(1.0_dl / a) + sin(1.0_dl)))
    end function w_sine_osc0
    
    function w_sine_osc1(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
        if(a < 1) then
            res = w0 * (1.0_dl + pm * (1.0_dl - a) * sin(1.0_dl / (1.0_dl - a)))
        else
            res = w0
        endif
    end function w_sine_osc1
    
    function w_neg1_sine(a, pm, C) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, C
        real(dl) :: res
    
        res = -1.0_dl + pm * (1.0_dl / C) * sin(C * (1.0_dl - a))
    end function w_neg1_sine
    
    function w_neg1_sine_incr(a, pm, C) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, C
        real(dl) :: res
    
        res = -1.0_dl + pm * (a / C) * sin(C * (1.0_dl - a))
    end function w_neg1_sine_incr
    
    function w_neg1_sine_decr(a, pm, C) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, C
        real(dl) :: res
    
        res = -1.0_dl + pm * (1.0_dl - a) / C * sin(C * a)
    end function w_neg1_sine_decr

    function w_cases(dde_model, a, pm, param) result(w_value)
        implicit none
        integer, intent(in) :: dde_model, pm
        real(dl), intent(in) :: a, param
        real(dl) :: w_value
    
        select case(dde_model)
            case (0)
                w_value = w_lcdm(a, pm, param)
            case (1)
                w_value = w_sine(a, pm, param)
            case (2)
                w_value = w_ba(a, pm, param)
            case (3)
                w_value = w_sine_osc0(a, pm, param)
            case (4)
                w_value = w_sine_osc1(a, pm, param)
            case (5)
                w_value = w_neg1_sine(a, pm, param)
            case (6)
                w_value = w_neg1_sine_incr(a, pm, param)
            case (7)
                w_value = w_neg1_sine_decr(a, pm, param)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
        end select
    end function w_cases

    function w0wa_sine_osc0(a, w0_val, wa_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: res
        
        res = w0_val + wa_val * (- a * sin(1.0_dl / a) + sin(1.0_dl))
    end function w0wa_sine_osc0

    function w0wa_sine_osc1(a, w0_val, wa_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: res
    
        res = w0_val + wa_val * (1.0_dl - a) * sin(1.0_dl / (1.0_dl - a))
    end function w0wa_sine_osc1

    function w0wa_cases(dde_model, a, w0_val, wa_val) result(w_value)
        implicit none
        integer, intent(in) :: dde_model
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: w_value
    
        select case(dde_model)
            case (8)
                w_value = w0_val + wa_val*(1._dl - a)
            case (9)
                w_value = w0wa_sine_osc0(a, w0_val, wa_val)
            case (10)
                w_value = w0wa_sine_osc1(a, w0_val, wa_val)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
        end select
    end function w0wa_cases

    function w0wa_sine_osc1_freq(a, w0_val, wa_val, freq_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val, freq_val
        real(dl) :: res
    
        res = w0_val + wa_val * (1.0_dl - a) * sin(freq_val * 1.0_dl / (1.0_dl - a))
    end function w0wa_sine_osc1_freq

    function w0wa_freq_cases(dde_model, a, w0_val, wa_val, freq_val) result(w_value)
        implicit none
        integer, intent(in) :: dde_model
        real(dl), intent(in) :: a, w0_val, wa_val, freq_val
        real(dl) :: w_value
    
        select case(dde_model)
            case (11)
                w_value = w0wa_sine_osc1_freq(a, w0_val, wa_val, freq_val)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
        end select
    end function w0wa_freq_cases

    function TDarkEnergyEqnOfState_w_de(this, a)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: TDarkEnergyEqnOfState_w_de, al
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        if (this%dde_model < 8) then
            TDarkEnergyEqnOfState_w_de = w_cases(this%dde_model, a, this%pm_in, this%w_lam)
        elseif (this%dde_model < 11) then
            TDarkEnergyEqnOfState_w_de = w0wa_cases(this%dde_model, a, this%w_lam, this%wa)
        else
            TDarkEnergyEqnOfState_w_de = w0wa_freq_cases(this%dde_model, a, this%w_lam, this%wa, this%dde_freq)
        endif
    else
        al=dlog(a)
        if(al <= this%equation_of_state%Xmin_interp) then
            TDarkEnergyEqnOfState_w_de= this%equation_of_state%F(1)
        elseif(al >= this%equation_of_state%Xmax_interp) then
            TDarkEnergyEqnOfState_w_de= this%equation_of_state%F(this%equation_of_state%n)
        else
            TDarkEnergyEqnOfState_w_de = this%equation_of_state%Value(al)
        endif
    endif

    end function TDarkEnergyEqnOfState_w_de  ! equation of state of the PPF DE


    subroutine TDarkEnergyEqnOfState_Effective_w_wa(this, w, wa)
    class(TDarkEnergyEqnOfState), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    w = this%w_lam
    wa = this%wa

    end subroutine TDarkEnergyEqnOfState_Effective_w_wa

    ! Helper function for Si, Ci evaluation; good to 1e-16
    ! Source: https://arxiv.org/pdf/1407.7676
    function f_pade(x) result(f)
        implicit none
        real(dl), intent(in) :: x
        real(dl) :: y, f
      
        y = 1.0_dl / (x * x)
      
        f = (1.0_dl + &
             y * (7.44437068161936700618e2_dl + &
                  y * (1.96396372895146869801e5_dl + &
                       y * (2.37750310125431834034e7_dl + &
                            y * (1.43073403821274636888e9_dl + &
                                 y * (4.33736238870432522765e10_dl + &
                                      y * (6.40533830574022022911e11_dl + &
                                           y * (4.20968180571076940208e12_dl + &
                                                y * (1.00795182980368574617e13_dl + &
                                                     y * (4.94816688199951963482e12_dl + &
                                                          y * (-4.94701168645415959931e11_dl))))))))))) &
             / (x * (1.0_dl + &
                       y * (7.46437068161927678031e2_dl + &
                            y * (1.97865247031583951450e5_dl + &
                                 y * (2.41535670165126845144e7_dl + &
                                      y * (1.47478952192985464958e9_dl + &
                                           y * (4.58595115847765779830e10_dl + &
                                                y * (7.08501308149515401563e11_dl + &
                                                     y * (5.06084464593475076774e12_dl + &
                                                          y * (1.43468549171581016479e13_dl + &
                                                               y * (1.11535493509914254097e13_dl)))))))))))
      
    end function f_pade

    ! Helper function for Si, Ci evaluation; good to 1e-16
    ! Source: https://arxiv.org/pdf/1407.7676
    function g_pade(x) result(g)
        implicit none
        real(dl), intent(in) :: x
        real(dl) :: y, g
    
        y = 1.0_dl / (x * x)
    
        g = (y * (1.0_dl + &
            y * (8.1359520115168615e2 + &
                y * (2.35239181626478200e5 + &
                    y * (3.12557570795778731e7 + &
                        y * (2.06297595146763354e9 + &
                            y * (6.83052205423625007e10 + &
                                y * (1.09049528450362786e12 + &
                                    y * (7.57664583257834349e12 + &
                                        y * (1.81004487464664575e13 + &
                                            y * (6.43291613143049485e12 + &
                                                y * (-1.36517137670871689e12)))))))))))) / &
            (1.0_dl + &
                y * (8.19595201151451564e2 + &
                    y * (2.40036752835578777e5 + &
                        y * (3.26026661647090822e7 + &
                            y * (2.23355543278099360e9 + &
                                y * (7.87465017341829930e10 + &
                                    y * (1.39866710696414565e12 + &
                                        y * (1.17164723371736605e13 + &
                                            y * (4.01839087307656620e13 + &
                                                y * (3.99653257887490811e13))))))))))
    
    end function g_pade    

    ! Sine integral (https://mathworld.wolfram.com/SineIntegral.html) using Pade approximants; good to 1e-16
    ! Source: https://arxiv.org/pdf/1407.7676
    function Si_pade(x) result(res)
        implicit none
        real(dl), intent(in) :: x
        real(dl) :: res, x2
      
        ! Approximation for 0 <= x <= 4
        if (x <= 4.0_dl) then
          x2 = x * x
      
          res = x * (1.0_dl + &
                     x2 * (-4.54393409816329991e-2_dl + &
                           x2 * (1.15457225751016682e-3_dl + &
                                 x2 * (-1.41018536821330254e-5_dl + &
                                       x2 * (9.43280809438713025e-8_dl + &
                                             x2 * (-3.53201978997168357e-10_dl + &
                                                   x2 * (7.08240282274875911e-13_dl + &
                                                         x2 * (-6.05338212010422477e-16_dl)))))))) &
                / (1.0_dl + &
                   x2 * (1.01162145739225565e-2_dl + &
                         x2 * (4.99175116169755106e-5_dl + &
                               x2 * (1.55654986308745614e-7_dl + &
                                     x2 * (3.28067571055789734e-10_dl + &
                                           x2 * (4.5049097575386581e-13_dl + &
                                                 x2 * (3.21107051193712168e-16_dl)))))))
        ! Approximation for x > 4
        else
          res = (acos(-1.0_dl) / 2.0_dl) - f_pade(x) * cos(x) - g_pade(x) * sin(x)
        end if
      
    end function Si_pade

    ! Cosine integral (https://mathworld.wolfram.com/CosineIntegral.html) using Pade approximants; good to 1e-16
    ! Source: https://arxiv.org/pdf/1407.7676
    function Ci_pade(x) result(res)
        implicit none
        real(dl), intent(in) :: x
        real(dl) :: res, x2
        real(dl), parameter :: euler_gamma = 0.5772156649015329_dl  ! Euler-Mascheroni constant (16 decimal places)
      
        ! Approximation for 0 <= x <= 4
        if (x <= 4.0_dl) then
          x2 = x * x
      
          res = euler_gamma + log(x) + &
                x2 * (-0.25_dl + &
                      x2 * (7.51851524438898291e-3_dl + &
                            x2 * (-1.27528342240267686e-4_dl + &
                                  x2 * (1.05297363846239184e-6_dl + &
                                        x2 * (-4.68889508144848019e-9_dl + &
                                              x2 * (1.06480802891189243e-11_dl + &
                                                    x2 * (-9.93728488857585407e-15_dl))))))) &
                / (1.0_dl + &
                   x2 * (1.1592605689110735e-2_dl + &
                         x2 * (6.72126800814254432e-5_dl + &
                               x2 * (2.55533277086129636e-7_dl + &
                                     x2 * (6.97071295760958946e-10_dl + &
                                           x2 * (1.38536352772778619e-12_dl + &
                                                 x2 * (1.89106054713059759e-15_dl + &
                                                       x2 * (1.39759616731376855e-18_dl))))))))
        ! Approximation for x > 4
        else
          res = f_pade(x) * sin(x) - g_pade(x) * cos(x)
        end if
      
    end function Ci_pade
    
    ! EQUATION OF STATE INTEGRALS

    function Ia_lcdm(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * log(1.0_dl / a)
    end function Ia_lcdm

    function Ia_sine(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (log(1.0_dl / a) + pm * (sin(1.0_dl) * (Ci_pade(1.0_dl) - Ci_pade(a)) &
                                    - cos(1.0_dl) * (Si_pade(1.0_dl) - Si_pade(a))))
    end function Ia_sine
    
    function Ia_ba(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (log(1.0_dl / a) + pm * 0.5_dl * log(1.0_dl + (1.0_dl / a - 1.0_dl)**2))
    end function Ia_ba
    
    function Ia_sine_osc0(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res
    
        res = w0 * (log(1.0_dl / a) + pm * (sin(1.0_dl) * log(1.0_dl / a) &
                                    - (sin(1.0_dl) - a * sin(1.0_dl / a)) &
                                    + (Ci_pade(1.0_dl) - Ci_pade(1.0_dl / a))))
    end function Ia_sine_osc0
    
    function Ia_sine_osc1(a, pm, w0) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, w0
        real(dl) :: res

        if(a < 1) then
            res = w0 * (log(1.0_dl / a) + pm * (Ci_pade(1.0_dl / (1.0_dl - a)) &
                                    + Si_pade(1.0_dl / (1.0_dl - a)) &
                                    - sin(1.0_dl) * Ci_pade(a / (1.0_dl - a)) &
                                    - cos(1.0_dl) * Si_pade(a / (1.0_dl - a)) &
                                    - (1.0_dl - a) * sin(1.0_dl / (1.0_dl - a)) &
                                    - (0.5_dl * acos(-1.0_dl) * (1.0_dl - cos(1.0_dl)))))
        else
            res = 0 
        endif
    end function Ia_sine_osc1
    
    function Ia_neg1_sine(a, pm, C) result(res)
        implicit none
        real(dl), intent(in) :: a, C
        integer, intent(in) :: pm
        real(dl) :: res
    
        res = log(a) + pm * (1.0_dl / C) * (sin(C) * (Ci_pade(C) - Ci_pade(C * a)) &
                                            - cos(C) * (Si_pade(C) - Si_pade(C * a)))
    end function Ia_neg1_sine
    
    function Ia_neg1_sine_incr(a, pm, C) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, C
        real(dl) :: res
    
        res = log(a) + pm * (1.0_dl / C**2) * (1.0_dl - cos(C * (1.0_dl - a)))
    end function Ia_neg1_sine_incr
    
    function Ia_neg1_sine_decr(a, pm, C) result(res)
        implicit none
        integer, intent(in) :: pm
        real(dl), intent(in) :: a, C
        real(dl) :: res
    
        res = log(a) + pm * (1.0_dl / C**2) * (C * (Si_pade(C) - Si_pade(C * a)) &
                                            + cos(C) - cos(C * a))
    end function Ia_neg1_sine_decr

    function Ia_cases(dde_model, a, pm, param) result(Ia_value)
        implicit none
        integer, intent(in) :: dde_model, pm
        real(dl), intent(in) :: a, param
        real(dl) :: Ia_value
    
        select case(dde_model)
            case (0)
                Ia_value = Ia_lcdm(a, pm, param)
            case (1)
                Ia_value = Ia_sine(a, pm, param)
            case (2)
                Ia_value = Ia_ba(a, pm, param)
            case (3)
                Ia_value = Ia_sine_osc0(a, pm, param)
            case (4)
                Ia_value = Ia_sine_osc1(a, pm, param)
            case (5)
                Ia_value = Ia_neg1_sine(a, pm, param)
            case (6)
                Ia_value = Ia_neg1_sine_incr(a, pm, param)
            case (7)
                Ia_value = Ia_neg1_sine_decr(a, pm, param)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
        end select
    end function Ia_cases

    function Ia_sine_osc0_w0wa(a, w0_val, wa_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: res
        
        res = w0_val * log(1.0_dl / a) + wa_val * (sin(1.0_dl) * log(1.0_dl / a) &
                                    - (sin(1.0_dl) - a * sin(1.0_dl / a)) &
                                    + (Ci_pade(1.0_dl) - Ci_pade(1.0_dl / a)))
    end function Ia_sine_osc0_w0wa

    function Ia_sine_osc1_w0wa(a, w0_val, wa_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: res
    
        if(a < 1) then
            res = w0_val * log(1.0_dl / a) + wa_val * (Ci_pade(1.0_dl / (1.0_dl - a)) &
                                    + Si_pade(1.0_dl / (1.0_dl - a)) &
                                    - sin(1.0_dl) * Ci_pade(a / (1.0_dl - a)) &
                                    - cos(1.0_dl) * Si_pade(a / (1.0_dl - a)) &
                                    - (1.0_dl - a) * sin(1.0_dl / (1.0_dl - a)) &
                                    - (0.5_dl * acos(-1.0_dl) * (1.0_dl - cos(1.0_dl))))
        else
            res = 0 
        endif
    end function Ia_sine_osc1_w0wa

    function Ia_w0wa_cases(dde_model, a, w0_val, wa_val) result(Ia_value)
        implicit none
        integer, intent(in) :: dde_model
        real(dl), intent(in) :: a, w0_val, wa_val
        real(dl) :: Ia_value
    
        select case(dde_model)
            case (9)
                Ia_value = Ia_sine_osc0_w0wa(a, w0_val, wa_val)
            case (10)
                Ia_value = Ia_sine_osc1_w0wa(a, w0_val, wa_val)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
            end select
    end function Ia_w0wa_cases

    function Ia_sine_osc1_w0wa_freq(a, w0_val, wa_val, freq_val) result(res)
        implicit none
        real(dl), intent(in) :: a, w0_val, wa_val, freq_val
        real(dl) :: res
    
        if(a < 1) then
            res = w0_val * log(1.0_dl / a) + wa_val * (freq_val * Ci_pade(freq_val / (1.0_dl - a)) &
                                    + Si_pade(freq_val / (1.0_dl - a)) &
                                    - sin(freq_val) * Ci_pade(freq_val * a / (1.0_dl - a)) &
                                    - cos(freq_val) * Si_pade(freq_val * a / (1.0_dl - a)) &
                                    - (1.0_dl - a) * sin(freq_val / (1.0_dl - a)) &
                                    - (0.5_dl * acos(-1.0_dl) * (1.0_dl - cos(freq_val))))
        else
            res = 0 
        endif
    end function Ia_sine_osc1_w0wa_freq

    function Ia_w0wa_freq_cases(dde_model, a, w0_val, wa_val, freq_val) result(Ia_value)
        implicit none
        integer, intent(in) :: dde_model
        real(dl), intent(in) :: a, w0_val, wa_val, freq_val
        real(dl) :: Ia_value
    
        select case(dde_model)
            case (11)
                Ia_value = Ia_sine_osc1_w0wa_freq(a, w0_val, wa_val, freq_val)
            case default
                print *, 'Error: Unknown dde_model.'
                stop
            end select
    end function Ia_w0wa_freq_cases
      

    function TDarkEnergyEqnOfState_grho_de(this, a) result(grho_de) !relative density (8 pi G a^4 rho_de /grhov)
    class(TDarkEnergyEqnOfState) :: this
    real(dl) :: grho_de, al, fint
    real(dl), intent(IN) :: a

    if(.not. this%use_tabulated_w) then
        if (this%dde_model < 8) then
            grho_de = a ** (1._dl) * exp( 3. * Ia_cases(this%dde_model, a, this%pm_in, this%w_lam) )
        elseif (this%dde_model == 8) then
            grho_de = a ** (1._dl - 3. * this%w_lam - 3. * this%wa)
            if (this%wa/=0) grho_de=grho_de*exp(-3. * this%wa * (1._dl - a))
        elseif (this%dde_model < 11) then
            grho_de = a ** (1._dl) * exp( 3. * Ia_w0wa_cases(this%dde_model, a, this%w_lam, this%wa) )
        else
            grho_de = a ** (1._dl) * exp( 3. * Ia_w0wa_freq_cases(this%dde_model, a, this%w_lam, this%wa, this%dde_freq) )
        endif
    else
        if(a == 0.d0)then
            grho_de = 0.d0      !assume rho_de*a^4-->0, when a-->0, OK if w_de always <0.
        else
            if (a>=1) then
                fint= 1
            else
                al = dlog(a)
                if(al <= this%logdensity%X(1)) then
                    ! assume here w=w_de(a_min)
                    fint = exp(this%logdensity%F(1) + (1. - 3. * this%equation_of_state%F(1))*(al - this%logdensity%X(1)))
                else
                    fint = exp(this%logdensity%Value(al))
                endif
            end if
            grho_de = fint
        endif
    endif

    end function TDarkEnergyEqnOfState_grho_de

    subroutine TDarkEnergyEqnOfState_PrintFeedback(this, FeedbackLevel)
    class(TDarkEnergyEqnOfState) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa

    end subroutine TDarkEnergyEqnOfState_PrintFeedback

    subroutine TDarkEnergyEqnOfState_ReadParams(this, Ini)
    use IniObjects
    use FileUtils
    class(TDarkEnergyEqnOfState) :: this
    class(TIniFile), intent(in) :: Ini
    real(dl), allocatable :: table(:,:)

    this%use_tabulated_w = Ini%Read_Logical('use_tabulated_w', .false.)
    if(.not. this%use_tabulated_w)then
        this%w_lam = Ini%Read_Double('w', -1.d0)
        this%wa = Ini%Read_Double('wa', 0.d0)
        ! trap dark energy becoming important at high redshift 
        ! (will still work if this test is removed in some cases)
        if (this%w_lam + this%wa > 0) &
             error stop 'w + wa > 0, giving w>0 at high redshift'
    else
        call File%LoadTxt(Ini%Read_String('wafile'), table)
        call this%SetwTable(table(:,1),table(:,2), size(table(:,1)))
    endif

    end subroutine TDarkEnergyEqnOfState_ReadParams


    subroutine TDarkEnergyEqnOfState_Init(this, State)
    use classes
    class(TDarkEnergyEqnOfState), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_cosmological_constant = .not. this%use_tabulated_w .and. &
        &  abs(this%w_lam + 1._dl) < 1.e-6_dl .and. this%wa==0._dl

    end subroutine TDarkEnergyEqnOfState_Init


    end module DarkEnergyInterface
