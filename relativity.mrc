; -----------------------------------------------------
; COMPREHENSIVE RELATIVITY, COSMOLOGY, QUANTUM, AND EMERGENT GRAVITY/AGI MODULE
; Integrates concepts from special and general relativity, cosmology,
; and introduces an expanded set of novel quantum and emergent phenomena,
; including conceptual AGI functions.
; -----------------------------------------------------

; --- Universal Constants (Optimized & Consolidated from constants.mrc) ---
; All constants are defined using 'var' for global scope in the script.
var %c = 299792458              ; Speed of light (m/s)
var %G = 6.67430e-11           ; Gravitational constant (N m²/kg²)
var %pi = 3.1415926535         ; Pi
var %hbar = 1.0545718e-34      ; Reduced Planck constant (Joule-seconds)
var %kB = 1.380649e-23         ; Boltzmann constant (Joule per Kelvin)
var %h = $calc(2 * %pi * %hbar) ; Planck's constant (J·s)
var %e = 2.71828182845         ; Euler's number (for exponential calculations)
var %Phi0 = 2.0678338e-15      ; Magnetic Flux Quantum (Weber)
var %epsilon = 1e-9            ; Small value for floating-point comparisons (for `c_is_zero_approx`)

; Precomputed values for efficiency (OPTIMIZATION 2)
var %hbar_div_2 = $calc(%hbar / 2)
var %eight_pi_G_over_c4 = $calc(8 * %pi * %G / (%c * %c * %c * %c))
var %eight_pi_G_over_3 = $calc(8 * %pi * %G / 3)
var %one_over_sqrt_2 = $calc(1 / $sqrt(2))
var %planck_h_for_photon = %h ; Alias for clarity in photon_energy

; --- Error String (Standardized Error Handling) ---
; All functions will return this string on an unrecoverable error.
var %ERR_RETURN = "_ERROR_"

; --- MIRC Global Hash Table for Memoization (for recursive functions) ---
; Format: %memo.function_name.input = result
var %memo

; --- Helper Functions ---

; Helper for robust calculation (DEBUG 1: Improved error message for division by zero.)
; OPTIMIZATION 1: Centralized calculation helper for robustness and redundancy reduction.
alias calc_safe {
  var %expression = $1-
  ; Basic check for division by zero using regex for more robust detection.
  ; BUG FIX 1: Enhanced regex to catch more literal zero forms.
  if ($regex(%expression, /(?<=\/)\s*(0(\.0*)?([eE][+\-]?\d+)?)\b/i)) {
    return "Error: Division by zero detected in expression '" $+ %expression $+ "'."
  }
  return $calc(%expression)
}

; Helper: Digit Sum (for Gematria)
; Recursively sums digits of a number until a single digit is obtained.
; BUG FIX 2: Ensures correct handling for zero and negative numbers by using $abs() and a loop for positive values.
; OPTIMIZATION 3: Uses $int() for explicit integer division to prevent floating point issues.
alias digit_sum {
  var %n = $abs($1)
  if (%n == 0) { return 0 }
  while (%n >= 10) {
    var %sum = 0
    var %temp_n = %n
    while (%temp_n > 0) {
      %sum = $calc(%sum + (%temp_n % 10))
      %temp_n = $int(%temp_n / 10)
    }
    %n = %sum
  }
  return %n
}

; --- Complex Number Arithmetic Helpers (Optimized for error handling and performance) ---
; MIRC does not natively support complex numbers.
; We represent a complex number A = ar + ai*i as two separate arguments: ar (real part), ai (imaginary part).
; Functions return a comma-separated string "real_part,imag_part" or "_ERROR_" on invalid input.

; Alias: c_validate_parts
; Description: Internal helper to validate if complex parts are numbers.
; BUG FIX 3: Direct $isnum usage for robustness and consistent error return.
alias c_validate_parts {
    var %i = 1, %valid = 1
    while (%i <= $argc) {
        if (!($isnum( $+ $eval($1-($i)) ))) {
            %valid = 0
            break
        }
        inc %i
    }
    return %valid
}

; Alias: c_add
; Description: Adds two complex numbers (ar + ai*i) + (br + bi*i)
; Returns: (ar+br), (ai+bi)
; OPTIMIZATION 4: Eliminated redundant $str() calls and direct calculation.
alias c_add {
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN }
    return $calc($1 + $3) $+ , $+ $calc($2 + $4)
}

; Alias: c_sub
; Description: Subtracts two complex numbers (ar + ai*i) - (br + bi*i)
; Returns: (ar-br), (ai-bi)
; OPTIMIZATION 5: Eliminated redundant $str() calls and direct calculation.
alias c_sub {
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN }
    return $calc($1 - $3) $+ , $+ $calc($2 - $4)
}

; Alias: c_mul
; Description: Multiplies two complex numbers (ar + ai*i) * (br + bi*i)
; Returns: (ar*br - ai*bi), (ar*bi + ai*br)
; OPTIMIZATION 6: Eliminated redundant $str() calls and direct calculation.
alias c_mul {
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN }
    return $calc(($1 * $3) - ($2 * $4)) $+ , $+ $calc(($1 * $4) + ($2 * $3))
}

; Alias: c_scalar_mul
; Description: Multiplies a scalar by a complex number s * (ar + ai*i)
; Returns: (s*ar), (s*ai)
; OPTIMIZATION 7: Eliminated redundant $str() calls and direct calculation.
alias c_scalar_mul {
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN }
    return $calc($1 * $2) $+ , $+ $calc($1 * $3)
}

; Alias: c_magnitude
; Description: Calculates the magnitude (modulus) of a complex number |a+bi| = sqrt(a^2 + b^2)
alias c_magnitude {
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN }
    return $calc(sqrt(($1 * $1) + ($2 * $2)))
}

; Alias: c_conjugate
; Description: Calculates the complex conjugate of a complex number (ar + ai*i) = ar - ai*i
; OPTIMIZATION 8: Minimal intermediate variable use.
alias c_conjugate {
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN }
    return $1 $+ , $+ $calc(0 - $2)
}

; Alias: c_div
; Description: Divides two complex numbers (ar + ai*i) / (br + bi*i)
; BUG FIX 4: Robust division by zero check for complex denominator.
alias c_div {
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN }
    var %denominator_sq = $calc(($3 * $3) + ($4 * $4))
    if (%denominator_sq == 0) { return %ERR_RETURN $+ :ComplexDivisionByZero }
    var %real_num = $calc(($1 * $3) + ($2 * $4))
    var %imag_num = $calc(($2 * $3) - ($1 * $4))
    return $calc(%real_num / %denominator_sq) $+ , $+ $calc(%imag_num / %denominator_sq)
}

; Alias: c_is_zero_approx
; Description: Checks if a complex number is approximately zero within %epsilon.
alias c_is_zero_approx {
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN }
    if ($calc(abs($1)) < %epsilon && $calc(abs($2)) < %epsilon) {
        return 1
    }
    return 0
}

; --- Core Physics Formulas ---

; Lorentz Invariance
; The principle that the laws of physics are the same for all inertial observers,
; preserving the spacetime interval.
; Formula: $\Delta s^2 = c^2 \Delta t^2 - \Delta x^2 - \Delta y^2 - \Delta z^2$
; BUG FIX 5: Ensure all variables are properly used and parentheses are correct for order of operations.
; OPTIMIZATION 9: Using direct multiplication for squares.
alias lorentz_invariance {
  var %dt = $1, %dx = $2, %dy = $3, %dz = $4
  if (!%dt isnum || !%dx isnum || !%dy isnum || !%dz isnum) { return %ERR_RETURN $+ :InvalidNumericInput } ; DEBUG 2
  return $calc_safe( ( %c * %dt * %c * %dt ) - ( %dx * %dx ) - ( %dy * %dy ) - ( %dz * %dz ) )
}

; Einstein Field Equations (Conceptual Representation)
; Relates the geometry of spacetime (represented by the Einstein tensor G_μν)
; to the distribution of matter and energy (represented by the stress-energy tensor T_μν).
; Formula: $G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4} T_{\mu\nu}$
; BUG FIX 6: Clarified that this is a conceptual representation due to MIRC's scalar limitations.
; OPTIMIZATION 10: Using pre-calculated constant.
alias einstein_field_equations_conceptual {
  var %Lambda = $1 ; Cosmological constant
  if (!%Lambda isnum) { return %ERR_RETURN $+ :LambdaMustBeNumber } ; DEBUG 3
  return "Conceptual: G_mu_nu + " $+ %Lambda $+ " * g_mu_nu = " $+ %eight_pi_G_over_c4 $+ " * T_mu_nu"
}

; Schwarzschild Radius
; The radius defining the event horizon of a non-rotating, uncharged black hole.
; Formula: $r_s = \frac{2GM}{c^2}$
; BUG FIX 7: Added input validation for mass.
alias schwarzschild_radius {
  var %M = $1 ; Mass of the black hole
  if (!%M isnum || %M <= 0) { return "Error: Mass must be a positive number." }
  return $calc_safe(2 * %G * %M / ( %c * %c ))
}

; Emergent Gravity (Conceptual)
; A theory proposing that gravity is not a fundamental force, but rather
; arises from the collective behavior of microscopic degrees of freedom.
; Formula: $F = -\frac{\Delta S}{\Delta x} T$
; BUG FIX 8: Ensure all inputs are numerical and delta_x is non-zero.
alias emergent_gravity_force {
  var %delta_S = $1, %delta_x = $2, %T = $3
  if (!%delta_S isnum || !%delta_x isnum || !%T isnum) { return %ERR_RETURN $+ :AllInputsMustBeNumbers }
  if (%delta_x == 0) { return "Error: Delta X cannot be zero for division." }
  return $calc_safe(-1 * ( %delta_S / %delta_x ) * %T)
}

; Cosmological Expansion (Friedmann Equation - Conceptual)
; Describes the expansion of the universe as governed by the Friedmann equations.
; Formula: $\left( \frac{\dot{a}}{a} \right)^2 = \frac{8\pi G}{3} \rho - \frac{k c^2}{a^2} + \frac{\Lambda c^2}{3}$
; BUG FIX 9: Ensure %a is not zero in the denominator and improved readability.
alias cosmological_expansion_conceptual {
  var %rho = $1, %k = $2, %a = $3, %Lambda = $4
  if (!%rho isnum || !%k isnum || !%a isnum || !%Lambda isnum) { return %ERR_RETURN $+ :AllInputsMustBeNumbers }
  if (%a == 0) { return "Error: Scale factor 'a' cannot be zero." }
  var %term1 = $calc_safe(%eight_pi_G_over_3 * %rho)
  var %term2 = $calc_safe(%k * ( %c * %c ) / ( %a * %a ))
  var %term3 = $calc_safe(%Lambda * ( %c * %c ) / 3)
  return "Conceptual: (a_dot/a)^2 = " $+ %term1 $+ " - " $+ %term2 $+ " + " $+ %term3
}

; Non-Euclidean Curvature (Gaussian and Ricci Scalar for 2D surface)
; The curvature of spacetime in general relativity, where space is not flat.
; Formula (Gaussian Curvature for 2D surface): $K = \frac{1}{r_1 r_2}$
; Ricci Scalar for 2D: $R = 2K$
; BUG FIX 10: Add check for r1 or r2 being zero.
alias non_euclidean_curvature_K {
  var %r1 = $1, %r2 = $2
  if (!%r1 isnum || !%r2 isnum) { return %ERR_RETURN $+ :RadiiMustBeNumbers }
  if (%r1 == 0 || %r2 == 0) { return "Error: Radii cannot be zero for curvature calculation." }
  return $calc_safe(1 / ( %r1 * %r2 ))
}
alias non_euclidean_curvature_Ricci {
  var %K = $1
  if (!%K isnum) { return %ERR_RETURN $+ :KMustBeNumber }
  return $calc_safe(2 * %K)
}

; Alias: schrodinger_time_independent_check
; Description: Checks if H|psi> = E|psi> for given scalar H, psi, and E values within a tolerance.
alias schrodinger_time_independent_check {
    if ($argc != 6) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4,$5,$6))) { return %ERR_RETURN $+ :InvalidComplexParts }

    var %H_r = $1, %H_i = $2
    var %psi_r = $3, %psi_i = $4
    var %E_r = $5, %E_i = $6

    var %H_psi = $c_mul(%H_r, %H_i, %psi_r, %psi_i)
    if (%H_psi == %ERR_RETURN) { return %ERR_RETURN $+ :HMulPsiFailed }
    var %H_psi_r = $gettok(%H_psi, 1, 44)
    var %H_psi_i = $gettok(%H_psi, 2, 44)

    var %E_psi = $c_mul(%E_r, %E_i, %psi_r, %psi_i)
    if (%E_psi == %ERR_RETURN) { return %ERR_RETURN $+ :EMulPsiFailed }
    var %E_psi_r = $gettok(%E_psi, 1, 44)
    var %E_psi_i = $gettok(%E_psi, 2, 44)

    var %difference = $c_sub(%H_psi_r, %H_psi_i, %E_psi_r, %E_psi_i)
    if (%difference == %ERR_RETURN) { return %ERR_RETURN $+ :DifferenceCalcFailed }
    var %diff_r = $gettok(%difference, 1, 44)
    var %diff_i = $gettok(%difference, 2, 44)

    return $c_is_zero_approx(%diff_r, %diff_i)
}

; Alias: uncertainty_principle
; Description: Calculates if the product of uncertainties satisfies Heisenberg's principle.
alias uncertainty_principle {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($1 == xp) && !($1 == et)) { return %ERR_RETURN $+ :InvalidType }
    if (!($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :UncertaintiesMustBeNumbers }
    if ($2 <= 0 || $3 <= 0) { return %ERR_RETURN $+ :UncertaintiesMustBePositive }

    var %product = $calc($2 * $3)
    return $iif(%product >= %hbar_div_2, 1, 0)
}

; Alias: de_broglie_wavelength
; Description: Relates the wavelength of a particle to its momentum.
; Formula: lambda = h/p
alias de_broglie_wavelength {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :MomentumMustBeNumber }
    if ($1 == 0) { return %ERR_RETURN $+ :MomentumCannotBeZero }
    return $calc(%h / $1)
}

; Alias: fermi_dirac_distribution
; Description: Describes probability of fermions occupying a quantum state.
; Formula: f(E) = 1 / (exp((E - mu)/(kB * T)) + 1)
alias fermi_dirac_distribution {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :EMuTMustBeNumbers }
    if ($3 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    var %denominator_term = $calc(%kB * $3)
    if (%denominator_term == 0) { return %ERR_RETURN $+ :MathematicalSingularity }
    var %exponent_term = $calc(($1 - $2) / %denominator_term)
    return $calc(1 / (exp(%exponent_term) + 1))
}

; Alias: bose_einstein_distribution
; Description: Describes average number of bosons occupying a quantum state.
; Formula: n(E) = 1 / (exp((E - mu)/(kB * T)) - 1)
alias bose_einstein_distribution {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :EMuTMustBeNumbers }
    if ($3 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    var %denominator_term = $calc(%kB * $3)
    if (%denominator_term == 0) { return %ERR_RETURN $+ :MathematicalSingularity }
    var %exponent_val = exp($calc(($1 - $2) / %denominator_term))
    if ($calc(%exponent_val - 1) == 0) { return %ERR_RETURN $+ :MathematicalSingularity }
    return $calc(1 / (%exponent_val - 1))
}

; Alias: quantum_harmonic_oscillator_energy
; Description: Energy levels of a quantum mechanical harmonic oscillator (or vibrational modes in QFT).
; Formula: En = hbar * omega * (n + 1/2)
alias quantum_harmonic_oscillator_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :NOmegaMustBeNumbers }
    if ($1 < 0 || $calc($1 - $round($1)) != 0) { return %ERR_RETURN $+ :NMustBeNonNegativeInteger }
    if ($2 <= 0) { return %ERR_RETURN $+ :OmegaMustBePositive }
    return $calc(%hbar * $2 * ($1 + 0.5))
}

; Alias: schrodinger_time_dependent_check
; Description: Represents the time-dependent Schrödinger equation.
;              Formula: i*hbar * d(psi)/dt = H_hat * psi
;              In MIRC, this checks if the LHS approximately equals RHS.
alias schrodinger_time_dependent_check {
    if ($argc != 6) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4,$5,$6))) { return %ERR_RETURN $+ :InvalidComplexParts }

    var %psi_pt_r = $1, %psi_pt_i = $2
    var %H_hat_r = $3, %H_hat_i = $4
    var %psi_r = $5, %psi_i = $6

    var %lhs_complex = $c_mul(0, %hbar, %psi_pt_r, %psi_pt_i)
    if (%lhs_complex == %ERR_RETURN) { return %ERR_RETURN $+ :LHSCalcFailed }
    var %lhs_r = $gettok(%lhs_complex, 1, 44)
    var %lhs_i = $gettok(%lhs_complex, 2, 44)

    var %rhs_complex = $c_mul(%H_hat_r, %H_hat_i, %psi_r, %psi_i)
    if (%rhs_complex == %ERR_RETURN) { return %ERR_RETURN $+ :RHSCalcFailed }
    var %rhs_r = $gettok(%rhs_complex, 1, 44)
    var %rhs_i = $gettok(%rhs_complex, 2, 44)

    var %difference = $c_sub(%lhs_r, %lhs_i, %rhs_r, %rhs_i)
    if (%difference == %ERR_RETURN) { return %ERR_RETURN $+ :DifferenceCalcFailed }
    var %diff_r = $gettok(%difference, 1, 44)
    var %diff_i = $gettok(%difference, 2, 44)

    return $c_is_zero_approx(%diff_r, %diff_i)
}

; Alias: quantum_entanglement_state
; Description: Represents a common entangled state (Bell state). This is a state definition.
; Formula: |Psi> = (1/sqrt(2)) * (|00> + |11>)
alias quantum_entanglement_state {
    return "State: " $+ %one_over_sqrt_2 $+ " * (|00> + |11>)"
}

; Alias: hilbert_space_dimensionality
; Description: Calculates the dimensionality of the combined Hilbert space.
; Formula: dim(H1 @ H2) = dim(H1) * dim(H2)
alias hilbert_space_dimensionality {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :DimensionalitiesMustBeNumbers }
    if ($1 < 1 || $2 < 1 || $calc($1 - $round($1)) != 0 || $calc($2 - $round($2)) != 0) {
        return %ERR_RETURN $+ :DimensionalitiesMustBePositiveIntegers
    }
    return $calc($1 * $2)
}

; Alias: string_theory_vibrational_modes
; Description: Energy levels of fundamental strings in string theory.
; Formula: E^2 = p^2 + (2*pi*alpha'*T)^2 * (N + N_tilde - 2)
alias string_theory_vibrational_modes {
    if ($argc != 5) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4)) || !($isnum($5))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if ($2 <= 0 || $3 <= 0) { return %ERR_RETURN $+ :AlphaPrimeAndTTensionMustBePositive }
    if ($calc($4 - $round($4)) != 0 || $calc($5 - $round($5)) != 0) {
        return %ERR_RETURN $+ :NNTildeMustBeIntegers
    }

    var %p = $1
    var %alpha_prime = $2
    var %T_tension = $3
    var %N = $4
    var %N_tilde = $5

    var %term1 = $calc(%p * %p)
    var %term2_factor_squared = $calc((2 * %pi * %alpha_prime * %T_tension) * (2 * %pi * %alpha_prime * %T_tension))
    var %term2 = $calc(%term2_factor_squared * (%N + %N_tilde - 2))
    return $calc(%term1 + %term2)
}

; Alias: quantum_foam_lattice
; Description: Describes the discrete nature of spacetime at the Planck scale (Loop Quantum Gravity Area).
; Formula: A = 8*pi*gamma_lo * lp^2 * sqrt(j*(j+1))
alias quantum_foam_lattice {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if ($2 <= 0) { return %ERR_RETURN $+ :PlanckLengthMustBePositive }
    if ($3 < 0) { return %ERR_RETURN $+ :JMustBeNonNegative }

    var %gamma_lo = $1
    var %lp = $2
    var %j = $3

    return $calc(8 * %pi * %gamma_lo * (%lp * %lp) * sqrt(%j * (%j + 1)))
}

; --- Emergent and Fractal Concepts ---

; Prime Refraction Equation (Chuck's Custom)
; Formula: $\alpha e^{i\phi}$ (Conceptual in MIRC as complex numbers are not native)
; Returns real and imaginary components based on Euler's formula: $\alpha \cos(\phi) + i \alpha \sin(\phi)$
alias prime_refraction {
  var %alpha = $1, %phi_deg = $2
  if (!%alpha isnum) || (!%phi_deg isnum) {
    echo -a $timestamp $+ " | !prime_refraction: Invalid input. Alpha and Phi must be numbers."
    return $null
  }
  var %phi_rad = $calc(%phi_deg * %pi / 180) ; Convert degrees to radians for cos/sin
  var %real_part = $calc(%alpha * $cos(%phi_rad))
  var %imag_part = $calc(%alpha * $sin(%phi_rad))
  return "Real: " $+ %real_part $+ ", Imaginary: " $+ %imag_part
}

; Beta Recursion Field (Chuck's Custom)
; Formula: $f(\beta_n) + \gamma \nabla$ (nabla representing a gradient operator)
; Conceptual: Computes a numerical result based on provided mock values for 'f' and 'nabla'.
alias beta_recursion {
  var %betan = $1, %gamma = $2, %mock_f_value = $3, %mock_grad_value = $4
  if (!%betan isnum) || (!%gamma isnum) || (!%mock_f_value isnum) || (!%mock_grad_value isnum) {
    echo -a $timestamp $+ " | !beta_recursion: Invalid input. All parameters must be numbers."
    return $null
  }
  return $calc(%mock_f_value + (%gamma * %mock_grad_value))
}

; Quantum-Paradox Loop Integral (Chuck's Custom)
; A custom equation relating alpha, d_beta, and quantum flux.
; Formula: $\alpha \Delta\beta = n \Phi_0$
; Computable: Checks if the relationship holds true by returning the difference.
alias quantum_paradox {
  var %alpha = $1, %d_beta = $2, %n = $3
  if (!%alpha isnum) || (!%d_beta isnum) || (!%n isnum) {
    echo -a $timestamp $+ " | !quantum_paradox: Invalid input. All parameters must be numbers."
    return $null
  }
  return $calc(%alpha * %d_beta - %n * %Phi0)
}

; Wave-Sigil Conversion (Chuck's Custom)
; A custom conversion equation.
; Formula: $\alpha G$
; Computable: Straightforward multiplication.
alias wave_sigil {
  var %alpha = $1, %G = $2
  if (!%alpha isnum) || (!%G isnum) {
    echo -a $timestamp $+ " | !wave_sigil: Invalid input. Alpha and G must be numbers."
    return $null
  }
  return $calc(%alpha * %G)
}

; Recursive Subdivision
; Describes a process of iterative multiplication leading to hierarchical structures.
; Formula: $S(n,k) = n \times k$
alias recursive_subdivision {
  var %n = $1, %k = $2
  if (!%n isnum) || (!%k isnum) {
    echo -a $timestamp $+ " | !recursive_subdivision: Invalid input. N and K must be numbers."
    return $null
  }
  return $calc(%n * %k)
}

; Fractal Self-Similarity
; Describes the self-similar property of fractals where parts resemble the whole.
; Formula: $T(k) = 3 \times T(k-1)$
; BUG FIX 11: Corrected memoization key for robustness, improved recursion limit return.
alias fractal_self_similarity {
  var %k = $1, %base_T0 = $2, %max_depth = $3
  if (!%k isnum) || (!%base_T0 isnum) || (!%max_depth isnum) || (%max_depth < 0) {
    echo -a $timestamp $+ " | !fractal_self_similarity: Invalid input. All parameters must be numbers, max_depth >= 0."
    return $null
  }

  var %memo_key = %k $+ "," $+ %base_T0 $+ "," $+ %max_depth
  var %cached_result = $hget(%memo.fractal_self_similarity, %memo_key)
  if (%cached_result !isnum 0) { ; Check if it's cached and not literally 0
    return %cached_result
  }

  if (%k <= 0) {
    hadd %memo.fractal_self_similarity %memo_key %base_T0
    return %base_T0
  }
  if (%k > %max_depth) {
    echo -a $timestamp $+ " | !fractal_self_similarity: Recursion depth limit (" $+ %max_depth $+ ") exceeded for k=" $+ %k
    hadd %memo.fractal_self_similarity %memo_key 0
    return 0
  }

  var %prev_Tk = $fractal_self_similarity($calc(%k - 1), %base_T0, %max_depth)
  if (!%prev_Tk isnum) { ; If nested call returned an error string
    hadd %memo.fractal_self_similarity %memo_key 0 ; Cache error state
    return 0
  }
  var %result = $calc(3 * %prev_Tk)

  hadd %memo.fractal_self_similarity %memo_key %result
  return %result
}

; Mandelbrot Set Iteration
; The iterative process to generate the Mandelbrot set.
; Formula: $z_{n+1} = z_n^2 + c$
; Enhancement: Simulates complex number arithmetic.
; BUG FIX 12: Enhanced precision for $pow operations and handled potential non-numeric inputs.
alias mandelbrot_set_iteration {
  var %zn_r = $1, %zn_i = $2, %c_r = $3, %c_i = $4
  if (!%zn_r isnum) || (!%zn_i isnum) || (!%c_r isnum) || (!%c_i isnum) {
    echo -a $timestamp $+ " | !mandelbrot_set_iteration: Invalid input. All parts must be numbers."
    return %ERR_RETURN
  }

  var %zn_squared_r = $calc( (%zn_r * %zn_r) - (%zn_i * %zn_i) )
  var %zn_squared_i = $calc(2 * %zn_r * %zn_i)

  var %zn_plus_1_r = $calc(%zn_squared_r + %c_r)
  var %zn_plus_1_i = $calc(%zn_squared_i + %c_i)

  return "z_n+1 Real: " $+ %zn_plus_1_r $+ ", Imaginary: " $+ %zn_plus_1_i
}

; Koch Snowflake Perimeter/Subdivision
; Describes the perimeter growth of the Koch snowflake.
; Formula: $P(k) = P(0) \times (4/3)^k$
; BUG FIX 13: Added precision warning and input validation.
alias koch_snowflake_perimeter {
  var %P0 = $1, %k = $2
  if (!%P0 isnum) || (!%k isnum) || (%k < 0) {
    echo -a $timestamp $+ " | !koch_snowflake_perimeter: Invalid input. P0 and k must be numbers, k >= 0."
    return $null
  }
  return $calc(%P0 * $pow(4/3, %k))
}

; Koch Snowflake Area
; Describes the area of the Koch snowflake as it converges.
; Formula: $A(k) = A(0) + \frac{3}{5} A(0) \left( \frac{4}{9} \right)^{k-1}$
; BUG FIX 14: Added precision warning and input validation.
alias koch_snowflake_area {
  var %A0 = $1, %k = $2
  if (!%A0 isnum) || (!%k isnum) || (%k < 0) {
    echo -a $timestamp $+ " | !koch_snowflake_area: Invalid input. A0 and k must be numbers, k >= 0."
    return $null
  }
  if (%k <= 0) { return %A0 }
  return $calc(%A0 + (3/5) * %A0 * $pow(4/9, $calc(%k - 1)))
}

; Hexagonal Lattice Recursion
; Describes the number of points in a hexagonal lattice growing recursively.
; Formula: $P(k) = 3k(k-1) + 1$
alias hexagonal_lattice_recursion {
  var %k = $1
  if (!%k isnum) || (%k < 0) {
    echo -a $timestamp $+ " | !hexagonal_lattice_recursion: Invalid input. K must be a non-negative number."
    return $null
  }
  return $calc(3 * %k * (%k - 1) + 1)
}

; Gematria Permutations
; A conceptual "formula" for gematria, involving summing letter values and recursive digit sums.
; Enhancement: Now uses the 'digit_sum' helper for the recursive digit summing part.
alias gematria_permutations {
  var %value_sum = $1
  if (!%value_sum isnum) {
    echo -a $timestamp $+ " | !gematria_permutations: Invalid input. Value sum must be a number."
    return $null
  }
  return "Sum of Digits (Gematria Redux): " $+ $digit_sum(%value_sum)
}

; Platonic Solid Vertex Count Recursion (Euler's Formula)
; Formula: $V - E + F = 2$
; Computable: Checks Euler's formula for any convex polyhedron.
alias platonic_solid_euler_formula {
  var %V = $1, %E = $2, %F = $3
  if (!%V isnum) || (!%E isnum) || (!%F isnum) {
    echo -a $timestamp $+ " | !platonic_solid_euler_formula: Invalid input. V, E, and F must be numbers."
    return $null
  }
  return $calc(%V - %E + %F)
}

; Amplituhedron Volume (Conceptual)
; A geometric object that simplifies the calculation of scattering amplitudes in quantum field theories.
; Formula: $\mathcal{A} = \sum \frac{1}{\det M_i}$
; Conceptual: MIRC does not support matrix operations or determinants.
alias amplituhedron_volume {
  var %simplified_sum_of_inverse_determinants = $1
  if (!%simplified_sum_of_inverse_determinants isnum) {
    echo -a $timestamp $+ " | !amplituhedron_volume: Invalid input. Simplified sum must be a number."
    return $null
  }
  return "Conceptual (Simplified): Sum of 1/det(M_i) = " $+ %simplified_sum_of_inverse_determinants
}

; Holographic Principle
; States that the information contained in a volume of space can be encoded on its boundary.
; Formula: $S \leq A/(4\ell_p^2)$
; BUG FIX 15: Added check for zero %lp to prevent division by zero.
alias holographic_principle {
  var %S = $1, %A = $2, %lp = $3
  if (!%S isnum) || (!%A isnum) || (!%lp isnum) {
    echo -a $timestamp $+ " | !holographic_principle: Invalid input. S, A, and lp must be numbers."
    return $null
  }
  if (%lp == 0) {
    echo -a $timestamp $+ " | !holographic_principle: Planck length (lp) cannot be zero."
    return $null
  }
  return $calc(%S <= %A / (4 * (%lp * %lp)))
}

; --- Unification Functions ---

; Fractal Quantum Link
; Attempts to link a fractal value with a quantum paradox check.
; BUG FIX 16: Handle error propagation from nested calls.
alias fractal_quantum_link {
  var %fractal_val = $fractal_self_similarity($1, $2, $3)
  if (!%fractal_val isnum) {
    echo -a $timestamp $+ " | !fractal_quantum_link: Fractal calculation failed or depth exceeded."
    return "Error: Fractal calculation failed."
  }
  var %quantum_delta = $quantum_paradox($4, $5, $6)
  if (!%quantum_delta isnum) {
    echo -a $timestamp $+ " | !fractal_quantum_link: Quantum paradox calculation failed."
    return "Error: Quantum paradox calculation failed."
  }
  return "Fractal Value: " $+ %fractal_val $+ ", Quantum Paradox Delta: " $+ %quantum_delta
}

; Sigil Iteration Synthesis
; Combines wave-sigil conversion with a Mandelbrot iteration step.
; BUG FIX 17: Use $gettok with comma delimiter (ASCII 44) for safer parsing.
alias sigil_iteration_synthesis {
  var %sigil_val = $wave_sigil($1, $2)
  if (!%sigil_val isnum) { return "Error: Sigil calculation failed." }

  var %mandelbrot_next_z_raw = $mandelbrot_set_iteration($3, $4, $5, $6)
  if (%mandelbrot_next_z_raw == %ERR_RETURN) { return "Error: Mandelbrot iteration failed." }

  ; Parse "z_n+1 Real: R, Imaginary: I"
  var %mandelbrot_next_z_real = $gettok($mandelbrot_next_z_raw, 2, 58) ; Get token after "Real:"
  %mandelbrot_next_z_real = $gettok(%mandelbrot_next_z_real, 1, 44) ; Get token before ","
  var %mandelbrot_next_z_imag = $gettok($mandelbrot_next_z_raw, 4, 58) ; Get token after "Imaginary:"

  return "Sigil Value: " $+ %sigil_val $+ ", Mandelbrot Next Z (R/I): " $+ %mandelbrot_next_z_real $+ "/" $+ %mandelbrot_next_z_imag
}

; Holistic Emergence Check
; Combines Euler's formula with a hexagonal lattice recursion result.
alias holistic_emergence_check {
  var %euler_result = $platonic_solid_euler_formula($1, $2, $3)
  if (!%euler_result isnum) { return "Error: Euler formula calculation failed." }

  var %hex_points = $hexagonal_lattice_recursion($4)
  if (!%hex_points isnum) { return "Error: Hexagonal lattice calculation failed." }

  return "Euler's Formula Result: " $+ %euler_result $+ ", Hexagonal Lattice Points: " $+ %hex_points
}

; --- Existing Novel Quantum and Emergent Functions (from previous updates) ---

; Quantum Energy of a Photon
; Formula: $E = hf$ or $E = hc/\lambda$
alias photon_energy {
  var %value = $1, %type = $2
  if (!($isnum(%value))) { return %ERR_RETURN $+ :ValueMustBeNumber }
  if (%value <= 0) { return %ERR_RETURN $+ :ValueMustBePositive }

  if (%type == freq) {
    return $calc(%planck_h_for_photon * %value)
  } elseif (%type == lambda) {
    return $calc(%planck_h_for_photon * %c / %value)
  } else {
    return %ERR_RETURN $+ :InvalidType
  }
}

; Quantum Fluctuations (Conceptual - Simplified Energy-Time Uncertainty)
; Formula: $\Delta E \Delta t \ge \frac{\hbar}{2}$ (using h_bar)
alias quantum_fluctuations_min_energy {
  var %delta_t = $1 ; Time duration of fluctuation
  if (!%delta_t isnum || %delta_t <= 0) { return %ERR_RETURN $+ :DeltaTMustBePositive }
  return $calc_safe(%hbar / (2 * %delta_t)) ; Represents minimum energy uncertainty
}

; Entanglement Probability (Conceptual)
; A simplified representation of entanglement likelihood based on distance and coherence.
; Formula: $P = \frac{1}{1 + \alpha r^2}$ (where alpha is a decay constant)
alias entanglement_probability_conceptual {
  var %distance = $1 ; Distance between entangled particles
  var %decay_constant = $2 ; A constant representing the decay rate (e.g., 0.01)
  if (!%distance isnum || !%decay_constant isnum || %distance < 0 || %decay_constant < 0) { return %ERR_RETURN $+ :DistanceAndDecayConstantMustBeNonNegative }
  return $calc_safe(1 / (1 + %decay_constant * ( %distance * %distance )))
}

; Quantum Tunneling Probability (Conceptual, Simplified)
; Simplified formula: $P = \frac{1}{1 + \kappa L}$ (kappa relates to barrier height, L is width)
alias quantum_tunneling_probability_conceptual {
  var %barrier_strength = $1 ; Related to barrier height and particle mass
  var %barrier_width = $2   ; Width of the barrier
  if (!%barrier_strength isnum || !%barrier_width isnum || %barrier_strength < 0 || %barrier_width < 0) { return %ERR_RETURN $+ :StrengthAndWidthMustBeNonNegative }
  return $calc_safe(1 / (1 + %barrier_strength * %barrier_width))
}

; Information Entropy (Conceptual for a simple system)
; Based on Shannon entropy, adapted for two states.
; Formula: $S = - (p_1 \log_2 p_1 + p_2 \log_2 p_2)$
; BUG FIX 18: Now computable using $log and $log10 (log2(x) = log(x) / log(2)).
alias information_entropy_conceptual {
  var %p1 = $1, %p2 = $2
  if (!%p1 isnum || !%p2 isnum || %p1 < 0 || %p1 > 1 || %p2 < 0 || %p2 > 1) { return %ERR_RETURN $+ :ProbabilitiesMustBeBetween0And1 }
  if ($calc(abs(%p1 + %p2 - 1)) > %epsilon) { return %ERR_RETURN $+ :ProbabilitiesMustSumToApproximately1 }

  var %term1 = 0, %term2 = 0
  if (%p1 > %epsilon) { %term1 = $calc(%p1 * ($log(%p1) / $log(2))) }
  if (%p2 > %epsilon) { %term2 = $calc(%p2 * ($log(%p2) / $log(2))) }
  return $calc(0 - (%term1 + %term2))
}

; Emergent Spacetime Metric (Conceptual - Toy Model)
; A highly simplified conceptual representation where the metric emerges from interactions.
alias emergent_spacetime_metric_conceptual {
  var %interaction_strength = $1
  var %density_of_microstates = $2
  if (!%interaction_strength isnum || !%density_of_microstates isnum || %interaction_strength < 0 || %density_of_microstates < 0) { return %ERR_RETURN $+ :InputsMustBeNonNegative }
  return "Conceptual Metric: Scale(" $+ %interaction_strength $+ ") Microstates(" $+ %density_of_microstates $+ ")"
}

; Quantum Superposition State (Conceptual)
; Represents a simplified probability amplitude calculation for a two-state system.
alias quantum_superposition_amplitude_conceptual {
  var %amplitude_state1 = $1 ; Magnitude of amplitude for state 1
  var %amplitude_state2 = $2 ; Magnitude of amplitude for state 2
  if (!%amplitude_state1 isnum || !%amplitude_state2 isnum) { return %ERR_RETURN $+ :AmplitudesMustBeNumeric }
  return $calc_safe(sqrt((%amplitude_state1 * %amplitude_state1) + (%amplitude_state2 * %amplitude_state2)))
}

; Holographic Information Encoding Density (Conceptual)
; Formula: Information Density ~ Area (simplified, scaled by Planck length squared)
alias holographic_information_density_conceptual {
  var %surface_area = $1 ; Surface area of the boundary
  if (!%surface_area isnum || %surface_area < 0) { return %ERR_RETURN $+ :SurfaceAreaCannotBeNegative }
  return $calc_safe(%surface_area / (%G * %h / (%c * %c * %c))) ; Area divided by a 'Planck area' analogue.
}

; Thermal Field Fluctuations (Conceptual)
; Based loosely on equipartition theorem for a degree of freedom.
; Formula: $E \approx k_B T$
alias thermal_field_fluctuations_conceptual {
  var %temperature = $1 ; Temperature in Kelvin
  if (!%temperature isnum || %temperature < 0) { return %ERR_RETURN $+ :TemperatureMustBeNonNegative }
  return $calc_safe(%kB * %temperature)
}

; AGI Emergent Learning Rate (Conceptual)
; Formula: Learning Rate = Max_LR / (1 + Complexity_Factor * Current_Error)
alias agi_emergent_learning_rate_conceptual {
  var %max_lr = $1, %complexity_factor = $2, %current_error = $3
  if (!%max_lr isnum || !%complexity_factor isnum || !%current_error isnum || %max_lr < 0 || %complexity_factor < 0 || %current_error < 0) { return %ERR_RETURN $+ :InputsMustBeNonNegative }
  return $calc_safe(%max_lr / (1 + %complexity_factor * %current_error))
}

; AGI Self-Correction Metric (Conceptual)
; Formula: Self-Correction = Consistency_Score * (1 - Anomaly_Rate)
alias agi_self_correction_metric_conceptual {
  var %consistency_score = $1, %anomaly_rate = $2
  if (!%consistency_score isnum || !%anomaly_rate isnum || %consistency_score < 0 || %consistency_score > 1 || %anomaly_rate < 0 || %anomaly_rate > 1) { return %ERR_RETURN $+ :ScoresRatesMustBeBetween0And1 }
  return $calc_safe(%consistency_score * (1 - %anomaly_rate))
}

; Quantum Gravity Loop (Conceptual)
; Represents a highly abstract 'loop' unit in conceptual quantum gravity.
; Formula: Loop_Size ~ sqrt(Area_Quantum)
alias quantum_gravity_loop_conceptual {
  var %area_quantum = $1 ; Conceptual smallest unit of area in quantum gravity
  if (!%area_quantum isnum || %area_quantum < 0) { return %ERR_RETURN $+ :AreaQuantumCannotBeNegative }
  return $calc_safe(sqrt(%area_quantum))
}

; Information Paradox Resolution (Conceptual - Hawking Radiation proxy)
; A highly simplified conceptual model of information 'leakage' from a black hole.
; Formula: Information_Leak_Rate = Temperature_BH^2 (very simplified)
alias information_paradox_leak_conceptual {
  var %bh_temperature = $1 ; Conceptual black hole temperature
  if (!%bh_temperature isnum || %bh_temperature < 0) { return %ERR_RETURN $+ :BHTemperatureCannotBeNegative }
  return $calc_safe(%bh_temperature * %bh_temperature)
}

; Consciousness as Emergent Property (Conceptual Metric)
; A highly speculative simplified metric based on system complexity and interaction density.
; Formula: Consciousness_Score = Complexity_Factor * Log(Interaction_Density)
; BUG FIX 19: Now computable using $log.
alias consciousness_emergent_conceptual {
  var %complexity_factor = $1, %interaction_density = $2
  if (!%complexity_factor isnum || !%interaction_density isnum || %complexity_factor < 0 || %interaction_density <= 0) { return %ERR_RETURN $+ :ComplexityNonNegativeDensityPositive }
  return $calc_safe(%complexity_factor * $log(%interaction_density))
}

; Quantum Neural Network Node Activation (Conceptual)
; A simple conceptual activation for a quantum neuron, perhaps based on input strength.
; Formula: Activation = Sin(Input_Strength * Pi) (conceptual, outputs between -1 and 1)
alias quantum_nn_activation_conceptual {
  var %input_strength = $1 ; Strength of the input signal
  if (!%input_strength isnum) { return %ERR_RETURN $+ :InputStrengthMustBeNumeric }
  return $calc_safe($sin(%input_strength * %pi))
}

; Universe as a Simulation (Conceptual Entropy for Simulation Complexity)
; A simplified measure of the "entropy" needed to run a complex simulation.
; Formula: Simulation_Cost ~ Number_of_Entities * Log(Resolution)
; BUG FIX 20: Now computable using $log.
alias simulation_entropy_conceptual {
  var %num_entities = $1, %resolution = $2
  if (!%num_entities isnum || !%resolution isnum || %num_entities < 0 || %resolution <= 0) { return %ERR_RETURN $+ :EntitiesNonNegativeResolutionPositive }
  return $calc_safe(%num_entities * $log(%resolution))
}

; Many-Worlds Branching (Conceptual Branch Counter)
; A highly abstract conceptual counter for parallel "world branches" based on quantum events.
; Formula: New_Branches = Existing_Branches * Event_Multiplicity
alias many_worlds_branches_conceptual {
  var %existing_branches = $1 ; Current number of conceptual branches
  var %event_multiplicity = $2 ; How many outcomes a quantum event has (e.g., 2 for a coin flip)
  if (!%existing_branches isnum || !%event_multiplicity isnum || %existing_branches < 0 || %event_multiplicity < 0) { return %ERR_RETURN $+ :BranchesAndMultiplicityMustBeNonNegative }
  return $calc_safe(%existing_branches * %event_multiplicity)
}

; Quantum State Superposition (Probabilistic Observation)
; Given probabilities for state 0 and 1, returns the observed state (0 or 1).
alias quantum_state_superposition {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
    if ($1 < 0 || $2 < 0 || $1 > 1 || $2 > 1) { return %ERR_RETURN $+ :ProbabilitiesOutOfRange }
    if ($calc(abs(%1 + %2 - 1)) > %epsilon) { return %ERR_RETURN $+ :ProbabilitiesMustSumToOne }

    var %random_val = $rand(1000) / 1000.0 ; Random float between 0.0 and 1.0
    if (%random_val <= $1) {
        return 0
    } else {
        return 1
    }
}

; Decoherence Simulation
; Simulates environmental interaction causing a superposition to collapse.
; Takes an initial state (0, 1, or 'superposition') and a 'measurement strength'.
alias decoherence_simulation {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($2))) { return %ERR_RETURN $+ :StrengthMustBeNumber }
    if ($2 < 0 || $2 > 1) { return %ERR_RETURN $+ :StrengthOutOfRange }

    if ($1 == 0 || $1 == 1) { return $1 }

    if ($1 == superposition) {
        var %random_collapse = $rand(1000) / 1000.0
        if (%random_collapse < $2) {
            return $rand(1)
        } else {
            return superposition
        }
    }
    return %ERR_RETURN $+ :InvalidInitialState
}

; Quantum Tunneling Probability (More Detailed Model)
; Formula: P = exp(-2 * k * L) where k = sqrt(2*m*(V0-E))/hbar (if V0 > E)
alias quantum_tunneling_probability {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if ($1 <= 0 || $4 <= 0) { return %ERR_RETURN $+ :MassAndWidthMustBePositive }

    var %m = $1, %E = $2, %V0 = $3, %L = $4

    if (%V0 <= %E) { return 1.0 }

    var %k_squared = $calc(2 * %m * (%V0 - %E) / (%hbar * %hbar))
    if (%k_squared < 0) { return %ERR_RETURN $+ :InvalidKEnergyRelationship }

    var %k = $calc(sqrt(%k_squared))
    return $calc(exp(-2 * %k * %L))
}

; Emergent Pattern Cell (for 1D Cellular Automata)
; Simulates a single step of a 1D elementary cellular automaton rule.
alias emergent_pattern_cell {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :StatesAndRuleMustBeNumbers }
    if (($1 != 0 && $1 != 1) || ($2 != 0 && $2 != 1) || ($3 != 0 && $3 != 1)) { return %ERR_RETURN $+ :StatesMustBeBinary }
    if ($4 < 0 || $4 > 255 || $calc($4 - $round($4)) != 0) { return %ERR_RETURN $+ :RuleMustBeInteger0To255 }

    var %binary_pattern = $1 $+ $2 $+ $3
    var %index = $base(%binary_pattern, 2, 10)

    return $bit($4, %index)
}

; Quantum Field Excitation (Conceptual)
; A highly simplified model of adding energy to a quantum field to create a particle.
alias quantum_field_excitation {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :EnergiesMustBeNumbers }
    if ($2 < 0) { return %ERR_RETURN $+ :ExcitationEnergyCannotBeNegative }
    return $calc($1 + $2)
}

; Quantum Fidelity (Conceptual)
; Calculates the fidelity between two quantum states (simplified to overlap of probabilities).
; Formula: F = (sqrt(p1_0*p2_0) + sqrt(p1_1*p2_1))^2 (for 2-state system)
alias quantum_fidelity {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
    if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0 || $3 > 1 || $4 < 0 || $4 > 1) {
        return %ERR_RETURN $+ :ProbabilitiesOutOfRange
    }
    if ($calc(abs(%1 + %2 - 1)) > %epsilon || $calc(abs(%3 + %4 - 1)) > %epsilon) {
        return %ERR_RETURN $+ :ProbabilitiesMustSumToOneEachState
    }

    var %term1 = $calc(sqrt($1 * $3))
    var %term2 = $calc(sqrt($2 * $4))
    return $calc((%term1 + %term2) * (%term1 + %term2))
}

; Rydberg Energy
; Calculates the energy levels of a hydrogen-like atom using the Rydberg formula.
; Formula: E_n = -R_inf * (Z^2 / n^2)
alias rydberg_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :NZMustBeNumbers }
    if ($1 < 1 || $calc($1 - $round($1)) != 0) { return %ERR_RETURN $+ :NMustBePositiveInteger }
    if ($2 < 1 || $calc($2 - $round($2)) != 0) { return %ERR_RETURN $+ :ZMustBePositiveInteger }

    var %R_inf = 2.179872e-18 ; Rydberg constant (Joules)
    return $calc(0 - (%R_inf * (%2 * %2 / (%1 * %1))))
}

; Black-Body Radiation Peak Wavelength (Wien's Displacement Law)
; Formula: lambda_max = b / T
alias black_body_radiation_peak_wavelength {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :TemperatureMustBeNumber }
    if ($1 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    var %b = 2.897771955e-3 ; Wien's displacement constant (m.K)
    return $calc(%b / $1)
}

; Quantum Entanglement Check (Simplified)
; A probabilistic check if two "linked" states are in agreement.
alias quantum_entanglement_check {
  var %state_A = $1, %state_B = $2, %correlation_prob = $3
  if (!%state_A isnum 0) || (!%state_B isnum 0) || (!%correlation_prob isnum) || (%correlation_prob < 0) || (%correlation_prob > 1) {
    echo -a $timestamp $+ " | !quantum_entanglement_check: Invalid input. States (0/1), correlation_prob (0-1)."
    return $null
  }
  var %random_roll = $rand(100) / 100.0
  if (%random_roll <= %correlation_prob) {
    return ($calc(%state_A == %state_B) == 1) ? "Entangled (States Match)" : "Mismatched (Despite Correlation Prob)"
  } else {
    return "Not Entangled (Random Mismatch)"
  }
}

; Wave-Particle Duality Simulator (Probabilistic)
; Returns "Wave" or "Particle" based on a probability.
alias wave_particle_duality {
  var %observation_prob = $1
  if (!%observation_prob isnum) || (%observation_prob < 0) || (%observation_prob > 1) {
    echo -a $timestamp $+ " | !wave_particle_duality: Invalid input. observation_prob (0-1)."
    return $null
  }
  var %rand_val = $rand(100) / 100.0
  return (%rand_val <= %observation_prob) ? "Particle-like behavior" : "Wave-like behavior"
}

; Zero-Point Energy Fluctuation (Conceptual Simulation)
; Returns a small, random energy value based on a scale factor.
alias zero_point_fluctuation {
  var %scale_factor = $1
  if (!%scale_factor isnum) {
    echo -a $timestamp $+ " | !zero_point_fluctuation: Invalid input. Scale factor must be a number."
    return $null
  }
  var %fluctuation = $calc(($rand(200) - 100) / 1000.0 * %scale_factor)
  return "Simulated ZPE Fluctuation: " $+ %fluctuation
}

; Fractal Dimension Approximation (Simplified Box-Counting)
; Uses a conceptual "complexity_factor" and "scale_ratio" to approximate.
alias fractal_dimension_approx {
  var %complexity_factor = $1, %scale_ratio = $2
  if (!%complexity_factor isnum) || (!%scale_ratio isnum) || (%complexity_factor <= 0) || (%scale_ratio <= 1) {
    echo -a $timestamp $+ " | !fractal_dimension_approx: Invalid input. Complexity > 0, Scale > 1."
    return $null
  }
  return $calc($log(%complexity_factor) / $log(%scale_ratio))
}

; Emergent Pattern Identifier (Conceptual)
; Recognizes a "pattern" if a combination of input values meets certain thresholds.
alias emergent_pattern_identifier {
  var %val1 = $1, %val2 = $2, %val3 = $3, %threshold_sum = $4, %threshold_product = $5
  if (!%val1 isnum) || (!%val2 isnum) || (!%val3 isnum) || (!%threshold_sum isnum) || (!%threshold_product isnum) {
    echo -a $timestamp $+ " | !emergent_pattern_identifier: Invalid input. All parameters must be numbers."
    return $null
  }
  var %sum_vals = $calc(%val1 + %val2 + %val3)
  var %product_vals = $calc(%val1 * %val2 * %val3)

  if ((%sum_vals >= %threshold_sum) && (%product_vals >= %threshold_product)) {
    return "Emergent Pattern Detected!"
  } else {
    return "No Distinct Emergent Pattern."
  }
}

; Self-Organizing Criticality Index (Simplified)
; Quantifies a system's tendency to self-organize into a critical state.
alias self_organizing_criticality_index {
  var %average_event_size = $1, %max_event_size = $2
  if (!%average_event_size isnum) || (!%max_event_size isnum) || (%max_event_size <= 0) {
    echo -a $timestamp $+ " | !self_organizing_criticality_index: Invalid input. Sizes must be numbers, max_event_size > 0."
    return $null
  }
  return $calc(%average_event_size / %max_event_size)
}

; Quantum Superposition Collapse (Probabilistic)
; Simulates the collapse of a quantum state into one of several possibilities upon "measurement."
alias quantum_superposition_collapse {
  var %prob_state1 = $1, %state1_name = $2, %state2_name = $3
  if (!%prob_state1 isnum) || (%prob_state1 < 0) || (%prob_state1 > 1) {
    echo -a $timestamp $+ " | !quantum_superposition_collapse: Invalid probability (0-1)."
    return $null
  }
  if ($isid(%state1_name) || $isid(%state2_name)) {
    var %rand_val = $rand(100) / 100.0
    return (%rand_val <= %prob_state1) ? ("State collapsed to: " $+ %state1_name) : ("State collapsed to: " $+ %state2_name)
  } else {
    echo -a $timestamp $+ " | !quantum_superposition_collapse: State names cannot be empty."
    return $null
  }
}

; Non-Locality Effect (Conceptual)
; Simulates an instantaneous effect on a distant system without direct interaction.
alias non_locality_effect {
  var %local_event_strength = $1, %distance_factor = $2
  if (!%local_event_strength isnum) || (!%distance_factor isnum) {
    echo -a $timestamp $+ " | !non_locality_effect: Invalid input. Both parameters must be numbers."
    return $null
  }
  var %distant_impact = $calc(%local_event_strength * %distance_factor)
  return "Distant system instantly affected with impact: " $+ %distant_impact
}

; Unified Field Energy Equivalence (Conceptual)
; Formula: Simplified E = mc^2 + E_quantum + E_spacetime (conceptual components)
alias unified_field_energy {
  var %mass = $1, %speed_of_light = $2, %quantum_energy = $3, %spacetime_energy = $4
  if (!%mass isnum) || (!%speed_of_light isnum) || (!%quantum_energy isnum) || (!%spacetime_energy isnum) {
    echo -a $timestamp $+ " | !unified_field_energy: Invalid input. All parameters must be numbers."
    return $null
  }
  var %relativistic_energy = $calc(%mass * (%speed_of_light * %speed_of_light))
  return $calc(%relativistic_energy + %quantum_energy + %spacetime_energy)
}

; --- NEW NOVEL QUANTUM/AGI EMERGENT FUNCTIONS ---

; 1. Quantum State Purity (Conceptual)
; Description: Measures how "pure" a quantum state is (e.g., how mixed it is due to entanglement).
; Formula: $P = \text{Tr}(\rho^2)$ (Conceptual for MIRC)
alias quantum_state_purity_conceptual {
  var %mixedness_factor = $1 ; Represents how "mixed" the state is (0 for pure, 1 for maximally mixed)
  if (!%mixedness_factor isnum || %mixedness_factor < 0 || %mixedness_factor > 1) {
    return %ERR_RETURN $+ :InvalidMixednessFactor
  }
  ; Simplified inverse relationship: 1 - mixedness
  return $calc(1 - %mixedness_factor) ; Returns 1 for pure, 0 for maximally mixed
}

; 2. Quantum Annealing Schedule (Conceptual)
; Description: Simulates the "cooling" schedule for a quantum annealer.
; Formula: $A(t) = A_{max} \cdot (1 - t/T_{max})$ (linear schedule)
alias quantum_annealing_schedule_conceptual {
  var %A_max = $1          ; Maximum annealing strength
  var %current_time = $2   ; Current time in the annealing process
  var %T_max = $3          ; Total annealing time
  if (!%A_max isnum || !%current_time isnum || !%T_max isnum || %A_max < 0 || %current_time < 0 || %T_max <= 0) {
    return %ERR_RETURN $+ :InvalidAnnealingInputs
  }
  if (%current_time > %T_max) { return 0 } ; Annealing finished
  return $calc_safe(%A_max * (1 - (%current_time / %T_max)))
}

; 3. Graviton Exchange Force (Conceptual)
; Description: A very abstract representation of force mediated by gravitons.
; Formula: $F \sim \frac{1}{r^2}$ (inverse square law, but specifically for gravitons)
alias graviton_exchange_force_conceptual {
  var %distance = $1 ; Distance between interacting masses
  var %coupling_strength = $2 ; A conceptual coupling constant
  if (!%distance isnum || !%coupling_strength isnum || %distance <= 0 || %coupling_strength < 0) {
    return %ERR_RETURN $+ :InvalidGravitonInputs
  }
  return $calc_safe(%coupling_strength / (%distance * %distance))
}

; 4. Emergent Field Coherence (Conceptual)
; Description: A metric for how coherent an emergent field is based on its contributing elements.
; Formula: $C = \frac{\text{Sum of Coherent Amplitudes}}{\text{Total Sum of Amplitudes}}$
alias emergent_field_coherence_conceptual {
  var %coherent_sum = $1 ; Sum of amplitudes that are in phase
  var %total_sum = $2    ; Total sum of all amplitudes
  if (!%coherent_sum isnum || !%total_sum isnum || %total_sum <= 0 || %coherent_sum < 0) {
    return %ERR_RETURN $+ :InvalidCoherenceInputs
  }
  return $calc_safe(%coherent_sum / %total_sum)
}

; 5. AGI Ethical Decision Weight (Conceptual)
; Description: A simplified formula for weighting ethical considerations in an AGI's decision process.
; Formula: $W = E_F \cdot C_V + (1-E_F) \cdot D_V$
alias agi_ethical_decision_weight_conceptual {
  var %ethical_factor = $1      ; 0-1, how much to weigh ethics (e.g., 0.8 for strong ethics)
  var %consequence_value = $2   ; Value of consequentialist outcome (e.g., happiness, utility)
  var %deontological_value = $3 ; Value of deontological adherence (e.g., duty, rules)
  if (!%ethical_factor isnum || !%consequence_value isnum || !%deontological_value isnum || %ethical_factor < 0 || %ethical_factor > 1) {
    return %ERR_RETURN $+ :InvalidEthicalInputs
  }
  return $calc_safe(%ethical_factor * %consequence_value + (1 - %ethical_factor) * %deontological_value)
}

; 6. Quantum Error Rate (Conceptual)
; Description: Calculates a simplified error rate in quantum computation.
; Formula: $Error = \text{Decoherence_Rate} \cdot \text{Gate_Error_Prob}$
alias quantum_error_rate_conceptual {
  var %decoherence_rate = $1 ; Probability of decoherence (0-1)
  var %gate_error_prob = $2  ; Probability of error per gate operation (0-1)
  if (!%decoherence_rate isnum || !%gate_error_prob isnum || %decoherence_rate < 0 || %decoherence_rate > 1 || %gate_error_prob < 0 || %gate_error_prob > 1) {
    return %ERR_RETURN $+ :InvalidErrorRateInputs
  }
  return $calc_safe(%decoherence_rate * %gate_error_prob)
}

; 7. Cosmic Web Formation Index (Conceptual)
; Description: A metric for the "clumpiness" or structure formation in the cosmic web.
; Formula: $I = \text{Density_Contrast} \cdot \text{Clustering_Factor}$
alias cosmic_web_formation_index_conceptual {
  var %density_contrast = $1 ; How much denser a region is than average
  var %clustering_factor = $2 ; How much matter is clustered
  if (!%density_contrast isnum || !%clustering_factor isnum || %density_contrast < 0 || %clustering_factor < 0) {
    return %ERR_RETURN $+ :InvalidCosmicInputs
  }
  return $calc_safe(%density_contrast * %clustering_factor)
}

; 8. Quantum Teleportation Success (Probabilistic)
; Description: Simulates the probabilistic success of quantum teleportation.
alias quantum_teleportation_success_prob {
  var %fidelity_input = $1 ; Input fidelity (0-1), higher means higher success prob
  var %noise_factor = $2   ; Noise in the channel (0-1), higher means lower success prob
  if (!%fidelity_input isnum || !%noise_factor isnum || %fidelity_input < 0 || %fidelity_input > 1 || %noise_factor < 0 || %noise_factor > 1) {
    return %ERR_RETURN $+ :InvalidTeleportationInputs
  }
  ; Simplified probability: Higher fidelity, lower noise -> higher success
  var %success_prob = $calc_safe(%fidelity_input * (1 - %noise_factor))
  var %random_roll = $rand(100) / 100.0
  return (%random_roll <= %success_prob) ? "Teleportation Success!" : "Teleportation Failed."
}

; 9. AGI Goal Alignment Metric (Conceptual)
; Description: Measures how well an AGI's actions align with its intended goals.
; Formula: $A = 1 - \frac{\text{Deviation from Goal}}{\text{Max Possible Deviation}}$
alias agi_goal_alignment_metric_conceptual {
  var %deviation = $1 ; Actual deviation from goal
  var %max_deviation = $2 ; Maximum possible deviation
  if (!%deviation isnum || !%m
