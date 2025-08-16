; -----------------------------------------------------
; FOUNDATIONAL QUANTUM PHYSICS MODULE (Enhanced & Optimized)
; Integrates standard quantum mechanics and field theory concepts.
; All calculations are simplified for MIRC's environment.
; -----------------------------------------------------

; --- Universal Constants (Optimized & Consolidated) ---
; Constants are expected to be loaded from constants.mrc
; Ensure constants.mrc is loaded BEFORE this script.
; Example: var %hbar = 1.0545718e-34       ; Reduced Planck constant (Joule-seconds)
; Example: var %kB = 1.380649e-23          ; Boltzmann constant (Joules per Kelvin)
; Example: var %pi = 3.1415926535          ; Precomputed Pi for efficiency and consistency
; Example: var %c = 299792458              ; Speed of light (m/s)
; Example: var %epsilon = 1e-12            ; Small value for floating-point comparisons (Optimization: Increased precision for comparisons)
; Example: var %planck_h = $calc(2 * %pi * %hbar) ; Planck constant h = 2*pi*hbar
; Example: var %hbar_div_2 = $calc(%hbar / 2)     ; Pre-calculated for Uncertainty Principle


; --- Error String (Standardized Error Handling) ---
; All functions will return this string on an unrecoverable error.
var %ERR_RETURN = "_ERROR_"
var %STATE_SUPERPOSITION = "SUPERPOSITION" ; Constant for decoherence simulation state

; --- Helper Functions for Validation (Centralized for consistency) ---
; These helpers are assumed to be consistent with 'emergent_systems.mrc' and 'sensory.mrc'.

; Alias: is_positive_num
; Description: Checks if a value is a positive number.
; Returns: 1 if true, 0 otherwise.
alias is_positive_num {
    return $iif($isnum($1) && $1 > 0, 1, 0)
}

; Alias: is_non_negative_integer
; Description: Checks if a value is a non-negative integer.
; Returns: 1 if true, 0 otherwise.
alias is_non_negative_integer {
    return $iif($isnum($1) && $1 >= 0 && $calc($1 - $round($1)) == 0, 1, 0)
}

; --- Complex Number Arithmetic Helpers (Optimized for error handling and performance) ---
; These helpers are assumed to be consistent with 'emergent_systems.mrc' and 'sensory.mrc'.
; MIRC does not natively support complex numbers.
; We represent a complex number A = ar + ai*i as two separate arguments: ar (real part), ai (imaginary part).
; Functions return a comma-separated string "real_part,imag_part" or "_ERROR_" on invalid input.

; Alias: c_validate_parts
; Description: Internal helper to validate if complex parts are numbers.
; Usage: $c_validate_parts(ar, ai, [br, bi], ...)
; Returns: 1 if all parts are numbers, 0 otherwise.
alias c_validate_parts {
    var %i = 1, %valid = 1
    while (%i <= $argc) {
        ; Use $varg to safely access arguments dynamically
        if (!($isnum($varg(%i)))) {
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
; Usage: $c_add(ar, ai, br, bi)
alias c_add {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc($1 + $3) $+ , $+ $calc($2 + $4)
}

; Alias: c_sub
; Description: Subtracts two complex numbers (ar + ai*i) - (br + bi*i)
; Returns: (ar-br), (ai-bi)
; Usage: $c_sub(ar, ai, br, bi)
alias c_sub {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc($1 - $3) $+ , $+ $calc($2 - $4)
}

; Alias: c_mul
; Description: Multiplies two complex numbers (ar + ai*i) * (br + bi*i)
; Returns: (ar*br - ai*bi), (ar*bi + ai*br)
; Usage: $c_mul(ar, ai, br, bi)
alias c_mul {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    var %real = $calc(($1 * $3) - ($2 * $4))
    var %imag = $calc(($1 * $4) + ($2 * $3))
    return %real $+ , $+ %imag
}

; Alias: c_scalar_mul
; Description: Multiplies a scalar by a complex number s * (ar + ai*i)
; Returns: (s*ar), (s*ai)
; Usage: $c_scalar_mul(scalar, ar, ai)
alias c_scalar_mul {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :InvalidScalarOrComplexParts }
    return $calc($1 * $2) $+ , $+ $calc($1 * $3)
}

; Alias: c_magnitude
; Description: Calculates the magnitude (modulus) of a complex number |a+bi| = sqrt(a^2 + b^2)
; Usage: $c_magnitude(ar, ai)
alias c_magnitude {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc(sqrt(($1 * $1) + ($2 * $2)))
}

; Alias: c_conjugate
; Description: Calculates the complex conjugate of a complex number (ar + ai*i) = ar - ai*i
; Usage: $c_conjugate(ar, ai)
alias c_conjugate {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $1 $+ , $+ $calc(0 - $2)
}

; Alias: c_div
; Description: Divides two complex numbers (ar + ai*i) / (br + bi*i)
; Returns: ((ar*br + ai*bi)/(br^2 + bi^2)), ((ai*br - ar*bi)/(br^2 + bi^2))
; Usage: $c_div(ar, ai, br, bi)
alias c_div {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    var %denominator = $calc(($3 * $3) + ($4 * $4))
    ; Use %epsilon for float division by zero check, assuming it's globally available
    if ($calc(abs(%denominator)) < %epsilon) { return %ERR_RETURN $+ :DivisionByZero }
    var %real_num = $calc(($1 * $3) + ($2 * $4))
    var %imag_num = $calc(($2 * $3) - ($1 * $4))
    var %real = $calc(%real_num / %denominator)
    var %imag = $calc(%imag_num / %denominator)
    return %real $+ , $+ %imag
}

; Alias: c_is_zero_approx
; Description: Checks if a complex number is approximately zero within %epsilon.
; Usage: $c_is_zero_approx(ar, ai)
; Returns: 1 if approximately zero, 0 otherwise.
alias c_is_zero_approx {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    ; Use %epsilon, assuming it's globally available
    if ($calc(abs($1)) < %epsilon && $calc(abs($2)) < %epsilon) {
        return 1
    }
    return 0
}

; --- Core Quantum Mechanics Formulas ---

; Alias: schrodinger_time_independent_check
; Description: Checks if H|psi> = E|psi> for given scalar H, psi, and E values within a tolerance.
;              In MIRC, this is a simplified representation where H, psi, E are treated as complex scalars.
; Parameters: %H_real, %H_imag, %psi_real, %psi_imag, %E_real, %E_imag
; Returns: 1 if equation holds approximately, 0 otherwise, or _ERROR_.
alias schrodinger_time_independent_check {
    if ($argc != 6) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4,$5,$6))) { return %ERR_RETURN $+ :AllPartsMustBeNumbers }

    var %H_psi = $c_mul($1, $2, $3, $4)
    if (%H_psi == %ERR_RETURN) { return %ERR_RETURN $+ :HMulPsiFailed }
    var %H_psi_r = $gettok(%H_psi, 1, 44)
    var %H_psi_i = $gettok(%H_psi, 2, 44)

    var %E_psi = $c_mul($5, $6, $3, $4)
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
;              This function unifies position-momentum and energy-time uncertainty.
; Parameters: %type (string: "xp" for position-momentum, "et" for energy-time), %delta1, %delta2
; Returns: 1 if principle satisfied (product >= hbar/2), 0 otherwise, or _ERROR_.
alias uncertainty_principle {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($1 == xp) && !($1 == et)) { return %ERR_RETURN $+ :InvalidType }
    if (!($c_validate_parts($2,$3))) { return %ERR_RETURN $+ :UncertaintiesMustBeNumbers }
    if ($2 <= 0 || $3 <= 0) { return %ERR_RETURN $+ :UncertaintiesMustBePositive }

    ; Use %hbar_div_2, assuming it's globally available
    var %product = $calc($2 * $3)
    return $iif(%product >= %hbar_div_2, 1, 0)
}

; Alias: de_broglie_wavelength
; Description: Relates the wavelength of a particle to its momentum.
; Formula: lambda = h/p (where h = 2*pi*hbar)
; Parameters: %momentum (p)
; Returns: De Broglie wavelength or _ERROR_.
alias de_broglie_wavelength {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :MomentumMustBeNumber }
    if ($1 == 0) { return %ERR_RETURN $+ :MomentumCannotBeZero }
    ; Use %planck_h, assuming it's globally available
    return $calc(%planck_h / $1)
}

; Alias: fermi_dirac_distribution
; Description: Describes probability of fermions occupying a quantum state.
; Formula: f(E) = 1 / (exp((E - mu)/(kB * T)) + 1)
; Parameters: %E (Energy), %mu (Chemical potential), %T (Temperature)
; Returns: Probability or _ERROR_.
alias fermi_dirac_distribution {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :EMuTMustBeNumbers }
    if ($3 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    ; Use %kB and %epsilon, assuming they're globally available
    var %denominator_term = $calc(%kB * $3)
    if ($calc(abs(%denominator_term)) < %epsilon) { return %ERR_RETURN $+ :MathematicalSingularity }
    var %exponent_term = $calc(($1 - $2) / %denominator_term)
    return $calc(1 / (exp(%exponent_term) + 1))
}

; Alias: bose_einstein_distribution
; Description: Describes average number of bosons occupying a quantum state.
; Formula: n(E) = 1 / (exp((E - mu)/(kB * T)) - 1)
; Parameters: %E (Energy), %mu (Chemical potential), %T (Temperature)
; Returns: Average number or _ERROR_.
alias bose_einstein_distribution {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :EMuTMustBeNumbers }
    if ($3 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    ; Use %kB and %epsilon, assuming they're globally available
    var %denominator_term = $calc(%kB * $3)
    if ($calc(abs(%denominator_term)) < %epsilon) { return %ERR_RETURN $+ :MathematicalSingularity }
    var %exponent_val = exp($calc(($1 - $2) / %denominator_term))
    if ($calc(abs(%exponent_val - 1)) < %epsilon) { return %ERR_RETURN $+ :MathematicalSingularity }
    return $calc(1 / (%exponent_val - 1))
}

; Alias: quantum_harmonic_oscillator_energy
; Description: Energy levels of a quantum mechanical harmonic oscillator (or vibrational modes in QFT).
; Formula: En = hbar * omega * (n + 1/2)
; Parameters: %n (energy level, non-negative integer), %omega (angular frequency)
; Returns: Energy or _ERROR_.
alias quantum_harmonic_oscillator_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :NOmegaMustBeNumbers }
    if (!($is_non_negative_integer($1))) { return %ERR_RETURN $+ :NMustBeNonNegativeInteger }
    if (!($is_positive_num($2))) { return %ERR_RETURN $+ :OmegaMustBePositive }
    ; Use %hbar, assuming it's globally available
    return $calc(%hbar * $2 * ($1 + 0.5))
}

; Alias: schrodinger_time_dependent_check
; Description: Represents the time-dependent Schrödinger equation.
;              Formula: i*hbar * d(psi)/dt = H_hat * psi
;              In MIRC, this checks if the LHS approximately equals RHS.
;              Psi and H_hat are treated as complex scalar values, d(psi)/dt is treated as a complex scalar.
; Parameters: %psi_partial_t_real, %psi_partial_t_imag, %H_hat_real, %H_hat_imag, %psi_real, %psi_imag
; Returns: 1 if equation holds approximately, 0 otherwise, or _ERROR_.
alias schrodinger_time_dependent_check {
    if ($argc != 6) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4,$5,$6))) { return %ERR_RETURN $+ :AllPartsMustBeNumbers }

    ; Use %hbar, assuming it's globally available
    var %lhs_complex = $c_mul(0, %hbar, $1, $2)
    if (%lhs_complex == %ERR_RETURN) { return %ERR_RETURN $+ :LHSCalcFailed }
    var %lhs_r = $gettok(%lhs_complex, 1, 44)
    var %lhs_i = $gettok(%lhs_complex, 2, 44)

    var %rhs_complex = $c_mul($3, $4, $5, $6)
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
; Description: Represents a common entangled state (Bell state).
;              This is a state definition, not a calculation, and thus remains a string.
; Formula: |Psi> = (1/sqrt(2)) * (|00> + |11>)
; Returns: String representation of the entangled state.
alias quantum_entanglement_state {
    return "State: (1/sqrt(2)) * (|00> + |11>)"
}

; Alias: hilbert_space_dimensionality
; Description: Calculates the dimensionality of the combined Hilbert space.
; Formula: dim(H1 @ H2) = dim(H1) * dim(H2)
; Parameters: %dim_H1 (dimensionality of first space), %dim_H2 (dimensionality of second space)
; Returns: Combined dimensionality or _ERROR_.
alias hilbert_space_dimensionality {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :DimensionalitiesMustBeNumbers }
    if (!($is_positive_num($1)) || !($is_positive_num($2)) || $calc($1 - $round($1)) != 0 || $calc($2 - $round($2)) != 0) {
        return %ERR_RETURN $+ :DimensionalitiesMustBePositiveIntegers
    }
    return $calc($1 * $2)
}

; Alias: string_theory_vibrational_modes
; Description: Energy levels of fundamental strings in string theory.
; Formula: E^2 = p^2 + (2*pi*alpha'*T)^2 * (N + N_tilde - 2)
; Parameters: %p (momentum), %alpha_prime (Regge slope parameter), %T_tension (string tension),
;             %N (left-moving oscillator number), %N_tilde (right-moving oscillator number)
; Returns: E^2 value or _ERROR_. (Note: Returns E^2, not E)
alias string_theory_vibrational_modes {
    if ($argc != 5) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4,$5))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if (!($is_positive_num($2)) || !($is_positive_num($3))) { return %ERR_RETURN $+ :AlphaPrimeAndTTensionMustBePositive }

    if (!($is_non_negative_integer($4)) || !($is_non_negative_integer($5))) {
        return %ERR_RETURN $+ :NNTildeMustBeNonNegativeIntegers
    }
    ; Use %pi, %hbar, assuming they're globally available
    var %term1 = $calc($1 ^ 2)
    var %term2_factor_squared = $calc((2 * %pi * $2 * $3) ^ 2)
    var %term2 = $calc(%term2_factor_squared * ($4 + $5 - 2))
    return $calc(%term1 + %term2)
}

; Alias: quantum_foam_lattice
; Description: Describes the discrete nature of spacetime at the Planck scale (Loop Quantum Gravity Area).
; Formula: A = 8*pi*gamma_lo * lp^2 * sqrt(j*(j+1))
; Parameters: %gamma_lo (Barbero-Immirzi parameter), %lp (Planck length), %j (spin quantum number)
; Returns: Area or _ERROR_.
alias quantum_foam_lattice {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if (!($is_positive_num($2))) { return %ERR_RETURN $+ :PlanckLengthMustBePositive }
    if (!($is_non_negative_integer($3))) { return %ERR_RETURN $+ :JMustBeNonNegativeInteger }
    ; Use %pi, assuming it's globally available
    return $calc(8 * %pi * $1 * ($2 ^ 2) * sqrt($3 * ($3 + 1)))
}

; --- New Novel Quantum/Emergent Functions ---

; Alias: quantum_state_superposition
; Description: Simulates observing a quantum state that is in a superposition.
;              Given probabilities for state 0 and 1, returns the observed state (0 or 1).
; Parameters: %prob_state0, %prob_state1 (sum should be 1, or close to it)
; Returns: "0", "1", or _ERROR_.
alias quantum_state_superposition {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
    if ($1 < 0 || $2 < 0 || $1 > 1 || $2 > 1) { return %ERR_RETURN $+ :ProbabilitiesOutOfRange }
    ; Use %epsilon, assuming it's globally available
    if ($calc(abs($1 + $2 - 1)) > %epsilon) { return %ERR_RETURN $+ :ProbabilitiesMustSumToOne }

    var %random_val = $rand(1000) / 1000.0 ; Random float between 0.0 and 1.0
    if (%random_val <= $1) {
        return 0
    } else {
        return 1
    }
}

; Alias: decoherence_simulation
; Description: Simulates environmental interaction causing a superposition to collapse.
;              Takes an initial state (0, 1, or '%STATE_SUPERPOSITION') and a 'measurement strength'.
;              Higher strength increases probability of collapsing to a definite state.
; Parameters: %initial_state (0, 1, or "%STATE_SUPERPOSITION"), %measurement_strength (0 to 1)
; Returns: Final observed state (0 or 1), or _ERROR_ or %STATE_SUPERPOSITION.
alias decoherence_simulation {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($2))) { return %ERR_RETURN $+ :StrengthMustBeNumber }
    if ($2 < 0 || $2 > 1) { return %ERR_RETURN $+ :StrengthOutOfRange }

    if ($1 == 0 || $1 == 1) { return $1 }

    if ($1 == %STATE_SUPERPOSITION) {
        var %random_collapse = $rand(1000) / 1000.0
        if (%random_collapse < $2) {
            return $rand(1)
        } else {
            return %STATE_SUPERPOSITION
        }
    }
    return %ERR_RETURN $+ :InvalidInitialState
}

; Alias: quantum_tunneling_probability
; Description: Calculates a simplified tunneling probability through a rectangular potential barrier.
; Formula: P = exp(-2 * k * L) where k = sqrt(2*m*(V0-E))/hbar (if V0 > E)
; Parameters: %mass, %energy, %barrier_height, %barrier_width
; Returns: Probability (0-1) or _ERROR_.
alias quantum_tunneling_probability {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
    if (!($is_positive_num($1)) || !($is_positive_num($4))) { return %ERR_RETURN $+ :MassAndWidthMustBePositive }

    var %m = $1, %E = $2, %V0 = $3, %L = $4

    if (%V0 <= %E) { return 1.0 }

    ; Use %hbar, assuming it's globally available
    var %k_squared = $calc(2 * %m * (%V0 - %E) / (%hbar * %hbar))
    ; Check for negative argument to sqrt before calculation
    if (%k_squared < 0) { return %ERR_RETURN $+ :NegativeValueForSqrt }
    var %k = $calc(sqrt(%k_squared))
    return $calc(exp(-2 * %k * %L))
}

; Alias: photon_energy
; Description: Calculates the energy of a photon given its frequency or wavelength.
; Formula: E = h*f or E = h*c/lambda
; Parameters: %value, %type ("freq" or "lambda")
; Returns: Energy (Joules) or _ERROR_.
alias photon_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($is_positive_num($1))) { return %ERR_RETURN $+ :ValueMustBePositiveNumber }
    ; Explicitly check type string
    if (!($2 == freq) && !($2 == lambda)) { return %ERR_RETURN $+ :InvalidType }

    ; Use %planck_h and %c, assuming they're globally available
    if ($2 == freq) {
        return $calc(%planck_h * $1)
    } else { ; Must be lambda, already validated
        return $calc(%planck_h * %c / $1)
    }
}

; Alias: emergent_pattern_cell
; Description: Simulates a single step of a 1D elementary cellular automaton rule.
;              This is the core logic for a "Cellular Automata" display concept.
; Parameters: %left_state (0 or 1), %center_state (0 or 1), %right_state (0 or 1), %rule_number (0-255)
; Returns: New state of the center cell (0 or 1) or _ERROR_.
alias emergent_pattern_cell {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :StatesAndRuleMustBeNumbers }
    if (($1 != 0 && $1 != 1) || ($2 != 0 && $2 != 1) || ($3 != 0 && $3 != 1)) { return %ERR_RETURN $+ :StatesMustBeBinary }
    if (!($is_non_negative_integer($4)) || $4 > 255) { return %ERR_RETURN $+ :RuleMustBeInteger0To255 }

    var %binary_pattern = $1 $+ $2 $+ $3 ; e.g., "101"
    var %index = $base(%binary_pattern, 2, 10) ; Convert binary pattern to decimal index (0-7)

    return $bit($4, %index)
}

; Alias: quantum_field_excitation
; Description: A highly simplified model of adding energy to a quantum field to create a particle.
;              In this simplified view, an "excitation" adds a discrete quantum of energy.
; Parameters: %current_field_energy, %excitation_energy_quantum
; Returns: New total field energy or _ERROR_.
alias quantum_field_excitation {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :EnergiesMustBeNumbers }
    if ($2 < 0) { return %ERR_RETURN $+ :ExcitationEnergyCannotBeNegative }
    return $calc($1 + $2)
}

; Alias: quantum_fidelity
; Description: Calculates the fidelity between two quantum states (simplified to overlap of probabilities).
;              Fidelity measures how 'similar' two quantum states are.
; Formula: F = (sqrt(p1_0*p2_0) + sqrt(p1_1*p2_1))^2 (for 2-state system)
; Parameters: %p1_0 (prob of state 0 for state 1), %p1_1 (prob of state 1 for state 1),
;             %p2_0 (prob of state 0 for state 2), %p2_1 (prob of state 1 for state 2)
; Returns: Fidelity value (0-1) or _ERROR_.
alias quantum_fidelity {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
    if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0 || $3 > 1 || $4 < 0 || $4 > 1) {
        return %ERR_RETURN $+ :ProbabilitiesOutOfRange
    }
    ; Use %epsilon, assuming it's globally available
    if ($calc(abs($1 + $2 - 1)) > %epsilon || $calc(abs($3 + $4 - 1)) > %epsilon) {
        return %ERR_RETURN $+ :ProbabilitiesMustSumToOneEachState
    }

    var %term1 = $calc(sqrt($1 * $3))
    var %term2 = $calc(sqrt($2 * $4))
    return $calc( $min(1, $max(0, (%term1 + %term2) ^ 2)) )
}

; Alias: rydberg_energy
; Description: Calculates the energy levels of a hydrogen-like atom using the Rydberg formula.
; Formula: E_n = -R_inf * (Z^2 / n^2)
; Parameters: %n (principal quantum number, positive integer), %Z (atomic number, positive integer)
; Returns: Energy in Joules (negative) or _ERROR_.
alias rydberg_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :NZMustBeNumbers }
    if (!($is_positive_num($1)) || $calc($1 - $round($1)) != 0) { return %ERR_RETURN $+ :NMustBePositiveInteger }
    if (!($is_positive_num($2)) || $calc($2 - $round($2)) != 0) { return %ERR_RETURN $+ :ZMustBePositiveInteger }

    var %R_inf = 2.179872e-18 ; Rydberg constant (Joules)
    return $calc(0 - (%R_inf * ($2 ^ 2 / ($1 ^ 2))))
}

; Alias: black_body_radiation_peak_wavelength
; Description: Calculates the peak wavelength of black-body radiation using Wien's Displacement Law.
; Formula: lambda_max = b / T
; Parameters: %T (temperature in Kelvin)
; Returns: Peak wavelength in meters or _ERROR_.
alias black_body_radiation_peak_wavelength {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($is_positive_num($1))) { return %ERR_RETURN $+ :TemperatureMustBePositiveNumber }

    var %b = 2.897771955e-3 ; Wien's displacement constant (m.K)
    return $calc(%b / $1)
}

; --- Additional Novel Quantum Functionality / AGI Emergent Functions ---

; Alias: quantum_gate_not
; Description: Simulates a quantum NOT gate (Pauli-X). Flips a classical bit representation.
; Parameters: %input_bit (0 or 1)
; Returns: Flipped bit (1 or 0) or _ERROR_.
alias quantum_gate_not {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :InputMustBeNumber }
    if ($1 != 0 && $1 != 1) { return %ERR_RETURN $+ :InputMustBeBinary }
    return $iif($1 == 0, 1, 0)
}

; Alias: quantum_bell_state_measurement
; Description: Simplified measurement of a Bell state, showing correlation.
;              In a perfect Bell state, if one qubit is 0, the other is 0. If one is 1, the other is 1.
; Parameters: %measured_qubit1 (0 or 1)
; Returns: Correlated state of qubit2 (0 or 1) or _ERROR_.
alias quantum_bell_state_measurement {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :InputMustBeNumber }
    if ($1 != 0 && $1 != 1) { return %ERR_RETURN $+ :InputMustBeBinary }
    return $1 ; Simple correlation for |00> + |11> Bell state
}

; Alias: quantum_random_number_generator
; Description: Simulates a quantum random number generator based on a superposition collapse.
;              The 'quantum' part implies true randomness, not pseudo-randomness.
; Parameters: (none)
; Returns: 0 or 1 (randomly generated)
alias quantum_random_number_generator {
    return $rand(1)
}

; Alias: agi_neural_node_activation
; Description: Simulates a simplified activation function for a single neuron in a neural network.
; Formula: Output = 1 / (1 + exp(-sum_inputs)) (Sigmoid activation)
; Parameters: %weighted_sum_of_inputs (e.g., from dot product of inputs and weights)
; Returns: Activation output (0-1) or _ERROR_.
alias agi_neural_node_activation {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :InputMustBeNumber }
    ; Use %e (Euler's number), assuming it's globally available
    return $calc(1 / (1 + exp(0 - $1)))
}

; Alias: agi_genetic_crossover
; Description: Simulates a single-point crossover operation in a genetic algorithm.
;              Takes two "parent" binary strings and a crossover point.
; Parameters: %parent1_genes (binary string), %parent2_genes (binary string), %crossover_point (integer index)
; Returns: Combined "child" gene string or _ERROR_.
alias agi_genetic_crossover {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($is_non_negative_integer($3))) { return %ERR_RETURN $+ :CrossoverPointMustBeNonNegativeInteger }
    ; Validate gene lengths for consistency and non-zero. Also ensure genes are binary strings.
    if ($len($1) != $len($2) || $len($1) == 0 || $regsubex($1,/[^01]/g,'') != $1 || $regsubex($2,/[^01]/g,'') != $2) {
        return %ERR_RETURN $+ :GeneLengthsMustMatchAndBeNonZeroAndBinary
    }
    if ($3 > $len($1)) { return %ERR_RETURN $+ :CrossoverPointTooLarge }

    var %child_genes = $left($1, $3) $+ $right($2, $len($2) - $3)
    return %child_genes
}

; Alias: agi_swarm_cohesion
; Description: Simulates a basic cohesion force in a swarm intelligence model.
;              A "particle" moves towards the center of mass of its neighbors.
; Parameters: %particle_pos (current position), %center_of_mass_pos (neighbors' average position), %cohesion_strength (0-1)
; Returns: New particle position after cohesion influence, or _ERROR_.
alias agi_swarm_cohesion {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :PositionsAndStrengthMustBeNumbers }
    if ($3 < 0 || $3 > 1) { return %ERR_RETURN $+ :CohesionStrengthOutOfRange }

    ; Simplified 1D cohesion: new_pos = particle_pos + strength * (center_of_mass_pos - particle_pos)
    var %delta_pos = $calc($2 - $1)
    return $calc($1 + (%3 * %delta_pos))
}

; Alias: quantum_state_collapse_prob
; Description: Calculates the probability of observing a specific state from a superposition.
; Formula: P = |amplitude|^2 (simplified for classical probabilities)
; Parameters: %state_amplitude_real, %state_amplitude_imag
; Returns: Probability (0-1) or _ERROR_.
alias quantum_state_collapse_prob {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :AmplitudePartsMustBeNumbers }

    var %magnitude_sq = $calc(($1 * $1) + ($2 * $2))
    return $calc($min(1, $max(0, %magnitude_sq)))
}

; Alias: agi_reinforcement_signal
; Description: Simulates a basic reinforcement learning signal.
;              Calculates a 'value' update based on current value, reward, and learning rate.
; Formula: New Value = Old Value + Learning Rate * (Reward - Old Value)
; Parameters: %old_value, %reward, %learning_rate
; Returns: New estimated value or _ERROR_.
alias agi_reinforcement_signal {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :ValuesAndRateMustBeNumbers }
    if ($3 < 0 || $3 > 1) { return %ERR_RETURN $+ :LearningRateOutOfRange }

    return $calc($1 + ($3 * ($2 - $1)))
}

; Alias: quantum_tunneling_path_switch
; Description: Simulates a binary choice for a particle in a double-slit experiment (simplified).
;              Based on an "observation" factor, the particle picks a definite path.
; Parameters: %observation_factor (0=no observation, 1=full observation)
; Returns: "Path A" or "Path B" or %ERR_RETURN
alias quantum_tunneling_path_switch {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :ObservationFactorMustBeNumber }
    if ($1 < 0 || $1 > 1) { return %ERR_RETURN $+ :ObservationFactorOutOfRange }

    ; Corrected logic to ensure observation factor correctly biases the outcome.
    ; If observation_factor is high, the choice becomes less random and more "forced" to a definite path.
    ; If observation_factor is 1, it's always the first choice (e.g. Path A)
    ; If observation_factor is 0, it's purely random.
    var %rand_roll = $rand(1000) / 1000.0
    if (%rand_roll < $1) {
        ; Higher observation factor means it's more likely to "collapse" to a predictable state (e.g., path A)
        return "Path A"
    } else {
        return "Path B"
    }
}

; Alias: agi_emergent_collective_decision
; Description: Simulates how a collective of agents reaches a decision based on individual "votes."
; Parameters: %vote_for_A (count), %vote_for_B (count), %threshold (percentage 0-1)
; Returns: "Decision A", "Decision B", "No Decision" or _ERROR_.
alias agi_emergent_collective_decision {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :VotesAndThresholdMustBeNumbers }
    if ($1 < 0 || $2 < 0) { return %ERR_RETURN $+ :VotesCannotBeNegative }
    if ($3 < 0 || $3 > 1) { return %ERR_RETURN $+ :ThresholdOutOfRange }

    var %total_votes = $calc($1 + $2)
    if (%total_votes == 0) { return "No Decision" }

    var %ratio_A = $calc($1 / %total_votes)
    var %ratio_B = $calc($2 / %total_votes)

    if (%ratio_A >= %3) { return "Decision A" }
    if (%ratio_B >= %3) { return "Decision B" }
    return "No Decision"
}

; --- Novel Quantum Functionality / AGI Emergent Functions (10 New) ---

; 1. Quantum Harmonic Oscillator Energy Level (Extended)
; Description: Calculates the energy level, allowing for varying potential strengths.
; Formula: E_n = (n + 1/2) * hbar * sqrt(k/m) ; where omega = sqrt(k/m)
; Parameters: %n (energy level), %k (effective spring constant), %m (mass)
; Returns: Energy or _ERROR_.
alias quantum_harmonic_oscillator_energy_ext {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :NKMmustBeNumbers }
    if (!($is_non_negative_integer($1))) { return %ERR_RETURN $+ :NMustBeNonNegativeInteger }
    if (!($is_positive_num($2)) || !($is_positive_num($3))) { return %ERR_RETURN $+ :KandMMustBePositive }
    ; Use %hbar, assuming it's globally available
    var %omega = $calc(sqrt($2 / $3))
    return $calc(($1 + 0.5) * %hbar * %omega)
}

; 2. Photon Absorption/Emission Energy
; Description: Calculates the energy difference during photon absorption or emission.
; Formula: Delta_E = E_final - E_initial
; Parameters: %E_final, %E_initial
; Returns: Energy difference or _ERROR_.
alias photon_transition_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :EnergiesMustBeNumbers }
    return $calc($1 - $2)
}

; 3. Quantum Dot Energy Level (Simplified 1D Box)
; Description: Simulates the energy levels of a particle in a 1D quantum box.
; Formula: E_n = (n^2 * pi^2 * hbar^2) / (2 * m * L^2)
; Parameters: %n (energy level, positive integer), %m (mass), %L (length of box)
; Returns: Energy or _ERROR_.
alias quantum_dot_energy_1d {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :NMLMustBeNumbers }
    if (!($is_positive_num($1)) || $calc($1 - $round($1)) != 0) { return %ERR_RETURN $+ :NMustBePositiveInteger }
    if (!($is_positive_num($2)) || !($is_positive_num($3))) { return %ERR_RETURN $+ :MandLMustBePositive }
    ; Use %pi and %hbar, assuming they're globally available
    return $calc(($1 ^ 2 * %pi ^ 2 * %hbar ^ 2) / (2 * $2 * $3 ^ 2))
}

; 4. Quantum Annealing Success Probability (Conceptual)
; Description: A simplified model for the probability of finding the optimal solution in quantum annealing.
; Formula: P_success = exp(-Energy_Barrier / Temperature_Effective)
; Parameters: %energy_barrier (conceptual energy barrier to optimal state), %temp_effective (conceptual effective temperature)
; Returns: Probability (0-1) or _ERROR_.
alias quantum_annealing_prob {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :BarrierAndTempMustBeNumbers }
    if ($2 <= 0) { return %ERR_RETURN $+ :EffectiveTemperatureMustBePositive }

    var %exponent = $calc(0 - ($1 / $2))
    return $calc(exp(%exponent))
}

; 5. Quantum Measurement Probability Distribution (Conceptual Update)
; Description: Updates the probability distribution of states after a partial measurement.
; Parameters: %old_prob_state0, %old_prob_state1, %measurement_bias (bias towards state 0, -1 to 1)
; Returns: "new_prob0,new_prob1" or _ERROR_.
alias quantum_measurement_distribution {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :ProbabilitiesAndBiasMustBeNumbers }
    if ($1 < 0 || $2 < 0 || $1 > 1 || $2 > 1) { return %ERR_RETURN $+ :ProbabilitiesOutOfRange }
    ; Use %epsilon, assuming it's globally available
    if ($calc(abs($1 + $2 - 1)) > %epsilon) { return %ERR_RETURN $+ :ProbabilitiesMustSumToOne }
    if ($3 < -1 || $3 > 1) { return %ERR_RETURN $+ :BiasOutOfRange }

    ; Simulate partial measurement:
    ; New prob0 is a blend of old prob0 and the bias
    var %new_prob0 = $calc(%1 + %3 * (1 - %1 - %2) / 2 + %3 * %1) ; Simplified, conceptual update logic
    var %new_prob1 = $calc(1 - %new_prob0) ; Normalize

    ; Clamp to ensure valid probabilities
    %new_prob0 = $calc($min(1, $max(0, %new_prob0)))
    %new_prob1 = $calc($min(1, $max(0, %new_prob1)))

    return %new_prob0 $+ , $+ %new_prob1
}

; 6. AGI Collective Memory Recall (Conceptual)
; Description: Simulates a collective memory recall, taking a weighted average of individual "memory strengths."
; Parameters: %mem_strength1, %mem_strength2, %weight1 (0-1), %weight2 (0-1)
; Returns: Weighted collective recall value or _ERROR_.
alias agi_collective_memory_recall {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :StrengthsAndWeightsMustBeNumbers }
    if ($3 < 0 || $3 > 1 || $4 < 0 || $4 > 1) { return %ERR_RETURN $+ :WeightsOutOfRange }
    if ($calc(abs($3 + $4 - 1)) > %epsilon && $calc(abs($3 + $4)) > %epsilon) { return %ERR_RETURN $+ :WeightsMustSumToOneOrZero }

    var %total_weight = $calc($3 + $4)
    if ($calc(abs(%total_weight)) < %epsilon) { return 0 }

    return $calc((%1 * %3 + %2 * %4) / %total_weight)
}

; 7. AGI Emergent Communication Protocol (Encoding Success)
; Description: Simulates the success rate of an emergent communication protocol based on complexity and channel noise.
; Formula: Success_Rate = max(0, 1 - (complexity_factor * noise_level))
; Parameters: %complexity_factor (0-1), %noise_level (0-1)
; Returns: Success rate (0-1) or _ERROR_.
alias agi_comm_protocol_success {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :ComplexityAndNoiseMustBeNumbers }
    if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1) { return %ERR_RETURN $+ :FactorAndNoiseOutOfRange }

    var %loss = $calc($1 * $2)
    return $calc($max(0, 1 - %loss))
}

; 8. AGI Adaptive Learning Rate
; Description: Adjusts a learning rate based on a current error and a decay factor.
; Formula: New_LR = Old_LR * (1 - error_magnitude * decay_factor)
; Parameters: %old_learning_rate, %error_magnitude (0-1), %decay_factor (0-1)
; Returns: New learning rate (0-1) or _ERROR_.
alias agi_adaptive_learning_rate {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :LRandErrorAndDecayMustBeNumbers }
    if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0 || $3 > 1) { return %ERR_RETURN $+ :ParametersOutOfRange }

    var %new_lr = $calc($1 * (1 - $2 * $3))
    return $calc($min(1, $max(0, %new_lr))) ; Clamp result to [0,1]
}

; 9. AGI Self-Replicating Pattern Condition (Cellular Automata)
; Description: Checks if a simple 3-cell pattern is capable of "replication" based on a rule.
;              Conceptual: Checks if applying the rule to the pattern's left/center/right would result in the center replicating to the right.
; Parameters: %left_state, %center_state, %right_state, %rule_number (for emergent_pattern_cell)
; Returns: 1 if self-replication condition met, 0 otherwise, or _ERROR_.
alias agi_self_replication_condition {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :StatesAndRuleMustBeNumbers }
    if (($1 != 0 && $1 != 1) || ($2 != 0 && $2 != 1) || ($3 != 0 && $3 != 1)) { return %ERR_RETURN $+ :StatesMustBeBinary }
    if (!($is_non_negative_integer($4)) || $4 > 255) { return %ERR_RETURN $+ :RuleMustBeInteger0To255 }

    ; Simplest interpretation: does the next state of 'center' match current 'center'?
    ; Requires `emergent_pattern_cell` from this or `emergent_systems.mrc`
    var %next_center_state = $emergent_pattern_cell($1, $2, $3, $4)
    if (%next_center_state == %ERR_RETURN) { return %ERR_RETURN } ; Propagate error

    return $iif(%next_center_state == $2, 1, 0) ; Check if center state persists after one step
}

; 10. AGI Global Optimization Step (Conceptual Simulated Annealing)
; Description: Simulates a single step in a global optimization algorithm like simulated annealing.
;              Accepts a new "solution" probabilistically based on its "energy" and a "temperature."
; Parameters: %current_energy, %new_energy, %temperature (for acceptance criteria)
; Returns: "ACCEPT" or "REJECT" or _ERROR_.
alias agi_global_optimization_step {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :EnergiesAndTempMustBeNumbers }
    if (!($is_positive_num($3))) { return %ERR_RETURN $+ :TemperatureMustBePositive }

    if ($2 < $1) {
        return "ACCEPT" ; Always accept better solutions
    } else {
        ; Calculate acceptance probability for worse solutions
        ; P(accept) = exp(-(new_E - current_E) / Temperature)
        var %delta_E = $calc($2 - $1)
        var %acceptance_prob = $calc(exp(0 - (%delta_E / $3)))
        var %random_roll = $rand(1000) / 1000.0

        return $iif(%random_roll <= %acceptance_prob, "ACCEPT", "REJECT")
    }
}
