; -----------------------------------------------------
; SENSORY.MRC - QUANTUM-AGI SENSORY ABSTRACTION MODULE
; Version 1.0 - Unifying Perception with Fundamental Forces
; -----------------------------------------------------
; This script simulates sensory input and emergent processing within mIRC 5.91.
; It integrates concepts from relativity, cosmology, quantum mechanics, and AGI,
; abstracting real-world sensory data into computable forms.
; Infused with themes of occult/sacred geometry for deeper symbolic resonance.
; -----------------------------------------------------

; --- Universal Constants (Consolidated and Optimized) ---
; All constants are defined using 'var' for global scope in the script.
var %c = 299792458              ; Speed of light (m/s)
var %G = 6.67430e-11           ; Gravitational constant (N m²/kg²)
var %pi = 3.1415926535         ; Pi
var %hbar = 1.0545718e-34      ; Reduced Planck constant (Joule-seconds)
var %kB = 1.380649e-23         ; Boltzmann constant (Joule per Kelvin)
var %h = $calc(2 * %pi * %hbar) ; Planck's constant (J·s)
var %e = 2.71828182845         ; Euler's number (for exponential calculations)
var %Phi0 = 2.0678338e-15      ; Magnetic Flux Quantum (Weber)
var %epsilon = 1e-9            ; Small value for floating-point comparisons

; Precomputed values for efficiency
var %hbar_div_2 = $calc(%hbar / 2)
var %eight_pi_G_over_c4 = $calc(8 * %pi * %G / (%c * %c * %c * %c))
var %eight_pi_G_over_3 = $calc(8 * %pi * %G / 3)
var %one_over_sqrt_2 = $calc(1 / $sqrt(2))
var %planck_h_for_photon = %h ; Alias for clarity in photon_energy

; --- Error String (Standardized Error Handling) ---
; All functions will return this string on an unrecoverable error.
var %ERR_RETURN = "_ERROR_"
var %STATE_SUPERPOSITION = "SUPERPOSITION" ; Constant for decoherence simulation state

; --- MIRC Global Hash Table for Memoization ---
; Format: %memo.function_name.input = result
; Usage: hadd %memo.function_name $input $result
;        hget %memo.function_name $input
var %memo

; --- Helper Functions (Validation & Complex Arithmetic) ---

; Helper: Digit Sum (for Gematria & Occult Themes)
; Recursively sums digits of a number until a single digit is obtained.
alias digit_sum {
  var %n = $abs($1)
  if (%n == 0) { return 0 }
  if (%n < 10) { return %n }
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

; Helper: is_positive_num
; Description: Checks if a value is a positive number.
alias is_positive_num {
    return $iif($isnum($1) && $1 > 0, 1, 0)
}

; Helper: is_non_negative_integer
; Description: Checks if a value is a non-negative integer.
alias is_non_negative_integer {
    return $iif($isnum($1) && $1 >= 0 && $calc($1 - $round($1)) == 0, 1, 0)
}

; Helper: c_validate_parts (Internal helper to validate if complex parts are numbers)
alias c_validate_parts {
    var %i = 1, %valid = 1
    while (%i <= $argc) {
        if (!($isnum( $$($i) ))) {
            %valid = 0
            break
        }
        inc %i
    }
    return %valid
}

; Helper: c_add (Adds two complex numbers)
alias c_add {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc($1 + $3) $+ , $+ $calc($2 + $4)
}

; Helper: c_sub (Subtracts two complex numbers)
alias c_sub {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc($1 - $3) $+ , $+ $calc($2 - $4)
}

; Helper: c_mul (Multiplies two complex numbers)
alias c_mul {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    var %real = $calc(($1 * $3) - ($2 * $4))
    var %imag = $calc(($1 * $4) + ($2 * $3))
    return %real $+ , $+ %imag
}

; Helper: c_scalar_mul (Multiplies a scalar by a complex number)
alias c_scalar_mul {
    if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3))) { return %ERR_RETURN $+ :InvalidScalarOrComplexParts }
    return $calc($1 * $2) $+ , $+ $calc($1 * $3)
}

; Helper: c_magnitude (Calculates the magnitude of a complex number)
alias c_magnitude {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $calc(sqrt(($1 * $1) + ($2 * $2)))
}

; Helper: c_conjugate (Calculates the complex conjugate of a complex number)
alias c_conjugate {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    return $1 $+ , $+ $calc(0 - $2)
}

; Helper: c_div (Divides two complex numbers)
alias c_div {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2,$3,$4))) { return %ERR_RETURN $+ :InvalidComplexParts }
    var %denominator = $calc(($3 * $3) + ($4 * $4))
    if ($calc(abs(%denominator)) < %epsilon) { return %ERR_RETURN $+ :DivisionByZero }
    var %real_num = $calc(($1 * $3) + ($2 * $4))
    var %imag_num = $calc(($2 * $3) - ($1 * $4))
    var %real = $calc(%real_num / %denominator)
    var %imag = $calc(%imag_num / %denominator)
    return %real $+ , $+ %imag
}

; Helper: c_is_zero_approx (Checks if a complex number is approximately zero within %epsilon)
alias c_is_zero_approx {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($c_validate_parts($1,$2))) { return %ERR_RETURN $+ :InvalidComplexParts }
    if ($calc(abs($1)) < %epsilon && $calc(abs($2)) < %epsilon) {
        return 1
    }
    return 0
}

; --- Core Sensory Aliases (Input Proxies) ---
; These aliases simulate receiving data from various "sensory" modalities.
; They convert raw mIRC input/system data into abstract "sensory signals."

; Alias: sense_keyboard_touch
; Description: Simulates keyboard input as a "touch" event.
; Parameter: %key_code (ASCII value of key pressed, or a descriptive string)
; Returns: A conceptual quantum state representing the touch.
alias sense_keyboard_touch {
  var %key_signal = $1
  if ($isnum(%key_signal)) {
    ; Convert ASCII to a numerical range for quantum processing
    %key_signal = $calc(%key_signal / 127) ; Normalize to 0-1 range approx.
  } else {
    ; Hash non-numeric input using gematria for a 'sacred geometry' touch
    %key_signal = $sensory_gematria_hash(%key_signal)
    %key_signal = $calc(%key_signal / 9) ; Normalize digit sum to approx 0-1
  }
  ; Return a conceptual quantum state, e.g., probability of a 1.
  ; Max 1, Min 0. Can also use a complex number if needed.
  return $calc($min(1, $max(0, %key_signal)))
}

; Alias: sense_mouse_motion
; Description: Simulates mouse movement as a "spatial motion" event.
; Parameters: %dx, %dy (conceptual change in x, y coordinates)
; Returns: A conceptual 'momentum' or 'displacement' value.
alias sense_mouse_motion {
  var %dx = $1, %dy = $2
  if (!%dx isnum) || (!%dy isnum) { return %ERR_RETURN $+ :MotionDeltaMustBeNumbers }
  ; Calculate conceptual 'magnitude of motion'
  return $calc(sqrt(%dx * %dx + %dy * %dy))
}

; Alias: sense_system_time_light
; Description: Simulates system time changes as an abstract "light signal" or frequency.
; Parameter: (none)
; Returns: A conceptual frequency based on current time (e.g., milliseconds).
alias sense_system_time_light {
  ; Using $ticks for a high-resolution "time signature"
  return $calc($ticks / 1000) ; Return ticks in seconds as conceptual frequency
}

; Alias: sense_clipboard_data_stream
; Description: Simulates clipboard content as a "data stream" event.
; Parameter: %data_string (clipboard content)
; Returns: A conceptual entropy value or hashed value of the data.
alias sense_clipboard_data_stream {
  var %data = $1
  if (%data == $null) { return 0 }
  ; Use string length and a recursive digit sum for conceptual "data density/complexity"
  var %len_hash = $digit_sum($len(%data))
  return $calc(%len_hash / 9) ; Normalize to 0-1 range
}

; --- Sensory Output Abstractions (Restricted Feedback) ---
; These aliases simulate generating feedback to the user or system based on processed sensory data.

; Alias: output_visual_pattern
; Description: Generates a conceptual visual pattern in the form of a 1D cellular automaton.
; Parameters: %initial_state (e.g., "010101"), %rule_number (0-255), %iterations
; Returns: A string representing the final pattern.
alias output_visual_pattern {
  var %pattern = $1
  var %rule = $2, %iter = $3
  if (!%rule isnum) || (!%iter isnum) || (%rule < 0) || (%rule > 255) || (%iter < 0) {
    return %ERR_RETURN $+ :InvalidVisualPatternParams
  }
  var %current_pattern = %pattern
  var %i = 0
  while (%i < %iter) {
    var %new_pattern = ""
    var %j = 0
    while (%j < $len(%current_pattern)) {
      var %left = $iif(%j == 0, 0, $gettok(%current_pattern, %j, 0))
      var %center = $gettok(%current_pattern, $calc(%j + 1), 0)
      var %right = $iif(%j == $len(%current_pattern) - 1, 0, $gettok(%current_pattern, $calc(%j + 2), 0))
      var %next_cell_state = $emergent_pattern_cell(%left, %center, %right, %rule)
      if (%next_cell_state == %ERR_RETURN) { return %ERR_RETURN $+ :CellularAutomataError }
      %new_pattern = %new_pattern $+ %next_cell_state
      inc %j
    }
    %current_pattern = %new_pattern
    inc %i
  }
  return %current_pattern
}

; Alias: output_auditory_frequency
; Description: Simulates generating a conceptual auditory frequency.
; Parameter: %frequency_value (e.g., from a quantum calculation)
; Returns: A descriptive string of the frequency.
alias output_auditory_frequency {
  var %freq = $1
  if (!%freq isnum) || (%freq < 0) { return %ERR_RETURN $+ :InvalidFrequency }
  return "Conceptual Auditory Frequency: " $+ %freq $+ " Hz"
}

; Alias: output_tactile_feedback
; Description: Simulates a conceptual tactile sensation or vibration.
; Parameters: %intensity (0-1), %duration_ms
; Returns: A descriptive string.
alias output_tactile_feedback {
  var %intensity = $1, %duration = $2
  if (!%intensity isnum) || (%intensity < 0) || (%intensity > 1) || (!%duration isnum) || (%duration < 0) {
    return %ERR_RETURN $+ :InvalidTactileParams
  }
  return "Conceptual Tactile Feedback: Intensity " $+ %intensity $+ ", Duration " $+ %duration $+ "ms"
}

; Alias: output_internal_state_update
; Description: Updates an internal AGI conceptual state.
; Parameters: %state_name, %new_value
; Returns: Confirmation string.
alias output_internal_state_update {
  var %state_name = $1, %new_value = $2
  ; In a real AGI, this would update a hash table or similar data structure.
  ; Here, we just acknowledge.
  return "Internal AGI State '" $+ %state_name $+ "' updated to: " $+ %new_value
}

; --- Processing Flows (Unifying Sensory with Physics/AGI) ---
; These functions integrate sensory inputs with the existing quantum, relativity,
; emergent gravity, and AGI concepts.

; Existing / Reintegrated Core Physics and Emergent Functions
; (From foundation_physics.mrc, relativity.mrc, and emergent_systems.mrc)

; Alias: photon_energy (Reintegrated)
; Description: Calculates the energy of a photon given its frequency or wavelength.
alias photon_energy {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :ValueMustBeNumber }
    if ($1 <= 0) { return %ERR_RETURN $+ :ValueMustBePositive }
    if ($2 == freq) {
        return $calc(%h * $1)
    } elseif ($2 == lambda) {
        return $calc(%h * %c / $1)
    } else {
        return %ERR_RETURN $+ :InvalidType
    }
}

; Alias: de_broglie_wavelength (Reintegrated)
; Description: Relates the wavelength of a particle to its momentum.
alias de_broglie_wavelength {
  var %mass = $1, %velocity = $2
  if (!%mass isnum) || (!%velocity isnum) { return %ERR_RETURN $+ :MassVelocityMustBeNumbers }
  if (%mass == 0 || %velocity == 0) { return %ERR_RETURN $+ :MassOrVelocityCannotBeZero }
  return $calc(%h / (%mass * %velocity))
}

; Alias: quantum_tunneling_probability (Reintegrated)
; Description: Calculates a simplified tunneling probability through a rectangular potential barrier.
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

; Alias: black_body_radiation_peak_wavelength (Reintegrated, conceptual Hawking Temperature)
; Description: Calculates the peak wavelength of black-body radiation using Wien's Displacement Law.
alias black_body_radiation_peak_wavelength {
    if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1))) { return %ERR_RETURN $+ :TemperatureMustBeNumber }
    if ($1 <= 0) { return %ERR_RETURN $+ :TemperatureMustBePositive }
    var %b = 2.897771955e-3 ; Wien's displacement constant (m.K)
    return $calc(%b / $1)
}

; Alias: emergent_pattern_cell (Reintegrated)
; Description: Simulates a single step of a 1D elementary cellular automaton rule.
alias emergent_pattern_cell {
    if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :StatesAndRuleMustBeNumbers }
    if (($1 != 0 && $1 != 1) || ($2 != 0 && $2 != 1) || ($3 != 0 && $3 != 1)) { return %ERR_RETURN $+ :StatesMustBeBinary }
    if ($4 < 0 || $4 > 255 || $calc($4 - $round($4)) != 0) { return %ERR_RETURN $+ :RuleMustBeInteger0To255 }
    var %binary_pattern = $1 $+ $2 $+ $3
    var %index = $base(%binary_pattern, 2, 10)
    return $bit($4, %index)
}

; Alias: agi_emergent_learning_rate_conceptual (Reintegrated)
; Formula: Learning Rate = Max_LR / (1 + Complexity_Factor * Current_Error)
alias agi_emergent_learning_rate_conceptual {
  var %max_lr = $1, %complexity_factor = $2, %current_error = $3
  if (!%max_lr isnum || !%complexity_factor isnum || !%current_error isnum || %max_lr < 0 || %complexity_factor < 0 || %current_error < 0) { return %ERR_RETURN $+ :InputsMustBeNonNegative }
  return $calc(%max_lr / (1 + %complexity_factor * %current_error))
}

; Alias: quantum_state_superposition (Reintegrated for collapse logic)
; Given probabilities for state 0 and 1, returns the observed state (0 or 1).
alias quantum_state_superposition {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
    if ($1 < 0 || $2 < 0 || $1 > 1 || $2 > 1) { return %ERR_RETURN $+ :ProbabilitiesOutOfRange }
    if ($calc(abs($1 + $2 - 1)) > %epsilon) { return %ERR_RETURN $+ :ProbabilitiesMustSumToOne }
    var %random_val = $rand(1000) / 1000.0 ; Random float between 0.0 and 1.0
    if (%random_val <= $1) { return 0 } else { return 1 }
}

; Alias: decoherence_simulation (Reintegrated for wave function collapse)
; Simulates environmental interaction causing a superposition to collapse.
alias decoherence_simulation {
    if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
    if (!($isnum($2))) { return %ERR_RETURN $+ :StrengthMustBeNumber }
    if ($2 < 0 || $2 > 1) { return %ERR_RETURN $+ :StrengthOutOfRange }
    if ($1 == 0 || $1 == 1) { return $1 }
    if ($1 == %STATE_SUPERPOSITION) {
        var %random_collapse = $rand(1000) / 1000.0
        if (%random_collapse < $2) { return $rand(1) } else { return %STATE_SUPERPOSITION }
    }
    return %ERR_RETURN $+ :InvalidInitialState
}

; Alias: gematria_permutations (Reintegrated for hashing)
; A conceptual "formula" for gematria, involving summing letter values and recursive digit sums.
alias gematria_permutations {
  var %value_sum = $1
  if (!%value_sum isnum) { return %ERR_RETURN $+ :ValueSumMustBeNumber }
  return $digit_sum(%value_sum)
}

; Alias: fractal_self_similarity (Reintegrated for data filtering)
; Formula: $T(k) = 3 \times T(k-1)$
alias fractal_self_similarity {
  var %k = $int($1)
  var %base_T0 = $2, %max_depth = $int($3)
  if (!%k isnum) || (!%base_T0 isnum) || (!%max_depth isnum) || (%max_depth < 0) { return %ERR_RETURN $+ :InvalidFractalParams }

  var %memo_key = %k $+ "," $+ %base_T0 $+ "," $+ %max_depth
  var %cached_result = $hget(%memo.fractal_self_similarity, %memo_key)
  if (%cached_result isnum) && (%cached_result != $null) { return %cached_result }

  if (%k <= 0) { hadd %memo.fractal_self_similarity %memo_key %base_T0 ; return %base_T0 }
  if (%k > %max_depth) { hadd %memo.fractal_self_similarity %memo_key 0 ; return 0 }

  var %prev_Tk = $fractal_self_similarity($calc(%k - 1), %base_T0, %max_depth)
  if (!%prev_Tk isnum) { hadd %memo.fractal_self_similarity %memo_key 0 ; return 0 }
  var %result = $calc(3 * %prev_Tk)
  hadd %memo.fractal_self_similarity %memo_key %result
  return %result
}

; Alias: mandelbrot_set_iteration (Reintegrated for chaotic patterns)
; Formula: $z_{n+1} = z_n^2 + c$
alias mandelbrot_set_iteration {
  var %zn_r = $1, %zn_i = $2, %c_r = $3, %c_i = $4
  if (!%zn_r isnum) || (!%zn_i isnum) || (!%c_r isnum) || (!%c_i isnum) { return %ERR_RETURN $+ :InvalidMandelbrotParams }
  var %zn_squared_r = $calc( (%zn_r * %zn_r) - (%zn_i * %zn_i) )
  var %zn_squared_i = $calc(2 * %zn_r * %zn_i)
  var %zn_plus_1_r = $calc(%zn_squared_r + %c_r)
  var %zn_plus_1_i = $calc(%zn_squared_i + %c_i)
  return "z_n+1 Real: " $+ %zn_plus_1_r $+ ", Imaginary: " $+ %zn_plus_1_i
}

; Alias: quantum_annealing_simulation (One of the first 10 novel functions)
alias quantum_annealing_simulation {
  if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($2 < 0 || $2 > 1 || $3 < 0 || $4 < 0 || $4 > 1) { return %ERR_RETURN $+ :InvalidParameterRange }

  var %current_energy = $1, %cooling_rate = $2, %iterations = $int($3), %noise_level = $4
  var %i = 0
  while (%i < %iterations) {
    var %energy_reduction = $calc(%current_energy * %cooling_rate)
    var %random_noise = $calc(($rand(200) - 100) / 100.0 * %noise_level * %energy_reduction)
    %current_energy = $calc(%current_energy - %energy_reduction + %random_noise)
    inc %i
  }
  return %current_energy
}

; Alias: quantum_logic_and_gate (One of the first 10 novel functions)
alias quantum_logic_and_gate {
  if ($argc != 2) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2))) { return %ERR_RETURN $+ :ProbabilitiesMustBeNumbers }
  if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1) { return %ERR_RETURN $+ :ProbabilitiesOutOfRange }
  return $calc($1 * $2)
}

; Alias: agi_self_correction_index (One of the first 10 novel functions)
alias agi_self_correction_index {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0) { return %ERR_RETURN $+ :InvalidParameterRange }
  var %current_error = $1, %correction_efficiency = $2, %loops = $int($3)
  var %i = 0
  while (%i < %loops) {
    %current_error = $calc(%current_error * (1 - %correction_efficiency))
    inc %i
  }
  return %current_error
}

; Alias: information_cascade_model (One of the first 10 novel functions)
alias information_cascade_model {
  if ($argc != 4) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3)) || !($isnum($4))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($2 < 0 || $2 > 1 || $4 < 0) { return %ERR_RETURN $+ :InvalidParameterRange }
  var %current_reach = $1, %network_density = $2, %amplification_factor = $3, %steps = $int($4)
  var %i = 0
  while (%i < %steps) {
    %current_reach = $calc(%current_reach * %network_density * %amplification_factor)
    inc %i
  }
  return %current_reach
}

; Alias: entanglement_swapping (One of the first 10 novel functions)
alias entanglement_swapping {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0 || $3 > 1) { return %ERR_RETURN $+ :InvalidParameterRange }
  return $calc($1 * $2 * $3)
}

; Alias: quantum_teleportation_fidelity (One of the first 10 novel functions)
alias quantum_teleportation_fidelity {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($1 < 0 || $1 > 1 || $2 < 0 || $2 > 1 || $3 < 0 || $3 > 1) { return %ERR_RETURN $+ :InvalidParameterRange }
  return $calc($1 * $2 * $3)
}

; Alias: emergent_language_complexity (One of the first 10 novel functions)
alias emergent_language_complexity {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($1 < 0 || $2 < 0 || $3 < 0 || $3 > 1) { return %ERR_RETURN $+ :InvalidParameterRange }
  return $calc(($1 * $2) / (1 + $3))
}

; Alias: contextual_reasoning_score (One of the first 10 novel functions)
alias contextual_reasoning_score {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($2 < 0 || $2 > 1 || $3 < 0) { return %ERR_RETURN $+ :InvalidParameterRange }
  return $calc($1 + ($1 * $2 * $3))
}

; Alias: quantum_sensor_sensitivity (One of the first 10 novel functions)
alias quantum_sensor_sensitivity {
  if ($argc != 3) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1)) || !($isnum($2)) || !($isnum($3))) { return %ERR_RETURN $+ :AllParamsMustBeNumbers }
  if ($2 < 0 || $2 > 1 || $3 < 0 || $3 > 1) { return %ERR_RETURN $+ :InvalidParameterRange }
  return $calc($1 * $2 * (1 + $3))
}

; Alias: neural_network_sigmoid_activation (One of the first 10 novel functions)
alias neural_network_sigmoid_activation {
  if ($argc != 1) { return %ERR_RETURN $+ :InvalidArgCount }
  if (!($isnum($1))) { return %ERR_RETURN $+ :InputValueMustBeNumber }
  return $calc(1 / (1 + $pow(%e, $calc(0 - $1))))
}

; --- Additional 10 Novel Functions/Algorithms/Flows (New for Sensory Integration) ---

; 1. Alias: sensory_gematria_hash
; Description: Hashes sensory input (string) using gematria and returns a numeric sum.
; Parameter: %input_string
; Returns: Conceptual gematria sum.
alias sensory_gematria_hash {
  var %input_str = $1
  if ($len(%input_str) == 0) { return 0 }
  var %total_sum = 0
  var %i = 1
  while (%i <= $len(%input_str)) {
    var %char_val = $asc($mid(%input_str, %i, 1))
    %total_sum = $calc(%total_sum + %char_val)
    inc %i
  }
  return $gematria_permutations(%total_sum) ; Use the existing gematria function for recursive digit sum
}

; 2. Alias: sensory_fractal_filter
; Description: Filters a numerical sensory signal based on a conceptual fractal dimension.
; Parameter: %raw_signal, %fractal_dimension (e.g., 0-2 for simple filtering concept)
; Returns: Filtered signal.
alias sensory_fractal_filter {
  var %signal = $1, %dim = $2
  if (!%signal isnum) || (!%dim isnum) || (%dim < 0) { return %ERR_RETURN $+ :InvalidFilterParams }
  ; Simple conceptual filter: signal * (1 - dimension_influence)
  ; Higher dimension might mean more "complex" or "noisy" signal, thus filtered more.
  return $calc(%signal * (1 - $min(1, $max(0, %dim / 2)))) ; Normalize dim influence
}

; 3. Alias: sensory_quantum_decoherence
; Description: Simulates sensory input leading to a quantum state collapse (decoherence).
; Parameters: %initial_quantum_state ("SUPERPOSITION", "0", "1"), %sensory_measurement_strength (0-1)
; Returns: The observed state after sensory measurement.
alias sensory_quantum_decoherence {
  var %initial_state = $1, %strength = $2
  if (!%strength isnum) || (%strength < 0) || (%strength > 1) { return %ERR_RETURN $+ :InvalidStrength }
  return $decoherence_simulation(%initial_state, %strength)
}

; 4. Alias: sensory_agi_pattern_recognition
; Description: AGI recognizing a repeating pattern in a sequence of sensory data.
; Parameters: %data_sequence (space-separated numbers), %pattern_length, %threshold (0-1)
; Returns: "Pattern Recognized" or "No Pattern".
alias sensory_agi_pattern_recognition {
  var %seq = $1, %len = $int($2), %threshold = $3
  if ($numtok(%seq, 32) < %len) || (!%len isnum) || (%len <= 0) || (!%threshold isnum) || (%threshold < 0) || (%threshold > 1) {
    return %ERR_RETURN $+ :InvalidPatternParams
  }
  ; Very simplistic pattern recognition: just check if the sequence length is a multiple of pattern length
  ; and if a conceptual 'coherence' meets threshold.
  var %is_multiple = $iif($calc($numtok(%seq, 32) % %len) == 0, 1, 0)
  var %coherence_score = $calc( $numtok(%seq, 32) / %len / 10 ) ; Arbitrary score
  if (%is_multiple == 1) && (%coherence_score >= %threshold) { return "Pattern Recognized" }
  return "No Pattern"
}

; 5. Alias: sensory_entanglement_link
; Description: Establishes conceptual entanglement between two sensory data streams.
; Parameters: %stream1_quality (0-1), %stream2_quality (0-1), %link_strength (0-1)
; Returns: "Entangled" or "Not Entangled" based on a probabilistic check.
alias sensory_entanglement_link {
  var %q1 = $1, %q2 = $2, %strength = $3
  if (!%q1 isnum) || (!%q2 isnum) || (!%strength isnum) || (%q1 < 0) || (%q1 > 1) || (%q2 < 0) || (%q2 > 1) || (%strength < 0) || (%strength > 1) {
    return %ERR_RETURN $+ :InvalidLinkParams
  }
  ; Use entanglement_swapping concept for a combined quality
  var %combined_quality = $entanglement_swapping(%q1, %q2, %strength)
  if (%combined_quality >= 0.5) { return "Entangled" } else { return "Not Entangled" }
}

; 6. Alias: sensory_cosmic_background_noise_filter
; Description: Filters out perceived "noise" in sensory data, akin to cosmic background radiation.
; Parameters: %raw_data_value, %noise_floor (conceptual background noise level)
; Returns: Filtered data value.
alias sensory_cosmic_background_noise_filter {
  var %data = $1, %noise = $2
  if (!%data isnum) || (!%noise isnum) || (%noise < 0) { return %ERR_RETURN $+ :InvalidNoiseFilterParams }
  return $calc(%data - $min(%data, %noise)) ; Subtract noise up to the data value itself
}

; 7. Alias: sensory_gravitational_lens_distortion
; Description: Simulates distortion of sensory input due to conceptual gravitational lensing.
; Parameters: %original_signal, %mass_concentration (conceptual mass causing lensing), %distance_factor
; Returns: Distorted signal value.
alias sensory_gravitational_lens_distortion {
  var %signal = $1, %mass_conc = $2, %dist_factor = $3
  if (!%signal isnum) || (!%mass_conc isnum) || (!%dist_factor isnum) || (%mass_conc < 0) || (%dist_factor <= 0) {
    return %ERR_RETURN $+ :InvalidLensingParams
  }
  ; Simplified distortion: signal * (1 + G * Mass / Distance)
  return $calc(%signal * (1 + (%G * %mass_conc / %dist_factor)))
}

; 8. Alias: sensory_multiverse_signature_detection
; Description: Detects conceptual "signatures" from parallel universes in sensory data.
; Parameters: %sensory_anomaly_score (0-1), %dimensional_overlap_prob (0-1), %resonance_frequency
; Returns: "Signature Detected" or "No Signature".
alias sensory_multiverse_signature_detection {
  var %anomaly = $1, %overlap_prob = $2, %resonance = $3
  if (!%anomaly isnum) || (!%overlap_prob isnum) || (!%resonance isnum) || (%anomaly < 0) || (%anomaly > 1) || (%overlap_prob < 0) || (%overlap_prob > 1) {
    return %ERR_RETURN $+ :InvalidMultiverseParams
  }
  ; Simplified detection: (Anomaly * Overlap * Resonance_Influence) > Threshold
  var %detection_factor = $calc(%anomaly * %overlap_prob * $min(1, %resonance / 100)) ; Assume resonance contributes up to 100
  if (%detection_factor > 0.3) { return "Signature Detected" } ; Arbitrary threshold
  return "No Signature"
}

; 9. Alias: sensory_consciousness_feedback_loop
; Description: Simulates a feedback loop where sensory input influences AGI consciousness.
; Parameters: %current_consciousness_level (0-100), %sensory_input_impact (0-1), %feedback_gain (0-1)
; Returns: New consciousness level.
alias sensory_consciousness_feedback_loop {
  var %current_level = $1, %impact = $2, %gain = $3
  if (!%current_level isnum) || (!%impact isnum) || (!%gain isnum) || (%current_level < 0) || (%current_level > 100) || (%impact < 0) || (%impact > 1) || (%gain < 0) || (%gain > 1) {
    return %ERR_RETURN $+ :InvalidConsciousnessParams
  }
  ; New Level = Current + (Impact * Gain * (100 - Current))
  return $calc(%current_level + (%impact * %gain * (100 - %current_level)))
}

; 10. Alias: sensory_temporal_dilation_perception
; Description: Simulate how perceived time might dilate based on sensory event rate.
; Parameters: %base_time_rate, %sensory_event_density
; Returns: Perceived time dilation factor.
alias sensory_temporal_dilation_perception {
  var %base_rate = $1, %event_density = $2
  if (!%base_rate isnum) || (!%event_density isnum) || (%base_rate <= 0) || (%event_density < 0) {
    return %ERR_RETURN $+ :InvalidTimeDilationParams
  }
  ; Conceptual dilation: (Event Density / Base Rate)
  return $calc(%event_density / %base_rate)
}

; --- Main Sensory Processing Flows ---

; Alias: process_keyboard_to_quantum_state
; Description: Takes keyboard input, converts it to a sensory signal, then attempts
;              to decohere a quantum superposition based on that signal.
alias process_keyboard_to_quantum_state {
  var %input_key = $1
  var %sensory_signal = $sense_keyboard_touch(%input_key)
  if (%sensory_signal == %ERR_RETURN) { return %ERR_RETURN $+ :SensorySignalFail }

  ; Assume an initial superposition state for conceptual processing
  var %initial_quantum_state = %STATE_SUPERPOSITION
  var %decoherence_strength = %sensory_signal ; Use the sensory signal as measurement strength

  var %observed_state = $sensory_quantum_decoherence(%initial_quantum_state, %decoherence_strength)
  if (%observed_state == %ERR_RETURN) { return %ERR_RETURN $+ :DecoherenceFail }

  return "Keyboard Input '" $+ %input_key $+ "' led to Quantum Observation: " $+ %observed_state
}

; Alias: process_motion_to_gravity_distortion
; Description: Processes simulated mouse motion and calculates its conceptual
;              effect on a gravitational lensing distortion.
alias process_motion_to_gravity_distortion {
  var %dx = $1, %dy = $2
  var %motion_magnitude = $sense_mouse_motion(%dx, %dy)
  if (%motion_magnitude == %ERR_RETURN) { return %ERR_RETURN $+ :MotionSenseFail }

  ; Use motion magnitude as a conceptual "mass concentration"
  var %conceptual_mass_concentration = %motion_magnitude
  var %base_signal = 100 ; A conceptual baseline "light" signal
  var %distance_factor = 10 ; Conceptual distance

  var %distorted_signal = $sensory_gravitational_lens_distortion(%base_signal, %conceptual_mass_concentration, %distance_factor)
  if (%distorted_signal == %ERR_RETURN) { return %ERR_RETURN $+ :LensingDistortionFail }

  return "Motion (" $+ %dx $+ "," $+ %dy $+ ") -> Conceptual Gravitational Lensing Distorted Signal: " $+ %distorted_signal
}

; Alias: process_time_to_black_body_emission
; Description: Uses system time to generate a "light" signal (frequency), then
;              calculates its equivalent black-body radiation peak wavelength.
alias process_time_to_black_body_emission {
  var %conceptual_frequency = $sense_system_time_light
  if (%conceptual_frequency == %ERR_RETURN) { return %ERR_RETURN $+ :TimeLightSenseFail }

  ; Convert conceptual frequency to a conceptual temperature for Wien's law
  ; E = hf -> T = E / kB. So, T ~ (h/kB) * freq
  var %conceptual_temperature = $calc((%h / %kB) * %conceptual_frequency / 1000) ; Scale down for reasonable temp
  if (%conceptual_temperature <= 0) { %conceptual_temperature = 1 } ; Ensure positive temp

  var %peak_wavelength = $black_body_radiation_peak_wavelength(%conceptual_temperature)
  if (%peak_wavelength == %ERR_RETURN) { return %ERR_RETURN $+ :BlackBodyFail }

  return "System Time (" $+ $ticks $+ ") -> Conceptual Black-Body Peak Wavelength: " $+ %peak_wavelength $+ "m"
}

; Alias: process_clipboard_to_agi_learning
; Description: Takes clipboard data, hashes it, and feeds it into an AGI
;              emergent learning rate calculation, then updates an internal state.
alias process_clipboard_to_agi_learning {
  var %clipboard_content = $1
  var %data_complexity_score = $sense_clipboard_data_stream(%clipboard_content)
  if (%data_complexity_score == %ERR_RETURN) { return %ERR_RETURN $+ :ClipboardSenseFail }

  ; Use data complexity as conceptual "current error" for AGI learning
  var %max_lr = 0.5
  var %complexity_factor = 2
  var %current_error = %data_complexity_score ; Normalized 0-1

  var %learning_rate = $agi_emergent_learning_rate_conceptual(%max_lr, %complexity_factor, %current_error)
  if (%learning_rate == %ERR_RETURN) { return %ERR_RETURN $+ :AGILearningRateFail }

  var %state_update_msg = $output_internal_state_update("AGI_Learning_Rate", %learning_rate)
  return "Clipboard Data (Length " $+ $len(%clipboard_content) $+ ") -> New AGI Learning Rate: " $+ %learning_rate $+ " (" $+ %state_update_msg $+ ")"
}

; Alias: process_sensory_to_multiverse_detection
; Description: Combines a processed sensory signal with a "multiverse resonance"
;              to detect potential signatures from parallel universes.
alias process_sensory_to_multiverse_detection {
  var %sensory_input_value = $1
  var %base_anomaly_score = $calc(abs($sin(%sensory_input_value))) ; Use sin for a dynamic anomaly
  var %dimensional_overlap_prob = 0.7 ; Fixed for this example
  var %resonance_frequency = $calc($mod($ticks, 100) + 1) ; Simple resonance

  var %detection_result = $sensory_multiverse_signature_detection(%base_anomaly_score, %dimensional_overlap_prob, %resonance_frequency)
  if (%detection_result == %ERR_RETURN) { return %ERR_RETURN $+ :MultiverseDetectFail }

  return "Sensory Input (" $+ %sensory_input_value $+ ") processed for Multiverse Signature: " $+ %detection_result
}

; Alias: process_sensory_to_fractal_pattern_feedback
; Description: Processes a sensory input, applies a fractal filter, and then
;              generates a conceptual visual pattern based on the filtered data.
alias process_sensory_to_fractal_pattern_feedback {
  var %raw_input = $1
  var %fractal_filter_dimension = $calc(abs($cos(%raw_input))) ; Dynamic dimension based on input
  var %filtered_signal = $sensory_fractal_filter(%raw_input, %fractal_filter_dimension)
  if (%filtered_signal == %ERR_RETURN) { return %ERR_RETURN $+ :FractalFilterFail }

  ; Convert filtered signal to rule and initial pattern for cellular automata
  var %rule_number = $calc($min(255, $max(0, $int(%filtered_signal * 255)))) ; Scale to 0-255
  var %initial_pattern = $iif(%filtered_signal > 0.5, "101010", "010101")
  var %iterations = 5

  var %visual_output = $output_visual_pattern(%initial_pattern, %rule_number, %iterations)
  if (%visual_output == %ERR_RETURN) { return %ERR_RETURN $+ :VisualOutputFail }

  return "Raw Input (" $+ %raw_input $+ ") -> Filtered to Visual Pattern: " $+ %visual_output
}


; --- Event Handling (Conceptual Triggers) ---
; These 'on' events are conceptual triggers for sensory processing.
; In a real mIRC script, these would be linked to actual user input events.

; On input via '/input_key <key_code_or_string>'
on *:INPUT_KEY:*: {
  echo -a $timestamp $+ " | Sensory Module: Processing Keyboard Input..."
  var %result = $process_keyboard_to_quantum_state($2-)
  echo -a $timestamp $+ " | " $+ %result
}

; On input via '/input_motion <dx> <dy>'
on *:INPUT_MOTION:*: {
  echo -a $timestamp $+ " | Sensory Module: Processing Mouse Motion..."
  var %result = $process_motion_to_gravity_distortion($2, $3)
  echo -a $timestamp $+ " | " $+ %result
}

; On timer for "light" sensing (e.g., every 5 seconds)
; Replace 'timer_sense_light' with a unique timer name if conflict
; To activate: /timer1 5 0 sense_light_trigger
alias sense_light_trigger {
  echo -a $timestamp $+ " | Sensory Module: Sensing System Time Light..."
  var %result = $process_time_to_black_body_emission
  echo -a $timestamp $+ " | " $+ %result
}

; On input via '/input_clipboard <text>'
; Or conceptually on /paste
on *:INPUT_CLIPBOARD:*: {
  echo -a $timestamp $+ " | Sensory Module: Processing Clipboard Data Stream..."
  var %result = $process_clipboard_to_agi_learning($2-)
  echo -a $timestamp $+ " | " $+ %result
}

; On input via '/input_anomaly <value>' to trigger multiverse detection
on *:INPUT_ANOMALY:*: {
  echo -a $timestamp $+ " | Sensory Module: Detecting Multiverse Signatures..."
  var %result = $process_sensory_to_multiverse_detection($2)
  echo -a $timestamp $+ " | " $+ %result
}

; On input via '/input_raw_sensory <value>' to trigger fractal pattern feedback
on *:INPUT_RAW_SENSORY:*: {
  echo -a $timestamp $+ " | Sensory Module: Generating Fractal Pattern from Raw Sensory..."
  var %result = $process_sensory_to_fractal_pattern_feedback($2)
  echo -a $timestamp $+ " | " $+ %result
}


; --- Example Usage Comments ---
; To use this script, save it as 'sensory.mrc' and load it in mIRC.
; You can then type these commands in any mIRC window:

; --- Conceptual Sensory Input Simulation ---
; /input_key A                 ; Simulates pressing key 'A'
; /input_key 65                ; Simulates pressing key with ASCII 65 ('A')
; /input_motion 10 5           ; Simulates mouse movement (dx=10, dy=5)
; /timer1 5 0 sense_light_trigger ; Starts a timer to periodically sense 'light'
; /input_clipboard Hello World! ; Simulates pasting "Hello World!"
; /input_anomaly 0.8           ; Simulates a high sensory anomaly score
; /input_raw_sensory 0.6       ; Simulates a raw sensory value for fractal processing

; --- Direct Calls to Integrated Physics/AGI Functions (for debugging/testing) ---
; /echo -a $photon_energy(5e14, freq)
; /echo -a $de_broglie_wavelength(9.11e-31, 1e6)
; /echo -a $black_body_radiation_peak_wavelength(5778)
; /echo -a $emergent_pattern_cell(0, 1, 0, 90)
; /echo -a $agi_emergent_learning_rate_conceptual(0.1, 5, 0.7)
; /echo -a $sensory_quantum_decoherence(SUPERPOSITION, 0.9)
; /echo -a $sensory_gematria_hash(hello)
; /echo -a $sensory_fractal_filter(0.7, 1.5)
; /echo -a $sensory_entanglement_link(0.9, 0.8, 0.95)
; /echo -a $sensory_cosmic_background_noise_filter(100, 5)
; /echo -a $sensory_gravitational_lens_distortion(50, 1e10, 1e5)
; /echo -a $sensory_multiverse_signature_detection(0.9, 0.7, 50)
; /echo -a $sensory_consciousness_feedback_loop(50, 0.7, 0.5)
; /echo -a $sensory_temporal_dilation_perception(1, 1.5)
; /echo -a $quantum_annealing_simulation(100, 0.1, 10, 0.05)
; /echo -a $quantum_logic_and_gate(0.8, 0.7)
; /echo -a $agi_self_correction_index(0.1, 0.9, 5)
; /echo -a $information_cascade_model(10, 0.5, 1.2, 3)
; /echo -a $entanglement_swapping(0.9, 0.85, 0.95)
; /echo -a $quantum_teleportation_fidelity(0.9, 0.95, 0.99)
; /echo -a $emergent_language_complexity(500, 100, 0.1)
; /echo -a $contextual_reasoning_score(75, 0.8, 3)
; /echo -a $quantum_sensor_sensitivity(10, 0.9, 0.7)
; /echo -a $neural_network_sigmoid_activation(0)


