; -----------------------------------------------------
; GLOBAL UNIVERSAL CONSTANTS MODULE
; Contains fundamental physical constants used across all physics modules.
; -----------------------------------------------------

; This file should be loaded once before other physics modules.

; Universal Constants (Approximated for MIRC-like scripting)
; All constants are defined using 'var' for global scope in the script, as is standard practice for mIRC scripting.
var $c = 299792458 ; [cite_start]Speed of light (m/s) [cite: 3, 4]
var $G = 6.67430e-11 ; [cite_start]Gravitational constant (N m²/kg²) [cite: 5]
var $pi = 3.1415926535 ; [cite_start]Pi [cite: 6]
var $hbar = 1.0545718e-34 ; [cite_start]Reduced Planck constant (Joule-seconds) [cite: 7]
var $kB = 1.380649e-23 ; [cite_start]Boltzmann constant (Joules per Kelvin) [cite: 8]
var $h = $calc(2 * $pi * $hbar) ; Planck's constant (J·s)
var $e = 2.71828182845 ; Euler's number
var $Phi0 = 2.0678338e-15 ; Magnetic Flux Quantum (Weber)
; Additional constants can be added here.

; Precomputed values for efficiency
[cite_start]var $hbar_over_2 = $calc($hbar / 2) [cite: 9]
[cite_start]var $eight_pi_G_over_c4 = $calc(8 * $pi * $G / pow($c, 4)) [cite: 9]
[cite_start]var $eight_pi_G_over_3 = $calc(8 * $pi * $G / 3) [cite: 9]
[cite_start]var $one_over_sqrt_2 = $calc(1 / $sqrt(2)) [cite: 9]
; [cite_start]More precomputed values can be added based on frequently used expressions. [cite: 10]
