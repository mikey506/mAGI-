; quantum_bits.mrc for mIRC 5.91 - Balanced Trinary System
; Updated per mAGI research: https://github.com/mikey506/mAGI-/blob/main/research/trinarylogic.txt
; Requires quantum.dll (C++ serialization engine with trinary support)

; === BALANCED TERNARY SYSTEM ===
alias ternary {
  if ($1 == on) { 
    set %ternary_mode 1 
    echo -e * Quantum: Balanced ternary mode activated (5 trits/byte)
  }
  elseif ($1 == off) { 
    unset %ternary_mode 
    echo -e * Quantum: Balanced ternary mode deactivated
  }
  else { return %ternary_mode }
}

; Converts integer to balanced trits (N=-1, 0=0, P=1)
alias int_to_trits {
  var %n = $1, %trits
  while (%n != 0) {
    var %r = $calc(%n % 3)
    if (%r == 2) { 
      %trits = N $+ %trits
      %n = $calc((%n + 1) / 3)
    }
    elseif (%r == -2) {
      %trits = P $+ %trits
      %n = $calc((%n - 1) / 3)
    }
    else {
      %trits = $iif(%r < 0, P, $iif(%r > 0, N, 0)) $+ %trits
      %n = $floor($calc(%n / 3))
    }
  }
  return $left($str(0,5) $+ %trits, -5)
}

; Converts balanced trits to integer
alias trits_to_int {
  var %trits = $1, %n = 0, %len = $len(%trits)
  for (%i = 1; %i <= %len; %i++) {
    %n = %n * 3
    var %char = $mid(%trits, %i, 1)
    if (%char == P) { %n = %n + 1 }
    elseif (%char == N) { %n = %n - 1 }
  }
  return %n
}

; === TRINARY CORE OPERATIONS ===
alias qb_trit_get {
  if (!$1 || !$2) { return 0 }
  var %trit_pos = $2
  var %byte_pos = $ceil($calc(%trit_pos / 5))
  var %offset = $calc(5 - ((%trit_pos - 1) % 5))  ; 1-5 (LSB to MSB)
  
  var %byte = $bvar($1, %byte_pos)
  var %trit_value = $floor($calc(%byte / $pow(3, %offset - 1))) % 3
  
  return $iif(%trit_value == 0, 0, $iif(%trit_value == 1, 1, -1))
}

alias qb_trit_set {
  if (!$1 || !$2 || !$3) { return 0 }
  var %trit_pos = $2, %value = $3
  var %byte_pos = $ceil($calc(%trit_pos / 5))
  var %offset = $calc(5 - ((%trit_pos - 1) % 5))  ; 1-5 (LSB to MSB)
  
  ; Convert -1/0/1 to 0/1/2 storage values
  var %storage = $iif(%value == -1, 2, $iif(%value == 0, 0, 1))
  
  ; Read current byte and modify trit
  var %current = $bvar($1, %byte_pos)
  var %power = $pow(3, %offset - 1)
  var %current_trit = $floor($calc(%current / %power)) % 3
  var %new_byte = %current + $calc((%storage - %current_trit) * %power)
  
  bset &$1 %byte_pos %new_byte
  return 1
}

; === TRINARY ENCODING/DECODING ===
alias qb_ternary_encode {
  var %i = 1, %len = $bvar($1, 0), %result
  while (%i <= %len) {
    var %byte = $bvar($1, %i), %trit_str
    for (%j = 4; %j >= 0; %j--) {
      var %trit = $floor($calc(%byte / $pow(3, %j))) % 3
      %trit_str = %trit_str $+ $iif(%trit == 0, 0, $iif(%trit == 1, P, N))
      %byte = %byte % $pow(3, %j)
    }
    %result = %result $+ %trit_str
    inc %i
  }
  return %result
}

alias qb_ternary_decode {
  if (!$1 || !$2) { echo -e * Quantum Error: Missing parameters | halt }
  var %str = $remove($1, $chr(32)), %len = $len(%str), %binvar = $2
  if ($calc(%len % 5)) { 
    echo -e * Quantum Error: Ternary string length must be multiple of 5
    halt
  }
  bunset & %binvar
  var %i = 1, %byte_pos = 1
  while (%i <= %len) {
    var %trits = $mid(%str, %i, 5), %byte = 0
    for (%j = 0; %j < 5; %j++) {
      var %char = $mid(%trits, $calc(%j + 1), 1)
      var %value = $iif(%char == P, 1, $iif(%char == N, 2, 0))
      %byte = %byte * 3 + %value
    }
    bset -t & %binvar %byte_pos %byte
    inc %byte_pos
    inc %i 5
  }
}

; === OPTIMIZED PACKING/UNPACKING ===
alias qb_pack_trit {
  if ($2 == -1) { var %pos $calc($bvar($1, 0) * 5 + 1) } else { %pos = $2 }
  qb_trit_set $1 %pos $3
  return %pos + 1
}

alias qb_unpack_trit {
  return qb_trit_get $1 $2
}

alias qb_pack_int8 {
  if ($2 == -1) { var %pos $calc($bvar($1, 0) * 5 + 1) } else { %pos = $2 }
  var %trits = $int_to_trits($3)
  qb_trit_set $1 %pos $iif($mid(%trits,1,1) == P, 1, $iif($mid(%trits,1,1) == N, -1, 0))
  qb_trit_set $1 $calc(%pos + 1) $iif($mid(%trits,2,1) == P, 1, $iif($mid(%trits,2,1) == N, -1, 0))
  qb_trit_set $1 $calc(%pos + 2) $iif($mid(%trits,3,1) == P, 1, $iif($mid(%trits,3,1) == N, -1, 0))
  qb_trit_set $1 $calc(%pos + 3) $iif($mid(%trits,4,1) == P, 1, $iif($mid(%trits,4,1) == N, -1, 0))
  qb_trit_set $1 $calc(%pos + 4) $iif($mid(%trits,5,1) == P, 1, $iif($mid(%trits,5,1) == N, -1, 0))
  return %pos + 5
}

alias qb_unpack_int8 {
  var %trit1 = qb_trit_get $1 $2
  var %trit2 = qb_trit_get $1 $calc($2 + 1)
  var %trit3 = qb_trit_get $1 $calc($2 + 2)
  var %trit4 = qb_trit_get $1 $calc($2 + 3)
  var %trit5 = qb_trit_get $1 $calc($2 + 4)
  var %trits = $iif(%trit1 == -1, N, $iif(%trit1 == 1, P, 0)) 
              $+ $iif(%trit2 == -1, N, $iif(%trit2 == 1, P, 0))
              $+ $iif(%trit3 == -1, N, $iif(%trit3 == 1, P, 0))
              $+ $iif(%trit4 == -1, N, $iif(%trit4 == 1, P, 0))
              $+ $iif(%trit5 == -1, N, $iif(%trit5 == 1, P, 0))
  return $trits_to_int(%trits)
}

; === ADVANCED OPERATIONS ===
alias qb_trit_count {
  return $calc($bvar($1, 0) * 5)
}

alias qb_trit_compare {
  if ($qb_trit_count($1) != $qb_trit_count($2)) { return 0 }
  var %i = 1, %total = $qb_trit_count($1)
  while (%i <= %total) {
    if ($qb_trit_get($1, %i) != $qb_trit_get($2, %i)) { return 0 }
    inc %i
  }
  return 1
}

alias qb_trit_slice {
  var %start = $3, %end = $4, %len = $calc(%end - %start + 1)
  bunset &$1
  var %src_byte_start = $ceil($calc(%start / 5))
  var %src_trit_offset = $calc((%start - 1) % 5)
  var %pos = 1
  
  while (%len > 0) {
    var %trits_to_copy = $iif(%len > 5, 5, %len)
    var %value = 0
    
    for (%i = 0; %i < %trits_to_copy; %i++) {
      var %trit_pos = $calc(%start + %i)
      var %trit_val = qb_trit_get $2 %trit_pos
      %value = %value * 3 + $iif(%trit_val == -1, 2, %trit_val)
    }
    
    ; Pad remaining trits with zeros
    while ($calc(%trits_to_copy + %i) < 5) {
      %value = %value * 3
      inc %i
    }
    
    bset &$1 %pos %value
    inc %pos
    %start = %start + %trits_to_copy
    %len = %len - %trits_to_copy
  }
}

; === QUANTUM GATES (TRIT OPERATIONS) ===
alias qb_trit_not {
  return $iif($1 == 1, -1, $iif($1 == -1, 1, 0))
}

alias qb_trit_consensus {
  return $iif($1 == $2 && $2 == $3, $1, 0)
}

alias qb_trit_any {
  return $iif($1 || $2 || $3, $iif($1 == 1 || $2 == 1 || $3 == 1, 1, -1), 0)
}

; === INITIALIZATION ===
alias qb_init {
  if !$dll(quantum.dll, version) { 
    echo -e * Quantum Error: DLL not found! | halt
  }
  echo -e * Quantum Bits initialized (Trinary mode: $iif($ternary,ON,OFF))
}
qb_init
