
.CODE
; -------------------------------------------------------------------------------------	;
; void CalculateSpectrum(spectrum_type* Spectrum, signal_type* Signal)					;
;	Прямое преобразование Фурье. Вычисляет спектр Spectrum по сигналу Signal			;
;	Типы данных spectrum_type и signal_type, а так же размер сигнала					;
;	определяются в файле Tuning.h														;
; -------------------------------------------------------------------------------------	;
CalculateSpectrum PROC	; [RCX] - Spectrum
						; [RDX] - Signal

	push R10
	push R11
	push R12
	push R13
	push R14
	push R15
	push RDI
	push RSI

	mov R15, [RDX]		; loading signal into R15

		; R[k] = 0 for k from 8 to 14

	xor R8, R8
	xor R9, R9
	xor R10, R10
	xor R11, R11
	xor R12, R12
	xor R13, R13
	xor R14, R14

		; loading individual signals into registers by shifting the R15 by 8 bits >>

	mov R8b, R15b
	shr r15, 8
	mov R12b, R15b
	shr r15, 8
	mov R9b, R15b
	shr r15, 8
	mov R13b, R15b
	shr r15, 8
	mov R10b, R15b
	shr r15, 8
	mov R14b, R15b
	shr r15, 8
	mov R11b, R15b
	shr r15, 8

		; Registers:
		;	R8-  x[0]
		;	R9-  x[2]
		;	R10- x[4]
		;	R11- x[6]
		;	R12- x[1]
		;	R13- x[3]
		;	R14- x[5]
		;	R15- x[7]
		
	sal r8, 7 * 8
	sar r8, 7 * 8
	sal r9, 7 * 8
	sar r9, 7 * 8
	sal r10, 7 * 8
	sar r10, 7 * 8
	sal r11, 7 * 8
	sar r11, 7 * 8
	sal r12, 7 * 8
	sar r12, 7 * 8
	sal r13, 7 * 8
	sar r13, 7 * 8
	sal r14, 7 * 8
	sar r14, 7 * 8
	sal r15, 7 * 8
	sar r15, 7 * 8

		; processing even indexes
	mov RSI, R10	; loading x[4] into RSI
	add R10, R8		; calculate e[0] and write it to R10
	sub R8, RSI		; calculate e[1] and write it to R8
	mov RSI, R9		; loading x[2] into RSI
	add R9, R11		; calculate o[0] and write it to R9
	sub R11, RSI	; calculate Imaginary part of o[1] and write it to R11

		; Registers :
		;	R8-  e[1]	Has no imaginary part
		;	R9-  o[0]	Has no imaginary part
		;	R10- e[0]	Has no imaginary part
		;	R11- o[1]	Has no real part

		; E[0] = o[0] + e[0], therefore E[0] = Re(E[0])
		; E[2] = o[0] - e[0], therefore E[2] = Re(E[2])

	mov RSI, R9		; loading o[0] into RSI
	add r9, R10		; calculate E[0] and write it to R9
	sub r10, RSI	; calculate E[2] and write it to R10

		; E[1] = Re(e[1]) + Im(o[1]), therefore no calculations are performed
		; E[3] = Re(e[1]) - Im(o[1])

		; swap section where calculated values of Es are getting stored into R12-R15

	mov RSI, R12
	mov R12, R8
	mov R8, RSI

	mov RSI, R13
	mov R13, R9
	mov R9, RSI

	mov RSI, R14
	mov R14, R10
	mov R10, RSI

	mov RSI, R15
	mov R15, R11
	mov R11, RSI

		; Registers:
		;	R12- Re(e[1])	Re(E[1]) and Re(E[3])
		;	R13- E[0]		Only has real part
		;	R14- E[2]		Only has real part
		;	R15- Im(o[1])	Im(E[1]) and -Im(E[3])

		; processing odd indexes
	mov RSI, R10		; loading x[5] into RSI
	add R10, R8			; calculate e[0] and write it to R10
	sub R8, RSI			; calculate e[1] and write it to R8
	mov RSI, R9			; loading x[3] into RSI
	add R9, R11			; calculate o[0] and write it to R9
	sub R11, RSI		; calculate Imaginary part of o[1] and write it to R11
	
		; Registers :
		;	R8-  e[1]	Has no imaginary part
		;	R9-  o[0]	Has no imaginary part
		;	R10- e[0]	Has no imaginary part
		;	R11- o[1]	Has no real part

		; O is an array which contains a vector before it's multiplied by W8

		; O[0] = o[0] + e[0], therefore O[0] = Re(O[0])
		; O[2] = o[0] - e[0], therefore O[2] = Re(O[2])

	mov RSI, R9		; load o[0] to RSI
	add r9, R10		; calculate O[0] and write it to R9
	sub r10, RSI	; calculate O[2] and write it to R10

		; O[1] = Re(e[1]) + Im(o[1]), therefore no calculations are performed
		; O[3] = Re(e[1]) - Im(o[1])
		
	mov RDI, R15

		; Registers :
		;	R8-  Re(O[1]) = Re(O[3])
		;	R9-  O[0]					Has no imaginary part
		;	R10- O[2]					Has no imaginary part
		;	R11- Im(O[1]) = -Im(O[3])
		;	R12- Re(E[1]) = Re(E[3])
		;	R13- E[0]					Has no imaginary part
		;	R14- E[2]					Has no imaginary part
		;	RDI- Im(E[1]) = -Im(E[3])

	finit

		; x[0] = Re(E[0]) + Re(O[0])
	mov R15, R13
	add R15, R9
	mov dword ptr [	RCX			], R15d		; spectrum[0] = x[0]
	fild dword ptr [RCX			]			; convert spectrum[0] from int to float
	fstp dword ptr [RCX			]			; put float version back to spectrum[0]
	mov dword ptr [	RCX + 8 * 4	], 0		; Im(x[0]) = 0

		; x[4] = Re(E[0]) - Re(O[0])
	mov R15, R13
	sub R15, R9
	mov dword ptr [ RCX + 4 * 4], R15d		; spectrum[4] = x[4]
	fild dword ptr [RCX + 4 * 4]			; convert spectrum[4] from int to float
	fstp dword ptr [RCX + 4 * 4]			; put float version back to spectrum[4]
	mov dword ptr [	RCX + (4 + 8) * 4], 0	; Im(x[0]) = 0

		; R9 and R13 can be used for storing values as they don't contain any unprocessed data
		
		; x[2] = Re(E[2]) - Im(O[2])
	mov dword ptr [	RCX +	 2 * 4		], R14d		; Re(spectrum[2]) = E[2]
	fild dword ptr [RCX +	 2 * 4		]
	fstp dword ptr [RCX +	 2 * 4		]
	mov R15, 0
	sub R15, R10
	mov dword ptr [RCX +	(2 + 8) * 4	], R15d		; Im(spectrum[2]) = (-1) * O[2]
	fild dword ptr [RCX +	(2 + 8) * 4	]
	fstp dword ptr [RCX +	(2 + 8) * 4	]
	
		; x[6] = Re(E[2]) + Im(O[2])
	mov dword ptr [RCX +	 6 * 4		], R14d		; Re(spectrum[6]) = E[2]
	fild dword ptr [RCX +	 6 * 4		]
	fstp dword ptr [RCX +	 6 * 4		]
	mov dword ptr [RCX +	(6 + 8) * 4	], R10d		; Im(spectrum[6]) = O[2]
	fild dword ptr [RCX +	(6 + 8) * 4	]
	fstp dword ptr [RCX +	(6 + 8) * 4	]

		; storing sqrt(2) / 2 in ST(0)
	mov dword ptr [RCX + 1 * 4], 2
	fild dword ptr [RCX + 1 * 4]
	fsqrt
	fild dword ptr [RCX + 1 * 4]
	fdivp
	
		; x[1] =
		;	Re(E[1]) [R12] + Im(E[1]) [RDI] +
		;	T(1-i) * (Re(e[1])[R8] + Im(o[1])[R11])
		;	[ R12 + TR8 + TR11 ; RDI + TR11 - TR8 ]
		; TRegister := (sqrt(2) / 2) * Register

		; pushing R12, R8, R11 and RDI onto ST, multiplying R8 and R11 by sqrt(2) / 2
	mov dword ptr [RCX + 1 * 4], R12d
	mov R15, RDI
	mov dword ptr [RCX + 3 * 4], R15d
	mov dword ptr [RCX + 5 * 4], R8d
	mov dword ptr [RCX + 7 * 4], R11d

	fild dword ptr [RCX + 1 * 4]
	fild dword ptr [RCX + 3 * 4]
	fild dword ptr [RCX + 5 * 4]
	fmul ST(0), ST(3)
	fild dword ptr [RCX + 7 * 4]
	fmul ST(0), ST(4)
	fld ST(3)

		; Stack:
		;	ST(0)- Accumulator (R12)
		;	ST(1)- TR11
		;	ST(2)- TR8
		;	ST(3)- R15 (RDI)
		;	ST(4)- R12
		;	ST(5)- sqrt(2) / 2

		; x[1] =
		;	[ R12 + TR8 + TR11 ; RDI + TR11 - TR8 ]
		;	[ ST(4)+ST(2)+ST(1); ST(3)+ST(1)-ST(2)]
	fadd ST(0), ST(2)
	fadd ST(0), ST(1)
	fstp dword ptr [RCX + 1 * 4]

	fld ST(2)
	fadd ST(0), ST(1)
	fsub ST(0), ST(2)
	fstp dword ptr [RCX + (1 + 8) * 4]
		
		; x[5] =
		;	Re(E[1]) [R12] + Im(E[1]) [RDI] -
		;	T(1-i) * (Re(e[1])[R8] + Im(o[1])[R11])
		;	[ R12 - TR8 - TR11 ; RDI - TR11 + TR8 ]
		;	[ ST(4)-ST(2)-ST(1); ST(3)-ST(1)+ST(2)]
	fld ST(3)
	fsub ST(0), ST(2)
	fsub ST(0), ST(1)
	fstp dword ptr [RCX + 5 * 4]
	
	fld ST(2)
	fsub ST(0), ST(1)
	fadd ST(0), ST(2)
	fstp dword ptr [RCX + (5 + 8) * 4]
		
		; x[3] =
		;	Re(E[1]) [R12] + Im(E[1]) [RDI] +
		;	T(1+i) * (Re(e[1])[R8] + Im(o[1])[R11])
		;	[ R12 - TR8 - TR11 ; -RDI + TR11 - TR8 ]
		;	[ ST(4)-ST(2)-ST(1); -ST(3)+ST(1)-ST(2)]
	fld ST(3)
	fsub ST(0), ST(2)
	fsub ST(0), ST(1)
	fstp dword ptr [RCX + 3 * 4]
	
	fld ST(0)
	fsub ST(0), ST(2)
	fsub ST(0), ST(3)
	fstp dword ptr [RCX + (3 + 8) * 4]
		
		; x[7] =
		;	Re(E[1]) [R12] + Im(E[1]) [RDI] -
		;	T(1-i) * (Re(e[1])[R8] + Im(o[1])[R11])
		;	[ R12 + TR8 + TR11 ; -RDI - TR11 + TR8 ]
		;	[ ST(4)+ST(2)+ST(1); -ST(3)-ST(1)+ST(2)]
	fld ST(3)
	fadd ST(0), ST(2)
	fadd ST(0), ST(1)
	fstp dword ptr [RCX + 7 * 4]
	
	fld ST(1)
	fsub ST(0), ST(1)
	fsub ST(0), ST(3)
	fstp dword ptr [RCX + (7 + 8) * 4]

	pop RSI
	pop RDI
	pop R15
	pop R14
	pop R13
	pop R12
	pop R11
	pop R10

	ret
CalculateSpectrum ENDP
; -------------------------------------------------------------------------------------	;
; void RecoverSignal(signal_type* Signal, spectrum_type* Spectrum)						;
;	Обратное преобразование Фурье. Вычисляет сигнал Signal по спектру Spectrum			;
;	Типы данных spectrum_type и signal_type, а так же размер сигнала					;
;	определяются в файле Tuning.h														;
; -------------------------------------------------------------------------------------	;
RecoverSignal PROC	; [RCX] - Signal
					; [RDX] - Spectrum
	sub RSP, (8 * 2) * 8; allocating space for 8 complex numbers

	finit
	fld dword ptr [RDX				]	; loading Re(x[0])
	fld dword ptr [RDX + 4 * 4		]	; loading Re(x[4])
	fld ST(1)							; ST0 = Re(x[0])
	fadd ST(0), ST(1)					; calculating Re(e[0]) = Re(x[0]) + Re(x[4])
	fld ST(2)							; ST0 = Re(x[0])
	fsub ST(0), ST(2)					; calculating Re(e[1]) = Re(x[0]) - Re(x[4])
	fstp dword ptr [RSP + 1 * 4	]		; storing Re(e[1])
	fstp dword ptr [RSP			]		; storing Re(e[0])
	finit								; clearing ST

		; Im(e[0]), Im(e[1]) are optimized out

	fld dword ptr [RDX + 2 * 4		]	; loading Re(x[2])
	fld dword ptr [RDX + 6 * 4		]	; loading Re(x[6])
	fld ST(1)							; ST0 = Re(x[2])
	fadd ST(0), ST(1)					; calculating Re(o[0]) = Re(x[2]) + Re(x[6])
		; calculating Re(o[1]) = Re(x[2]) - Re(x[6])
		; o[1] * w(1, n) = o[1] * (-i)
		; Re(x[2]) - Re(x[6]) = -Im(o[1])
		; Im(o[1]) is therefore optimized out
	fstp dword ptr [RSP + 2 * 4 ]		; storing Re(o[0])
	finit								; clearing ST
	
	fld dword ptr [RDX + (2 + 8) * 4 ]	; loading Im(x[2])
	fld dword ptr [RDX + (6 + 8) * 4 ]	; loading Im(x[6])
	fld ST(1)
		; calculating Im(o[1]) = Im(x[2]) - Im(x[6])
		; o[1] * w(1, n) = o[1] * (-i)
		; Im(x[2]) - Im(x[6]) = Re(o[1])
	fsub ST(0), ST(2)
	fstp dword ptr [RSP + 3 * 4		  ]	; storing Re(o[1])
	finit								; clearing ST

		; [RSP + x * 4] + [RSP + (x + 8) * 4]*i :
		;	x = 0 -> e[0]
		;	x = 1 -> e[1]
		;	x = 2 -> o[0]
		;	x = 3 -> o[1]

	fld dword ptr [RSP		 ]		; loading Re(e[0])
	fld dword ptr [RSP + 2 * 4 ]	; loading Re(o[0])
	fld ST(1)						; ST0 = Re(e[0])
	fadd ST(0), ST(1)				; calculating Re(E[0])
	fld ST(2)						; ST0 = Re(e[0])
	fsub ST(0), ST(2)				; calculating Re(E[2])
	fstp dword ptr [RSP + 2 * 4 ]	; storing Re(E[2])
	fstp dword ptr [RSP			]	; storing Re(E[0])
	finit
	
	fld dword ptr [RSP + 1 * 4 ]	; loading Re(e[1])
	fld dword ptr [RSP + 3 * 4 ]	; loading Re(o[1])
	fld ST(1)						; ST0 = Re(e[1])
	fadd ST(0), ST(1)				; calculating Re(E[1])
	fld ST(2)						; ST0 = Re(e[1])
	fsub ST(0), ST(2)				; calculating Re(E[3])
	fstp dword ptr [RSP + 3 * 4 ]	; storing Re(E[3])
	fstp dword ptr [RSP + 1 * 4 ]	; storing Re(E[1])
	finit

		; Im(E[0..3]) have no effect on Re(X[k]) and are optimized out
	
		; processing odd indexes
	finit
	fld dword ptr [RDX + 1 * 4		]	; loading Re(x[1])
	fld dword ptr [RDX + 5 * 4		]	; loading Re(x[5])
	fld ST(1)							; ST0 = Re(x[1])
	fadd ST(0), ST(1)					; calculating Re(e[0]) = Re(x[1]) + Re(x[5])
	fld ST(2)							; ST0 = Re(x[1])
	fsub ST(0), ST(2)					; calculating Re(e[1]) = Re(x[1]) - Re(x[5])
	fstp dword ptr [RSP + 5 * 4	]		; storing Re(e[1])
	fstp dword ptr [RSP + 4 * 4]		; storing Re(e[0])
	finit								; clearing ST

	fld dword ptr [RDX + (1 + 8) * 4 ]	; loading Im(x[1])
	fld dword ptr [RDX + (5 + 8) * 4 ]	; loading Im(x[5])
	fld ST(1)							; ST0 = Im(x[1])
	fadd ST(0), ST(1)					; calculating Im(e[0]) = Im(x[1]) + Im(x[5])
	fld ST(2)							; ST0 = Im(x[1])
	fsub ST(0), ST(2)					; calculating Im(e[1]) = Im(x[1]) - Im(x[5])
	fstp dword ptr [RSP + (5 + 8) * 4 ] ; storing Im(e[1])
	fstp dword ptr [RSP + (4 + 8) * 4 ]	; storing Im(e[0])
	finit								; clearing ST

	fld dword ptr [RDX + 3 * 4		]	; loading Re(x[3])
	fld dword ptr [RDX + 7 * 4		]	; loading Re(x[7])
	fld ST(1)							; ST0 = Re(x[3])
	fadd ST(0), ST(1)					; calculating Re(o[0]) = Re(x[3]) + Re(x[7])
	fld ST(2)							; ST0 = Re(x[7])
		; calculating Re(o[1]) = Re(x[3]) - Re(x[7])
		; o[1] * w(1, n) = o[1] * (-i)
		; Re(x[3]) - Re(x[7]) = -Im(o[1]*w)
		; Reverse butterfly method -> Im(o[1]*w) = Im(x[3]) - Im(x[7])
	fsub ST(0), ST(2)
	fstp dword ptr [RSP + (7 + 8) * 4 ]	; storing Im(o[1]*w)
	fstp dword ptr [RSP + 6 * 4 ]		; storing Re(o[0]*w)
	finit								; clearing ST
	
	fld dword ptr [RDX + (3 + 8) * 4 ]	; loading Im(x[3])
	fld dword ptr [RDX + (7 + 8) * 4 ]	; loading Im(x[7])
	fld ST(1)							; ST0 = Im(x[3])
	fadd ST(0), ST(1)					; calculating Im(o[0]) = Im(x[3]) + Im(x[7])
	fld ST(1)							; ST0 = Im(x[7])
		; calculating Im(o[1]) = Im(x[3]) - Im(x[7])
		; o[1] * w(1, n) = o[1] * (-i)
		; Im(x[3]) - Im(x[7]) = Re(o[1]*w)
		; Reverse butterfly method -> Re(o[1]*w) = Im(x[7]) - Im(x[3])
	fsub ST(0), ST(3)
	fstp dword ptr [RSP + 7 * 4		  ]	; storing Re(o[1]*w)
	fstp dword ptr [RSP + (6 + 8) * 4 ]	; storing Im(o[0]*w)
	finit								; clearing ST

		; [RSP + x * 4] + [RSP + (x + 8) * 4]*i :
		;	x = 4 -> e[0]
		;	x = 5 -> e[1]
		;	x = 6 -> o[0]*w
		;	x = 7 -> o[1]*w

	fld dword ptr [RSP + 4 * 4 ]	; loading Re(e[0])
	fld dword ptr [RSP + 6 * 4 ]	; loading Re(o[0])
	fld ST(1)						; ST0 = Re(e[0])
	fadd ST(0), ST(1)				; calculating Re(O[0])
	fstp dword ptr [RSP + 4 * 4 ]	; storing Re(O[0])
	finit
		; O[0] is stored
	
	fld dword ptr [RSP + 5 * 4 ]	; loading Re(e[1])
	fld dword ptr [RSP + 7 * 4 ]	; loading Re(o[1])
	fld ST(1)						; ST0 = Re(e[1])
	fadd ST(0), ST(1)				; calculating Re(O[1])
	fld ST(2)						; ST0 = Re(e[1])
	fsub ST(0), ST(2)				; calculating Re(O[3])
	fstp dword ptr [RSP + 7 * 4 ]	; storing Re(O[3])
	fstp dword ptr [RSP + 5 * 4 ]	; storing Re(O[1])
	finit
		; Re(O[1]) and Re(O[3])
		; Im(O[2]) * (-i) = (e[0] - o[0]) * (-i) = Re(O[2]*w)
		; Reverse butterfly method -> Re(O[2]*w) = o[0] - e[0]
	fld dword ptr [RSP + (4 + 8) * 4 ]	; loading Im(e[0])
	fld dword ptr [RSP + (6 + 8) * 4 ]	; loading Im(o[0])
	fld ST(0)							; ST0 = Im(o[0])
	fsub ST(0), ST(2)					; calculating Re(O[2])
	fstp dword ptr [RSP + 6 * 4 ]		; storing Re(O[2])
	finit
	
	fld dword ptr [RSP + (5 + 8) * 4 ]	; loading Im(e[1])
	fld dword ptr [RSP + (7 + 8) * 4 ]	; loading Im(o[1])
	fld ST(1)							; ST0 = Im(e[1])
	fadd ST(0), ST(1)					; calculating Im(O[1])
	fld ST(2)							; ST0 = Im(e[1])
	fsub ST(0), ST(2)					; calculating Im(O[3])
	fstp dword ptr [RSP + (7 + 8) * 4 ]	; storing Im(O[3])
	fstp dword ptr [RSP + (5 + 8) * 4 ]	; storing Im(O[1])

		; Stack:
		; [RSP + x * 4] + [RSP + (x + 8) * 4] :
		;	[0..3] -> E[0..3]
		;	4 -> O[0]
		;	5 -> O[1] / W(1, 8)
		;	6 -> Re(O[2])
		;	7 -> O[3] / W(3, 8)

		; storing sqrt(2) / 2 in ST(0)
	finit
	mov dword ptr [RCX + 1 * 4 ], 2
	fild dword ptr [RCX + 1 * 4]
	fsqrt
	fild dword ptr [RCX + 1 * 4]
	fdivp
	fld dword ptr [RSP + (5 + 8) * 4]	; loading Im(O[1])
	fld dword ptr [RSP + 5 * 4]			; loading Re(O[1])
		; multiplying O[1] by W(1, 8) is the same as calculating
		;	Re(O[1]*w) = Im(O[1]) + Re(O[1])
		;	Im(O[1]*w) = Im(O[1]) - Re(O[1])- optimized out
		; Reverse butterfly method -> Re(O[1] * (1 + i)) = -Im(O[1]) + Re(O[1])
	fsub ST(0), ST(1)					; calculating Re(O[1]*w)
	fmul ST(0), ST(2)
	fstp dword ptr [RSP + 5 * 4]		; storing Re(O[1]*w)
	fstp ST(0)
	
	fld dword ptr [RSP + 7 * 4]			; loading Re(O[3])
	fld dword ptr [RSP + (7 + 8) * 4]	; loading Im(O[3])
		; multiplying O[3] by W(3, 8) is the same as calculating
		;	Re(O[3]*w) = Im(O[3]) - Re(O[3])
		;	Im(O[3]*w) = Im(O[3]) + Re(O[3])- optimized out
		; Reverse butterfly method -> Re(O[3] * (i - 1)) = -Im(O[1]) - Re(O[1])
	fldz
	fsub ST(0), ST(1)
	fsub ST(0), ST(2)					; calculating Re(O[3]*w)
	fmul ST(0), ST(3)
	fstp dword ptr [RSP + 7 * 4]		; storing Re(O[3]*w)
	
	finit
	fld1				; ST(0) = 1
	push 0				; M(0) = 8
	mov dword ptr [RSP], 8;
	fild dword ptr [RSP];
	fdivp ST(1), ST(0)	; ST(0) = 1 / 8

	fld dword ptr [RSP + 1 * 4 + 4]	; loading E[0]
	fld dword ptr [RSP + 5 * 4 + 4]	; loading O[0]
	fld ST(1)						; ST(0) = E[0]
	fadd ST(0), ST(1)				; X[0] = E[0] + O[0]
	fmul ST(0), ST(3)				; X[0] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[0])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 7], R8b			;
	fld ST(1)						; ST(0) = E[0]
	fsub ST(0), ST(1)				; X[4] = E[0] - O[0]
	fmul ST(0), ST(3)				; X[4] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[4])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 4], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 2 * 4 + 4]	; loading E[1]
	fld dword ptr [RSP + 6 * 4 + 4]	; loading O[1]
	fld ST(1)						; ST(0) = E[1]
	fadd ST(0), ST(1)				; X[1] = E[1] + O[1]
	fmul ST(0), ST(3)				; X[1] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[1])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 1], R8b		;
	fld ST(1)						; ST(0) = E[1]
	fsub ST(0), ST(1)				; X[5] = E[1] - O[1]
	fmul ST(0), ST(3)				; X[5] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[5])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 5], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 3 * 4 + 4]	; loading E[2]
	fld dword ptr [RSP + 7 * 4 + 4]	; loading O[2]
	fld ST(1)						; ST(0) = E[2]
	fadd ST(0), ST(1)				; X[2] = E[2] + O[2]
	fmul ST(0), ST(3)				; X[2] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[2])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 2], R8b		;
	fld ST(1)						; ST(0) = E[2]
	fsub ST(0), ST(1)				; X[6] = E[2] - O[2]
	fmul ST(0), ST(3)				; X[6] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[6])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 6], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 4 * 4 + 4]	; loading E[3]
	fld dword ptr [RSP + 8 * 4 + 4]	; loading O[3]
	fld ST(1)						; ST(0) = E[3]
	fadd ST(0), ST(1)				; X[3] = E[3] + O[3]
	fmul ST(0), ST(3)				; X[3] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[3])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 3], R8b		;
	fld ST(1)						; ST(0) = E[3]
	fsub ST(0), ST(1)				; X[7] = E[3] - O[3]
	fmul ST(0), ST(3)				; X[7] * (1 / N)
	fistp dword ptr [RSP]			; storing Re(X[7])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 7], R8b		;

	; fld ST(0)
	; fsub ST(0), ST(2)
	; fmul ST(0), ST(3)
	; fistp dword ptr [RCX + 4]
	; fstp ST(0)
	; fstp ST(0)
	; push 1

	add RSP, (8 * 2) * 8 + 8;

	ret
RecoverSignal ENDP
END
