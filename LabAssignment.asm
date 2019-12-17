; -------------------------------------------------------------------------------------	;
;	Лабораторная работа №2 по курсу "Программирование на языке ассемблера"				;
;	Вариант №1.2																		;
;	Выполнил студент Макеев Владислав													;
;																						;
;	Исходный модуль LabAssignment.asm													;
;	Содержит функции на языке ассемблера, разработанные в соответствии с заданием		;
; -------------------------------------------------------------------------------------	;
;	Задание: Реализовать прямое и обратное преобразования Фурье
;	Формат данных сигнала: __int8
;	Формат данных спектра: float
;	Размер (количество отсчетов) сигнала и спектра: 8
;	Способ реализации: DFT 2x2 + 2 бабочки
;	Отсчеты спектра являются комплексными числами. Причем действительные части хранятся
;	в первой половине массива, а мнимые - во второй

.DATA
	; Вектора-множители для выбора направления Spectrum <-> Signal
	NormalMV	dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
	ImMV		dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	3f800000h
				dword	0bf800000h
				dword	0bf800000h
				dword	0bf800000h
				dword	0bf800000h
	XImMV		dword	0bf800000h
				dword	0bf800000h
				dword	3f800000h
				dword	3f800000h
				dword	0bf800000h
				dword	0bf800000h
				dword	0bf800000h
				dword	0bf800000h
	; Вектора-множители для векторов oU (Вектор O верхней части бабочки), oB (Нижней бабочки) и O
		; ReRMV = < 0, 0, 0, 0,			1, 1, -1, -1 >
	ReRMV		dword 0bf800000h
				dword 0bf800000h
				dword 3f800000h
				dword 3f800000h
				dword 00000000h
				dword 00000000h
				dword 00000000h
				dword 00000000h
		; ImRMV = < -1, -1, 1, 1,			0, 0, 0, 0 >
	ImRMV		dword 00000000h
				dword 00000000h
				dword 00000000h
				dword 00000000h
				dword 3f800000h
				dword 3f800000h
				dword 0bf800000h
				dword 0bf800000h
		; ReIMV = < -1, -1, 1, 1,			0, 0, 0, 0 >
	ReIMV		dword 00000000h
				dword 00000000h
				dword 00000000h
				dword 00000000h
				dword 0bf800000h
				dword 0bf800000h
				dword 3f800000h
				dword 3f800000h
		; XReRMV = < 1, s, 0, -s, -1, -s, 0, s >
	XReRMV	dword	00000000h
				dword	00000000h
				dword	0bf800000h
				dword	3f800000h
				dword	3f3504f3h
				dword	0bf3504f7h
				dword	0bf3504f7h
				dword	3f3504f3h
		; XReIMV = < 0, s, 1, s, 0, -s, -1, -s >
	XReIMV	dword	0bf800000h
				dword	3f800000h
				dword	00000000h
				dword	00000000h
				dword	0bf3504f7h
				dword	3f3504f3h
				dword	0bf3504f7h
				dword	3f3504f3h
		; XImRMV = [ 1, -1, 0, 0, s, -s, s, -s ]
	XImRMV		dword	3f800000h
				dword	0bf800000h
				dword	00000000h
				dword	00000000h
				dword	3f3504f3h
				dword	0bf3504f7h
				dword	3f3504f3h
				dword	0bf3504f7h
	; Вектор чисел 1/8
	const18		dword	3e000000h
				dword	3e000000h
				dword	3e000000h
				dword	3e000000h
				dword	3e000000h
				dword	3e000000h
				dword	3e000000h
				dword	3e000000h

.CODE

ProcessLOandLE PROC							;
											; ymm0/1=Re/Im	@ x3/7,x2/6,x1/5,x0/4
	vaddps xmm2, xmm0, xmm1					; ymm2 = Re/Im	@ o0b, o0u, e0b, e0u,	0, 0, 0, 0
	vsubps xmm3, xmm0, xmm1					; ymm3 = Re/Im	@ o1b, o1u, e1b, e1u,	0, 0, 0, 0
	vperm2f128 ymm4, ymm2, ymm3, 00100000b	; ymm4 = Re/Im	@ o0B, o0U, e0B, e0U,	o1B, o1U, e1B, e1U
	vpermilps ymm5, ymm4, 00010001b			; ymm5 = Re/Im	@ o0U, o0B, o0U, o0B,	o1U, o1B, o1U, o1B
	vpermilps ymm6, ymm4, 10111011b			; ymm6 = Re/Im	@ e0U, e0B, e0U, e0B,	e1U, e1B, e1U, e1B
	ret
ProcessLOandLE ENDP

FourierTransform PROC	; Вспомогательная функция, содержащая общие для вычисления спектра и сигнала методы
		; [RCX]- результирующий вектор
		; Обозначения:
		;	e[i]- верхний вектор длины 4
		;	o[i]- нижний вектор длины 4
		;	E[i]- верхний вектор длины 8
		;	O[i]- нижний вектор длины 8
		
		; Символ @ обозначает то же, что и пара { [, ] } с тем различием,
		;	что элементы отображаются в обратном порядке
		; ymm0 = Re @ x3, x2, x1, x0, 0, 0, 0, 0
		; ymm1 = Re @ x7, x6, x5, x4, 0, 0, 0, 0
		; ymm12 = [ 0 .. Im(x3..x0)]
		; ymm14 = [ 0 .. Im(x7..x4)]
	sub RSP, 4 * 8 * 2 * 2;
	call ProcessLOandLE; Вызов функции вычисления Re(o-s, e-s)
	vmovups ymm7, ymm5	; ymm7 = [ Re(o1B, o1U, o1B, o1U,	o0B, o0U, o0B, o0U ) ]
	vmovups ymm8, ymm6	; ymm8 = [ Re(e1B, e1U, e1B, e1U,	e0B, e0U, e0B, e0U ) ]
	vmovups ymm0, ymm12;
	vmovups ymm1, ymm14;
	call ProcessLOandLE; Вызов функции вычисления Im(o-s, e-s)
		; ymm5 = Im(o-s)
		; ymm6 = Im(e-s)
		; ymm7 = Re(o-s)
		; ymm8 = Re(e-s)
	; Calculating o*w
	
	vmovups ymm0, [ ReRMV ]	; Загрузка векторов-множителей
	vmovups ymm1, [ ReIMV ]	;
	vmovups ymm2, [ ImRMV ]	;
	
	vmulps ymm1, ymm1, ymm15; Умножение множителей Im(MV) на Sign(Im(w(1, 8))*(-1)
	vmulps ymm2, ymm2, ymm15;
	
	vmulps ymm9, ymm7, ymm0		;
	vmulps ymm10, ymm5, ymm1	;
	vaddps ymm11, ymm9, ymm10	;
		; ymm11 = Re(o-s * w)
	vmulps ymm9, ymm7, ymm2		;
	vmulps ymm10, ymm5, ymm0	; ImIMV = ReRMV
	vaddps ymm12, ymm9, ymm10	;
		; ymm12 = Im(o-s * w)
		
		; ymm12 = Im(o-s * w)
		; ymm6 = Im(e-s)
		; ymm11 = Re(o-s * w)
		; ymm8 = Re(e-s)
	
	vaddps ymm9, ymm8, ymm11			; ymm9 =	Re(O[1], E[1], O[3], E[3], O[0], E[0], O[2], E[2])
	vaddps ymm10, ymm6, ymm12			; ymm10 =	Im(O[1], E[1], O[3], E[3], O[0], E[0], O[2], E[2])

	vpermilps ymm5, ymm9, 10100000b		; ymm5 =	Re(E[1], E[1], E[3], E[3], E[0], E[0], E[2], E[2])
	vpermilps ymm6, ymm10, 10100000b	; ymm6 =	Im(E[1], E[1], E[3], E[3], E[0], E[0], E[2], E[2])
	vpermilps ymm7, ymm9, 11110101b		; ymm7 =	Re(O[1], O[1], O[3], O[3], O[0], O[0], O[2], O[2])
	vpermilps ymm8, ymm10, 11110101b	; ymm8 =	Im(O[1], O[1], O[3], O[3], O[0], O[0], O[2], O[2])

	vmovups ymm0, [ XReRMV ]	;
	vmovups ymm1, [ XReIMV ]	;
	vmovups ymm2, [ XImRMV ]	;
	
	vmulps ymm1, ymm1, ymm13; Умножение множителей Im(MV) на Sign(Im(w(1, 8))*(-1)
	vmulps ymm2, ymm2, ymm13;

	vmulps ymm9, ymm7, ymm0		;
	vmulps ymm10, ymm8, ymm1	;
	vaddps ymm11, ymm9, ymm10	;
		; ymm11 = Re(O-s)

	vmulps ymm9, ymm7, ymm2		;
	vmulps ymm10, ymm8, ymm0	; ImIMV = ReRMV
	vaddps ymm12, ymm9, ymm10	;
		; ymm12 = Im(O-s)
		
		; ymm5 =	Re(E[1], E[1], E[3], E[3], E[0], E[0], E[2], E[2])
		; ymm6 =	Im(E[1], E[1], E[3], E[3], E[0], E[0], E[2], E[2])
		; ymm11 =	Re[ O1, O5, O3, O7, O0, O4, O2, O6 ]
		; ymm12 =	Im[ O1, O5, O3, O7, O0, O4, O2, O6 ]

	vaddps ymm0, ymm5, ymm11;
	vaddps ymm1, ymm6, ymm12;
	vmovups ymmword ptr [RSP], ymm0
	vmovups ymmword ptr [RSP + 4 * 8], ymm1

	; Запись результатов
	;	[ 0 1 2 3	4 5 6 7 ] -> [ 6 2 4 0	7 3 5 1 ]
	
	vmovss dword ptr [ RCX + 4 * 6 ], xmm0	;
	vperm2f128 ymm0, ymm0, ymm0, 00000001b	;
	vmovss dword ptr [ RCX + 4 * 7 ], xmm0	;
	vpermilps ymm0, ymm0, 10010011b			;
	vmovss dword ptr [ RCX + 4 * 1 ], xmm0	;
	vperm2f128 ymm0, ymm0, ymm0, 00000001b	;
	vmovss dword ptr [ RCX ], xmm0			;
	vpermilps ymm0, ymm0, 10010011b			;
	vmovss dword ptr [ RCX + 4 * 4 ], xmm0	;
	vperm2f128 ymm0, ymm0, ymm0, 00000001b	;
	vmovss dword ptr [ RCX + 4 * 5], xmm0	;
	vpermilps ymm0, ymm0, 10010011b			;
	vmovss dword ptr [ RCX + 4 * 3 ], xmm0	;
	vperm2f128 ymm0, ymm0, ymm0, 00000001b	;
	vmovss dword ptr [ RCX + 4 * 2], xmm0	;
	vpermilps ymm0, ymm0, 10010011b			;
	
	vmovss dword ptr [ RCX + 4 * (8 + 6) ], xmm1	;
	vperm2f128 ymm1, ymm1, ymm1, 00000001b			;
	vmovss dword ptr [ RCX + 4 * (8 + 7) ], xmm1	;
	vpermilps ymm1, ymm1, 10010011b					;
	vmovss dword ptr [ RCX + 4 * (8 + 1) ], xmm1	;
	vperm2f128 ymm1, ymm1, ymm1, 00000001b			;
	vmovss dword ptr [ RCX + 4 * 8 ], xmm0			;
	vpermilps ymm1, ymm1, 10010011b					;
	vmovss dword ptr [ RCX + 4 * (8 + 4) ], xmm1	;
	vperm2f128 ymm1, ymm1, ymm1, 00000001b			;
	vmovss dword ptr [ RCX + 4 * (8 + 5) ], xmm1	;
	vpermilps ymm1, ymm1, 10010011b					;
	vmovss dword ptr [ RCX + 4 * (8 + 3) ], xmm1	;
	vperm2f128 ymm1, ymm1, ymm1, 00000001b			;
	vmovss dword ptr [ RCX + 4 * (8 + 2) ], xmm1	;
	vpermilps ymm1, ymm1, 10010011b					;

	add RSP, 4 * 8 * 2 * 2;
	ret
FourierTransform ENDP

; -------------------------------------------------------------------------------------	;
; void CalculateSpectrum(spectrum_type* Spectrum, signal_type* Signal)					;
;	Прямое преобразование Фурье. Вычисляет спектр Spectrum по сигналу Signal			;
;	FT_SIGNAL_SIZE = 8																	;
;	spectrum_type = float																;
;	signal_type = __int8																;
; -------------------------------------------------------------------------------------	;
CalculateSpectrum PROC	; [RCX] - Spectrum
						; [RDX] - Signal
	push RCX
	push RDX
	mov eax, 1	; Проверка поддержки
	cpuid		;	AVX процессором
	bt ecx, 28
	jc AVX_SUPPORTED
	pop RDX
	pop RCX
	push R12		; Сохранение non-volatile регистров
	push R13		;	в стек перед началом работы
	push R14		;	подпрограммы с целью
	push R15		;	их дальнейшего восстановления
	push RDI		;	перед возвращением
	push RSI		;	в вызывающую функцию

	mov R15, [RDX]	; Загрузка сигнала в регистр R15

	xor R8, R8		; Обнуление регистров
	xor R9, R9		;	для дальнейшей записи
	xor R10, R10	;	элемента сигнала
	xor R11, R11	;	x[i] в байтовую часть
	xor R12, R12	;	регистра R8+i
	xor R13, R13	;	посредством переноса
	xor R14, R14	;	из R15b

	mov R8b, R15b	; Перенос байта x[0] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R12b, R15b	; Перенос байта x[1] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R9b, R15b	; Перенос байта x[2] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R13b, R15b	; Перенос байта x[3] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R10b, R15b	; Перенос байта x[4] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R14b, R15b	; Перенос байта x[5] из R15
	shr r15, 8		;	с последующим сдвигом >> 8
	mov R11b, R15b	; Перенос байта x[6] из R15
	shr r15, 8		;	с последующим сдвигом >> 8

		; Информация в регистрах:
		;	R8-  x[0]
		;	R9-  x[2]
		;	R10- x[4]
		;	R11- x[6]
		;	R12- x[1]
		;	R13- x[3]
		;	R14- x[5]
		;	R15- x[7]
		
	sal r8, 7 * 8	; В результате
	sar r8, 7 * 8	;	переноса из R15b
	sal r9, 7 * 8	;	в регистрах записаны
	sar r9, 7 * 8	;	значения x[i].
	sal r10, 7 * 8	; Следует, однако же,
	sar r10, 7 * 8	;	учитывать тот факт, что
	sal r11, 7 * 8	;	что перенос осуществлялся
	sar r11, 7 * 8	;	не с помощью MOVSX,
	sal r12, 7 * 8	;	а значит остальные биты
	sar r12, 7 * 8	;	в случае записи в регистр
	sal r13, 7 * 8	;	отрицательного числа
	sar r13, 7 * 8	;	следует инвертировать.
	sal r14, 7 * 8	; Для этого используется
	sar r14, 7 * 8	;	побитовый сдвиг влево
	sal r15, 7 * 8	;	с последующим знаковым
	sar r15, 7 * 8	;	сдвигом вправо

		; Обработка чётных индексов

	mov RSI, R10	; Загрузка x[4] на RSI
	add R10, R8		; Вычисление e[0] и его запись на R10
	sub R8, RSI		; Вычисление e[1] и его запись на R8
	mov RSI, R9		; Загрузка x[2] на RSI
	add R9, R11		; Вычисление o[0] и его запись на R9
	sub R11, RSI	; Вычисление мнимой части o[1] и её запись на R11

		; Информация в регистрах:
		;	R8-  Re(e[1])	Не имеет мнимой части
		;	R9-  Re(o[0])	Не имеет мнимой части
		;	R10- Re(e[0])	Не имеет мнимой части
		;	R11- Im(o[1])	Не имеет действительной части

		; E[0] = o[0] + e[0] => E[0] = Re(E[0])
		; E[2] = o[0] - e[0] => E[2] = Re(E[2])

	mov RSI, R9		; Загрузка o[0] на RSI
	add r9, R10		; Вычисление E[0] и его запись на R9
	sub r10, RSI	; Вычисление E[2] и его запись на R10

		; E[1] = Re(e[1]) + Im(o[1])
		; E[3] = Re(e[1]) - Im(o[1])

		; Информация в регистрах:
		;	R8-  Re(e[1])		Re(E[1]) и Re(E[3])
		;	R9-  Re(E[0])		Не имеет мнимой части
		;	R10- Re(E[2])		Не имеет мнимой части
		;	R11- Im(o[1])		Im(E[1]) и -Im(E[3])

		; Обработка нечётных индексов

	mov RSI, R14		; Загрузка x[5] into RSI
	add R14, R12		; Вычисление e[0] и его запись на R14
	sub R12, RSI		; Вычисление e[1] и его запись на R12
	mov RSI, R13		; Загрузка x[3] into RSI
	add R13, R15		; Вычисление o[0] и его запись на R13
	sub R15, RSI		; Вычисление мнимой части o[1] и её запись на R15
	
		; Информация в регистрах:
		;	R12- Re(e[1])	Не имеет мнимой части
		;	R13- Re(o[0])	Не имеет мнимой части
		;	R14- Re(e[0])	Не имеет мнимой части
		;	R15- Im(o[1])	Не имеет действительной части

		; O[0] = o[0] + e[0] => O[0] = Re(O[0])
		; O[2] = o[0] - e[0] => O[2] = Re(O[2])

	mov RSI, R13		; Загрузка o[0] на RSI
	add R13, R14		; Вычисление O[0] и его запись на R13
	sub R14, RSI		; Вычисление O[2] и его запись на R14

		; O[1] = Re(e[1]) + Im(o[1])
		; O[3] = Re(e[1]) - Im(o[1])
		
	mov RDI, R11; RDI не имеет постфикса -b для записи в младший байт, а для вычислений нужен регистр R-b

		; Информация в регистрах:
		;	R12-  Re(O[1]) = Re(O[3])
		;	R13-  Re(O[0])					Не имеет мнимой части
		;	R14-  Re(O[2])					Не имеет мнимой части
		;	R15-  Im(O[1]) = -Im(O[3])
		;	R8-   Re(E[1]) = Re(E[3])
		;	R9-   Re(E[0])					Не имеет мнимой части
		;	R10-  Re(E[2])					Не имеет мнимой части
		;	RDI-  Im(E[1]) = -Im(E[3])

	finit	; Инициализация FPU

		; x[0] = Re(E[0]) + Re(O[0])
	mov R11, R9
	add R11, R13
	mov dword ptr [	RCX			], R11d		; spectrum[0] = x[0]
	fild dword ptr [RCX			]			; Конвертация spectrum[0] int -> float
	fstp dword ptr [RCX			]			; Запись конвертированного числа в spectrum[0]
	mov dword ptr [	RCX + 8 * 4	], 0		; Im(x[0]) = 0

		; x[4] = Re(E[0]) - Re(O[0])
	mov R11, R9
	sub R11, R13
	mov dword ptr [ RCX + 4 * 4], R11d		; spectrum[4] = x[4]
	fild dword ptr [RCX + 4 * 4]			; Конвертация spectrum[4] int -> float
	fstp dword ptr [RCX + 4 * 4]			; Запись конвертированного числа в spectrum[4]
	mov dword ptr [	RCX + (4 + 8) * 4], 0	; Im(x[4]) = 0
		
		; x[2] = Re(E[2]) - Im(O[2])
	mov dword ptr [	RCX +	 2 * 4		], R10d		; Re(spectrum[2]) = E[2]
	fild dword ptr [RCX +	 2 * 4		]			; Загрузка Re(spectrum[2]) на FPU
	fstp dword ptr [RCX +	 2 * 4		]			; Запись конвертированного во float Re(spectrum[2])
	mov R11, 0
	sub R11, R14
	mov dword ptr [RCX +	(2 + 8) * 4	], R11d		; Im(spectrum[2]) = (-1) * O[2]
	fild dword ptr [RCX +	(2 + 8) * 4	]			; Загрузка Im(spectrum[2]) на FPU
	fstp dword ptr [RCX +	(2 + 8) * 4	]			; Запись конвертированного во float Im(spectrum[2])
	
		; x[6] = Re(E[2]) + Im(O[2])
	mov dword ptr [RCX +	 6 * 4		], R10d		; Re(spectrum[6]) = E[2]
	fild dword ptr [RCX +	 6 * 4		]			; Загрузка Re(spectrum[6]) на FPU
	fstp dword ptr [RCX +	 6 * 4		]			; Запись конвертированного во float Re(spectrum[6])
	mov dword ptr [RCX +	(6 + 8) * 4	], R14d		; Im(spectrum[6]) = O[2]
	fild dword ptr [RCX +	(6 + 8) * 4	]			; Загрузка Im(spectrum[2]) на FPU
	fstp dword ptr [RCX +	(6 + 8) * 4	]			; Запись конвертированного во float Im(spectrum[2])

		; Загрузка sqrt(2) / 2 на ST(0)
	mov dword ptr [RCX + 1 * 4], 2
	fild dword ptr [RCX + 1 * 4]
	fsqrt
	fild dword ptr [RCX + 1 * 4]
	fdivp
	
		; x[1] =
		;	Re(E[1]) [R8] + Im(E[1]) [RDI] +
		;	T(1-i) * (Re(e[1])[R12] + Im(o[1])[R15])
		;	[ R8 + TR12 + TR15 ; RDI + TR15 - TR12 ]
		; TRegister := (sqrt(2) / 2) * Register

		; Загрузка R8, R12, R15 и RDI в стек FPU, умножение R12 и R15 на (sqrt(2) / 2)
	mov dword ptr [RCX + 1 * 4], R8d
	mov R11, RDI
	mov dword ptr [RCX + 3 * 4], R11d
	mov dword ptr [RCX + 5 * 4], R12d
	mov dword ptr [RCX + 7 * 4], R15d

	fild dword ptr [RCX + 1 * 4]
	fild dword ptr [RCX + 3 * 4]
	fild dword ptr [RCX + 5 * 4]
	fmul ST(0), ST(3)
	fild dword ptr [RCX + 7 * 4]
	fmul ST(0), ST(4)
	fld ST(3)

		; Данные в стеке:
		;	ST(0)- Accumulator (R8)- вершина стека используется для вычисления
		;	ST(1)- TR15
		;	ST(2)- TR12
		;	ST(3)- R11 (RDI)
		;	ST(4)- R8
		;	ST(5)- sqrt(2) / 2

		; x[1] =
		;	[ R8 + TR12 + TR15 ; RDI + TR15 - TR12 ]
		;	[ ST(4)+ST(2)+ST(1); ST(3)+ST(1)-ST(2)]
	fadd ST(0), ST(2)
	fadd ST(0), ST(1)
	fstp dword ptr [RCX + 1 * 4]

	fld ST(2)
	fadd ST(0), ST(1)
	fsub ST(0), ST(2)
	fstp dword ptr [RCX + (1 + 8) * 4]
		
		; x[5] =
		;	Re(E[1]) [R8] + Im(E[1]) [RDI] -
		;	T(1-i) * (Re(e[1])[R12] + Im(o[1])[R15])
		;	[ R8 - TR12 - TR15 ; RDI - TR15 + TR12 ]
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
		;	Re(E[1]) [R8] + Im(E[1]) [RDI] +
		;	T(1+i) * (Re(e[1])[R12] + Im(o[1])[R15])
		;	[ R8 - TR12 - TR15 ; -RDI + TR15 - TR12 ]
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
		;	Re(E[1]) [R8] + Im(E[1]) [RDI] -
		;	T(1-i) * (Re(e[1])[R12] + Im(o[1])[R15])
		;	[ R8 + TR12 + TR15 ; -RDI - TR15 + TR12 ]
		;	[ ST(4)+ST(2)+ST(1); -ST(3)-ST(1)+ST(2)]
	fld ST(3)
	fadd ST(0), ST(2)
	fadd ST(0), ST(1)
	fstp dword ptr [RCX + 7 * 4]
	
	fld ST(1)
	fsub ST(0), ST(1)
	fsub ST(0), ST(3)
	fstp dword ptr [RCX + (7 + 8) * 4]

		; Восстановление значений non-volatile регистров

	pop RSI
	pop RDI
	pop R15
	pop R14
	pop R13
	pop R12

	ret
	AVX_SUPPORTED:

	pop RDX
	pop RCX

	sub RSP, 4 * 8 * 2 * 2;

	vzeroall;
	vmovups xmm1, dword ptr [RDX + 4]	; Получение нужных
	vmovups xmm0, dword ptr [RDX]		;	данных из памяти
	vpmovsxbd xmm0, xmm0				;	и их запись в
	vpmovsxbd xmm1, xmm1				;	необходимом для
	vcvtdq2ps xmm0, xmm0				;	корректной работы
	vcvtdq2ps xmm1, xmm1				;	функции виде
	vpermilps ymm0, ymm0, 00011011b;
	vpermilps ymm1, ymm1, 00011011b;
		; ymm0 = [ 0 .. Re(x3..x0) ]
		; ymm1 = [ 0 .. Re(x7..x4) ]
		; ymm12 = [ 0 .. Im(x3..x0) ]
		; ymm14 = [ 0 .. Im(x7..x4) ]
	
	vmovups xmmword ptr [RSP], xmm0;
	vmovups xmmword ptr [RSP + 4 * 4], xmm1;
	vmovups ymm13, ymmword ptr [NormalMV];
	vmovups ymm15, ymmword ptr [NormalMV];
	
	push RCX;
	push RDX;
	mov RDX, RSP;
	call FourierTransform
	pop RDX;
	pop RCX;
	add RSP, 4 * 8 * 2 * 2;

	ret
CalculateSpectrum ENDP
; -------------------------------------------------------------------------------------	;
; void RecoverSignal(signal_type* Signal, spectrum_type* Spectrum)						;
;	Обратное преобразование Фурье. Вычисляет сигнал Signal по спектру Spectrum			;
;	FT_SIGNAL_SIZE = 8																	;
;	spectrum_type = float																;
;	signal_type = __int8																;
; -------------------------------------------------------------------------------------	;
RecoverSignal PROC	; [RCX] - Signal
					; [RDX] - Spectrum

	sub RSP, 4 * 8 * 7; Выделение места
	
	push RCX
	push RDX
	mov eax, 1	; Проверка поддержки
	cpuid		;	AVX процессором
	bt ecx, 28
	jc AVX_SUPPORTED
	pop RDX
	pop RCX
	sub RSP, (8 * 2) * 8; Выделение места для 8 комплексных чисел

	finit
	fld dword ptr [RDX				]	; Загрузка Re(x[0])
	fld dword ptr [RDX + 4 * 4		]	; Загрузка Re(x[4])
	fld ST(1)							; ST0 = Re(x[0])
	fadd ST(0), ST(1)					; Вычисление Re(e[0]) = Re(x[0]) + Re(x[4])
	fld ST(2)							; ST0 = Re(x[0])
	fsub ST(0), ST(2)					; Вычисление Re(e[1]) = Re(x[0]) - Re(x[4])
	fstp dword ptr [RSP + 1 * 4	]		; Запись Re(e[1])
	fstp dword ptr [RSP			]		; Запись Re(e[0])
	finit								; Инициализация FPU (Применяется для обнуления указателя стека)

		; Im(e[0]), Im(e[1]) отбрасываются по причине отсутствия влияния на Re(X)

	fld dword ptr [RDX + 2 * 4		]	; Загрузка Re(x[2])
	fld dword ptr [RDX + 6 * 4		]	; Загрузка Re(x[6])
	fld ST(1)							; ST0 = Re(x[2])
	fadd ST(0), ST(1)					; Вычисление Re(o[0]) = Re(x[2]) + Re(x[6])
		; Вычисление Re(o[1]) = Re(x[2]) - Re(x[6])
		; o[1] * w(1, n) = o[1] * (-i)
		; Re(x[2]) - Re(x[6]) = -Im(o[1])
		; Im(o[1]) отбрасывается по причине отсутствия влияния на Re(X)
	fstp dword ptr [RSP + 2 * 4 ]		; Запись Re(o[0])
	finit								; Инициализация FPU (Применяется для обнуления указателя стека)
	
	fld dword ptr [RDX + (2 + 8) * 4 ]	; Загрузка Im(x[2])
	fld dword ptr [RDX + (6 + 8) * 4 ]	; Загрузка Im(x[6])
	fld ST(0)							; ST0 = Im(x[6])
		; Вычисление Im(o[1]) = Im(x[2]) - Im(x[6])
		; o[1] * w(1, n) = o[1] * (-i)
		; Im(x[2]) - Im(x[6]) = Re(o[1])
		; Reverse butterfly method -> Re(o[1]) = Im(x[6]) - Im(x[2])
	fsub ST(0), ST(2)
	fstp dword ptr [RSP + 3 * 4		  ]	; Запись Re(o[1])
	finit								; clearing ST

		; [RSP + x * 4] + [RSP + (x + 8) * 4]*i :
		;	x = 0 -> e[0]
		;	x = 1 -> e[1]
		;	x = 2 -> o[0]
		;	x = 3 -> o[1]

	fld dword ptr [RSP		 ]		; Загрузка Re(e[0])
	fld dword ptr [RSP + 2 * 4 ]	; Загрузка Re(o[0])
	fld ST(1)						; ST0 = Re(e[0])
	fadd ST(0), ST(1)				; Вычисление Re(E[0])
	fld ST(2)						; ST0 = Re(e[0])
	fsub ST(0), ST(2)				; Вычисление Re(E[2])
	fstp dword ptr [RSP + 2 * 4 ]	; Запись Re(E[2])
	fstp dword ptr [RSP			]	; Запись Re(E[0])
	finit
	
	fld dword ptr [RSP + 1 * 4 ]	; Загрузка Re(e[1])
	fld dword ptr [RSP + 3 * 4 ]	; Загрузка Re(o[1])
	fld ST(1)						; ST0 = Re(e[1])
	fadd ST(0), ST(1)				; Вычисление Re(E[1])
	fld ST(2)						; ST0 = Re(e[1])
	fsub ST(0), ST(2)				; Вычисление Re(E[3])
	fstp dword ptr [RSP + 3 * 4 ]	; Запись Re(E[3])
	fstp dword ptr [RSP + 1 * 4 ]	; Запись Re(E[1])
	finit

		; Im(E[0..3]) отбрасываются по причине отсутствия влияния на Re(X)
	
		; processing odd indexes
	finit
	fld dword ptr [RDX + 1 * 4		]	; Загрузка Re(x[1])
	fld dword ptr [RDX + 5 * 4		]	; Загрузка Re(x[5])
	fld ST(1)							; ST0 = Re(x[1])
	fadd ST(0), ST(1)					; Вычисление Re(e[0]) = Re(x[1]) + Re(x[5])
	fld ST(2)							; ST0 = Re(x[1])
	fsub ST(0), ST(2)					; Вычисление Re(e[1]) = Re(x[1]) - Re(x[5])
	fstp dword ptr [RSP + 5 * 4	]		; Запись Re(e[1])
	fstp dword ptr [RSP + 4 * 4]		; Запись Re(e[0])
	finit								; clearing ST

	fld dword ptr [RDX + (1 + 8) * 4 ]	; Загрузка Im(x[1])
	fld dword ptr [RDX + (5 + 8) * 4 ]	; Загрузка Im(x[5])
	fld ST(1)							; ST0 = Im(x[1])
	fadd ST(0), ST(1)					; Вычисление Im(e[0]) = Im(x[1]) + Im(x[5])
	fld ST(2)							; ST0 = Im(x[1])
	fsub ST(0), ST(2)					; Вычисление Im(e[1]) = Im(x[1]) - Im(x[5])
	fstp dword ptr [RSP + (5 + 8) * 4 ] ; Запись Im(e[1])
	fstp dword ptr [RSP + (4 + 8) * 4 ]	; Запись Im(e[0])
	finit								; clearing ST

	fld dword ptr [RDX + 3 * 4		]	; Загрузка Re(x[3])
	fld dword ptr [RDX + 7 * 4		]	; Загрузка Re(x[7])
	fld ST(1)							; ST0 = Re(x[3])
	fadd ST(0), ST(1)					; Вычисление Re(o[0]) = Re(x[3]) + Re(x[7])
	fld ST(2)							; ST0 = Re(x[7])
		; Вычисление Re(o[1]) = Re(x[3]) - Re(x[7])
		; o[1] * w(1, n) = o[1] * (-i)
		; Re(x[3]) - Re(x[7]) = -Im(o[1]*w)
		; Reverse butterfly method -> Im(o[1]*w) = Im(x[3]) - Im(x[7])
	fsub ST(0), ST(2)
	fstp dword ptr [RSP + (7 + 8) * 4 ]	; Запись Im(o[1]*w)
	fstp dword ptr [RSP + 6 * 4 ]		; Запись Re(o[0]*w)
	finit								; clearing ST
	
	fld dword ptr [RDX + (3 + 8) * 4 ]	; Загрузка Im(x[3])
	fld dword ptr [RDX + (7 + 8) * 4 ]	; Загрузка Im(x[7])
	fld ST(1)							; ST0 = Im(x[3])
	fadd ST(0), ST(1)					; Вычисление Im(o[0]) = Im(x[3]) + Im(x[7])
	fld ST(1)							; ST0 = Im(x[7])
		; Вычисление Im(o[1]) = Im(x[3]) - Im(x[7])
		; o[1] * w(1, n) = o[1] * (-i)
		; Im(x[3]) - Im(x[7]) = Re(o[1]*w)
		; Reverse butterfly method -> Re(o[1]*w) = Im(x[7]) - Im(x[3])
	fsub ST(0), ST(3)
	fstp dword ptr [RSP + 7 * 4		  ]	; Запись Re(o[1]*w)
	fstp dword ptr [RSP + (6 + 8) * 4 ]	; Запись Im(o[0]*w)
	finit								; clearing ST

		; [RSP + x * 4] + [RSP + (x + 8) * 4]*i :
		;	x = 4 -> e[0]
		;	x = 5 -> e[1]
		;	x = 6 -> o[0]*w
		;	x = 7 -> o[1]*w

	fld dword ptr [RSP + 4 * 4 ]	; Загрузка Re(e[0])
	fld dword ptr [RSP + 6 * 4 ]	; Загрузка Re(o[0])
	fld ST(1)						; ST0 = Re(e[0])
	fadd ST(0), ST(1)				; Вычисление Re(O[0])
	fstp dword ptr [RSP + 4 * 4 ]	; Запись Re(O[0])
	finit
	
	fld dword ptr [RSP + 5 * 4 ]	; Загрузка Re(e[1])
	fld dword ptr [RSP + 7 * 4 ]	; Загрузка Re(o[1])
	fld ST(1)						; ST0 = Re(e[1])
	fadd ST(0), ST(1)				; Вычисление Re(O[1])
	fld ST(2)						; ST0 = Re(e[1])
	fsub ST(0), ST(2)				; Вычисление Re(O[3])
	fstp dword ptr [RSP + 7 * 4 ]	; Запись Re(O[3])
	fstp dword ptr [RSP + 5 * 4 ]	; Запись Re(O[1])
	finit
		; Re(O[1]) and Re(O[3])
		; Im(O[2]) * (-i) = (e[0] - o[0]) * (-i) = Re(O[2]*w)
		; Обратная бабочка использует комплексно-сопряжённые w:
		;	Re(O[2]*w) = o[0] - e[0]
	fld dword ptr [RSP + (4 + 8) * 4 ]	; Загрузка Im(e[0])
	fld dword ptr [RSP + (6 + 8) * 4 ]	; Загрузка Im(o[0])
	fld ST(0)							; ST0 = Im(o[0])
	fsub ST(0), ST(2)					; Вычисление Re(O[2])
	fstp dword ptr [RSP + 6 * 4 ]		; Запись Re(O[2])
	finit
	
	fld dword ptr [RSP + (5 + 8) * 4 ]	; Загрузка Im(e[1])
	fld dword ptr [RSP + (7 + 8) * 4 ]	; Загрузка Im(o[1])
	fld ST(1)							; ST0 = Im(e[1])
	fadd ST(0), ST(1)					; Вычисление Im(O[1])
	fld ST(2)							; ST0 = Im(e[1])
	fsub ST(0), ST(2)					; Вычисление Im(O[3])
	fstp dword ptr [RSP + (7 + 8) * 4 ]	; Запись Im(O[3])
	fstp dword ptr [RSP + (5 + 8) * 4 ]	; Запись Im(O[1])

		; Stack:
		; [RSP + x * 4] + [RSP + (x + 8) * 4] :
		;	[0..3] -> E[0..3]
		;	4 -> O[0]
		;	5 -> O[1] / W(1, 8)
		;	6 -> Re(O[2])
		;	7 -> O[3] / W(3, 8)

		; Запись sqrt(2) / 2 in ST(0)
	finit
	mov dword ptr [RCX + 1 * 4 ], 2
	fild dword ptr [RCX + 1 * 4]
	fsqrt
	fild dword ptr [RCX + 1 * 4]
	fdivp
	fld dword ptr [RSP + (5 + 8) * 4]	; Загрузка Im(O[1])
	fld dword ptr [RSP + 5 * 4]			; Загрузка Re(O[1])
		; multiplying O[1] by W(1, 8) is the same as Вычисление
		;	Re(O[1]*w) = Im(O[1]) + Re(O[1])
		;	Im(O[1]*w) = Im(O[1]) - Re(O[1])- optimized out
		; Reverse butterfly method -> Re(O[1] * (1 + i)) = -Im(O[1]) + Re(O[1])
	fsub ST(0), ST(1)					; Вычисление Re(O[1]*w)
	fmul ST(0), ST(2)
	fstp dword ptr [RSP + 5 * 4]		; Запись Re(O[1]*w)
	fstp ST(0)
	
	fld dword ptr [RSP + 7 * 4]			; Загрузка Re(O[3])
	fld dword ptr [RSP + (7 + 8) * 4]	; Загрузка Im(O[3])
		; Умножение O[3] на W(3, 8) даёт тот же результат, что и вычисление
		;	Re(O[3]*w) = Im(O[3]) - Re(O[3])
		;	Im(O[3]*w) = Im(O[3]) + Re(O[3])- Отбрасывается по причине отсутствия влияния на Re(X)
		; Обратная бабочка использует комплексно-сопряжённые w:
		;	Re(O[3] * (i - 1)) = -Im(O[1]) - Re(O[1])
	fldz
	fsub ST(0), ST(1)
	fsub ST(0), ST(2)					; Вычисление Re(O[3]*w)
	fmul ST(0), ST(3)
	fstp dword ptr [RSP + 7 * 4]		; Запись Re(O[3]*w)
	
	finit
	fld1				; ST(0) = 1
	push 0				; M(0) = 8
	mov dword ptr [RSP], 8;
	fild dword ptr [RSP];
	fdivp ST(1), ST(0)	; ST(0) = 1 / 8

	fld dword ptr [RSP + 1 * 4 + 4]	; Загрузка E[0]
	fld dword ptr [RSP + 5 * 4 + 4]	; Загрузка O[0]
	fld ST(1)						; ST(0) = E[0]
	fadd ST(0), ST(1)				; X[0] = E[0] + O[0]
	fmul ST(0), ST(3)				; X[0] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[0])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX], R8b			;
	fld ST(1)						; ST(0) = E[0]
	fsub ST(0), ST(1)				; X[4] = E[0] - O[0]
	fmul ST(0), ST(3)				; X[4] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[4])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 4], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 2 * 4 + 4]	; Загрузка E[1]
	fld dword ptr [RSP + 6 * 4 + 4]	; Загрузка O[1]
	fld ST(1)						; ST(0) = E[1]
	fadd ST(0), ST(1)				; X[1] = E[1] + O[1]
	fmul ST(0), ST(3)				; X[1] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[1])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 1], R8b		;
	fld ST(1)						; ST(0) = E[1]
	fsub ST(0), ST(1)				; X[5] = E[1] - O[1]
	fmul ST(0), ST(3)				; X[5] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[5])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 5], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 3 * 4 + 4]	; Загрузка E[2]
	fld dword ptr [RSP + 7 * 4 + 4]	; Загрузка O[2]
	fld ST(1)						; ST(0) = E[2]
	fadd ST(0), ST(1)				; X[2] = E[2] + O[2]
	fmul ST(0), ST(3)				; X[2] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[2])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 2], R8b		;
	fld ST(1)						; ST(0) = E[2]
	fsub ST(0), ST(1)				; X[6] = E[2] - O[2]
	fmul ST(0), ST(3)				; X[6] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[6])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 6], R8b		;
	fstp ST(0)
	fstp ST(0)

	fld dword ptr [RSP + 4 * 4 + 4]	; Загрузка E[3]
	fld dword ptr [RSP + 8 * 4 + 4]	; Загрузка O[3]
	fld ST(1)						; ST(0) = E[3]
	fadd ST(0), ST(1)				; X[3] = E[3] + O[3]
	fmul ST(0), ST(3)				; X[3] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[3])
	mov R8b, byte ptr [RSP]			;
	mov byte ptr [RCX + 3], R8b		;
	fld ST(1)						; ST(0) = E[3]
	fsub ST(0), ST(1)				; X[7] = E[3] - O[3]
	fmul ST(0), ST(3)				; X[7] * (1 / N)
	fistp dword ptr [RSP]			; Запись Re(X[7])
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

	AVX_SUPPORTED:
	pop RDX
	pop RCX
	vzeroall;
	vmovups xmm0, xmmword ptr [RDX]					; Получение нужных
	vmovups xmm1, xmmword ptr [RDX + 4 * 4]			;	данных из памяти
	vmovups xmm12, xmmword ptr [RDX + 4 * 8]		;	и их запись
	vmovups xmm14, xmmword ptr [RDX + 4 * (8 + 4)]	;	в необходимом
	vpermilps ymm0, ymm0, 00011011b					;	для корректной работы
	vpermilps ymm1, ymm1, 00011011b					;	функции виде
	vpermilps ymm12, ymm12, 00011011b				;	в регистры
	vpermilps ymm14, ymm14, 00011011b				;	ymm(0, 1, 12, 14)

	push RCX
	push RDX
	mov RCX, RSP
	add RCX, 4 * 8 * 2
	vmovups ymm13, ymmword ptr [XImMV]	; Запись соответствующего
	vmovups ymm15, ymmword ptr [ImMV]	;	множителя (Spectrum -> Signal)
	call FourierTransform				; Вызов функции для вычисления X[0-7]
	vmovups ymm0, [RCX]					; Запись результата в регистр ymm0
	pop RDX
	pop RCX

	vmovups ymm2, [const18] ; Деление результата на 1/8
	vmulps ymm0, ymm0, ymm2
	vmulps ymm1, ymm1, ymm2
	vcvtps2dq ymm0, ymm0

	; vpmovdb xmm0, ymm0- AVX512
	; Запись результата в виде массива int8
	vmovups ymmword ptr [RSP], ymm0
	mov r8b, byte ptr [RSP]
	mov byte ptr [RCX], r8b;
	mov r8b, byte ptr [RSP + 4 * 1]
	mov byte ptr [RCX + 1], r8b;
	mov r8b, byte ptr [RSP + 4 * 2]
	mov byte ptr [RCX + 2], r8b;
	mov r8b, byte ptr [RSP + 4 * 3]
	mov byte ptr [RCX + 3], r8b;
	mov r8b, byte ptr [RSP + 4 * 4]
	mov byte ptr [RCX + 4], r8b;
	mov r8b, byte ptr [RSP + 4 * 5]
	mov byte ptr [RCX + 5], r8b;
	mov r8b, byte ptr [RSP + 4 * 6]
	mov byte ptr [RCX + 6], r8b;
	mov r8b, byte ptr [RSP + 4 * 7]
	mov byte ptr [RCX + 7], r8b;

	add RSP, 4 * 8 * 7

	ret
RecoverSignal ENDP
END
