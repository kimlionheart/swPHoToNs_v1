	.set noreorder
	.set noat
	.arch sw3
	#  /usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-500/be::5.421-sw-500

	#-----------------------------------------------------------
	# Compiling test_cpe.c (/tmp/ccI#.XlR7a4)
	#-----------------------------------------------------------

	#-----------------------------------------------------------
	# Options:
	#-----------------------------------------------------------
	#  Target:SW3, ISA:ISA_1, Endian:little, Pointer Size:64
	#  -O2	(Optimization level)
	#  -g0	(Debug level)
	#  -m2	(Report advisories)
	#-----------------------------------------------------------

	.file	1	"/home/export/online1/swmore/opensource/lwpf/lwpf2_release/test_cpe.c"
	.file	2	"/usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-500/include/simd.h"
	.file	3	"/home/export/online1/swmore/opensource/lwpf/lwpf2_release/lwpf2.h"

	.comm	lwpf_global_counter_TEST, 6144, 32

	.section .tdata, "waT", "progbits"
	.align	0

	.section .tdata_local_fix, "waT", "progbits"
	.align	5

	.section .tdata_local, "waT", "progbits"
	.align	5

	.section .text1, "ax", "progbits"
	.align	4

	.section .data, "wa", "progbits"
	.align	5
	.section .text1

	# Program Unit: lwpf_init_TEST
	.align 4
	.ent	lwpf_init_TEST#
	.weak	lwpf_init_TEST
lwpf_init_TEST:	# 0x0
	# anon14 = 64
	# return_address = 32
	.loc	1	98	0
#  94      lwpf_local_counter[kernel] = simd_vaddl(lwpf_local_counter[kernel], cntrs); \
#  95    }									\
#  96  
#  97  #define U(x) lwpf_init_ ## x
#  98  LWPF_WEAK void LWPF_UNIT(perf_config_t *conf){
.BB1_lwpf_init_TEST:
	.prologue
	ldih	$gp,0($27)               	!gpdisp!1	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!1	# [1] 0
	
$ng..lwpf_init_TEST:
	ldi	$sp,-96($sp)              	# [3] 
.LCFI_lwpf_init_TEST_adjustsp:
	stl	$26,32($sp)               	# [4] return_address
.LCFI_lwpf_init_TEST_store26:
	.loc	1	13	0
	ldl	$17,16($16)               	# [4] id:45
	.loc	1	11	0
	.frame $15,96,$26,0
	.mask 0x4008000,-96
	ldl	$8,0($16)                 	# [5] id:43
	.loc	1	12	0
	ldl	$19,8($16)                	# [5] id:44
	.loc	1	100	0
#  99    int256 lwpf_local_counter[LWPF_KERNELS_END];
# 100    set_perf_mode(conf->pcr0, conf->pcr1, conf->pcr2, conf->pcrc);
	ldl	$18,24($16)               	# [6] id:42
	.loc	1	13	0
	sll	$17,59,$17                	# [7] 
	.loc	1	11	0
	sll	$8,59,$8                  	# [8] 
	.loc	1	12	0
	sll	$19,59,$16                	# [8] 
	wcsr $8, 5
wcsr $16, 6
wcsr $17, 7
wcsr $18, 8

	.loc	1	57	0
	ldl	$7,lwpf_global_counter_TEST($gp)	!literal	# [11] lwpf_global_counter_TEST
	ldih	$5,lwpf_local_counter($31)	!tprelhi	# [11] lwpf_local_counter
	.globl	_MYID
	ldih	$0,_MYID($31)            	!tprelhi	# [11] _MYID
.BB2_lwpf_init_TEST:
	ldw	$18,_MYID($0)             	!tprello	# [0] _MYID
	.loc	1	56	0
	stw	$31,64($sp)               	# [0] anon14
	.loc	1	57	0
	mov	$31,$16                   	# [0] 
	mov	$31,$21                   	# [0] 
	ldi	$17,lwpf_local_counter($5)	!tprello	# [1] lwpf_local_counter
	mov	96,$19                    	# [1] 
	ldi	$20,64($sp)               	# [1] anon14
	sll	$18,5,$18                 	# [3] 
	stl	$31,0($sp)                	# [4] id:48
	s4subl	$18,$18,$18            	# [4] 
	addl	$7,$18,$18               	# [5] 
	.globl	athread_put
	bsr	$26,athread_put           	# [5] athread_put
.BB3_lwpf_init_TEST:
	ldw	$0,64($sp)                	# [0] anon14
	ldih	$gp,0($26)               	!gpdisp!2	# [0] 0
	ldi	$gp,0+4($gp)              	!gpdisp!2	# [1] 0
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.Lt_0_3842             	# [4] 
.BB8_lwpf_init_TEST:
.BB6_lwpf_init_TEST:
#<loop> Loop body line 57
#<loop> unrolled 3 times
	ldw	$0,64($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB7_lwpf_init_TEST    	# [4] 
.BB9_lwpf_init_TEST:
#<loop> Part of loop body line 196608, head labeled .BB8_lwpf_init_TEST
#<loop> unrolled 3 times
	ldw	$0,64($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB7_lwpf_init_TEST    	# [4] 
.Lt_0_4354:
#<loop> Part of loop body line 196608, head labeled .BB8_lwpf_init_TEST
	ldw	$0,64($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	beq	$0,.BB8_lwpf_init_TEST    	# [4] 
.Lt_0_3842:
.BB7_lwpf_init_TEST:
	.loc	1	106	0
# 102    int i;
# 103    for (i = 0; i < LWPF_KERNELS_END; i ++){
# 104      lwpf_local_counter[i] = v0;
# 105    }
# 106    lwpf_sync_counters_c2m();
	ldl	$26,32($sp)               	# [0] return_address
	ldi	$sp,96($sp)               	# [0] 
	ret	$31,($26),1               	# [3] 
.L_CC_lwpf_init_TEST:
# PU cycle count: 0.000000
	.end	lwpf_init_TEST

	.section .tdata_local
	.org 0x0
	.align	0
	.weak	lwpf_local_counter
	.type	lwpf_local_counter, @object
	.size	lwpf_local_counter, 96
lwpf_local_counter:	# 0x0
	.skip 96
	# end of initialization for lwpf_local_counter

	.section .data
	.org 0x0
	.align	0
	.weak	lwpf_kernel_names_TEST
	.type	lwpf_kernel_names_TEST, @object
	.size	lwpf_kernel_names_TEST, 32
lwpf_kernel_names_TEST:	# 0x0
	.quad	.rodata +0
	.quad	.rodata +2
	.quad	.rodata +4
	.quad	.rodata +8
	# end of initialization for lwpf_kernel_names_TEST
	.org 0x20
	.align	0
	.weak	lwpf_kernel_count_TEST
	.type	lwpf_kernel_count_TEST, @object
	.size	lwpf_kernel_count_TEST, 8
lwpf_kernel_count_TEST:	# 0x20
	# offset 32
	.quad	3
	# end of initialization for lwpf_kernel_count_TEST

	.section .rodata, "a", "progbits"
	.align	3
	.comm	lwpf_global_counter_TEST, 6144, 32
	.section .text1
	.align 4

	# Program Unit: test
	.align 4
	.ent	test#
	.globl	test
test:	# 0xb0
	# anon19 = 72
	# anon14 = 68
	# anon13 = 64
	# return_address = 32
	# _temp_gra_spill0 = 80
	# _temp_gra_spill1 = 88
	.loc	1	6	0
#   2  
#   3  #define LWPF_KERNELS K(A) K(B) K(C)
#   4  #define LWPF_UNIT U(TEST)
#   5  #include "lwpf2.h"
#   6  int test(){
.BB1_test:
	.prologue
	ldih	$gp,0($27)               	!gpdisp!3	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!3	# [1] 0
	
$ng..test:
	ldi	$sp,-96($sp)              	# [3] 
.LCFI_test_adjustsp:
	ldih	$28,_MYID($31)           	!tprelhi	# [3] _MYID
	ldih	$18,lwpf_local_counter($31)	!tprelhi	# [3] lwpf_local_counter
	.loc	1	50	0
	.frame $15,96,$26,0
	.mask 0x4008000,-96
	mov	96,$19                    	# [3] 
	.loc	1	6	0
	stl	$26,32($sp)               	# [4] return_address
.LCFI_test_store26:
	.loc	1	50	0
	ldw	$27,_MYID($28)            	!tprello	# [4] _MYID
	.loc	1	6	0
	ldi	$18,lwpf_local_counter($18)	!tprello	# [4] lwpf_local_counter
	ldi	$28,_MYID($28)            	!tprello	# [4] _MYID
	ldl	$17,lwpf_global_counter_TEST($gp)	!literal	# [5] lwpf_global_counter_TEST
	.loc	1	49	0
	stw	$31,64($sp)               	# [5] anon13
	.loc	1	50	0
	mov	$31,$21                   	# [5] 
	mov	$31,$16                   	# [5] 
	.loc	1	6	0
	stl	$18,88($sp)               	# [6] _temp_gra_spill1
	stl	$28,80($sp)               	# [6] _temp_gra_spill0
	.loc	1	50	0
	ldi	$20,64($sp)               	# [6] anon13
	sll	$27,5,$27                 	# [7] 
	s4subl	$27,$27,$27            	# [8] 
	stl	$31,0($sp)                	# [8] id:31
	stl	$31,8($sp)                	# [9] id:32
	addl	$17,$27,$17              	# [9] 
	.globl	athread_get
	bsr	$26,athread_get           	# [9] athread_get
.BB2_test:
	ldw	$0,64($sp)                	# [0] anon13
	ldih	$gp,0($26)               	!gpdisp!4	# [0] 0
	ldi	$gp,0+4($gp)              	!gpdisp!4	# [1] 0
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.Lt_1_4610             	# [4] 
.BB17_test:
.BB15_test:
#<loop> Loop body line 50
#<loop> unrolled 3 times
	ldw	$0,64($sp)                	# [0] anon13
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB16_test             	# [4] 
.BB18_test:
#<loop> Part of loop body line 196608, head labeled .BB17_test
#<loop> unrolled 3 times
	ldw	$0,64($sp)                	# [0] anon13
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB16_test             	# [4] 
.Lt_1_5122:
#<loop> Part of loop body line 196608, head labeled .BB17_test
	ldw	$0,64($sp)                	# [0] anon13
	cmpeq	$0,1,$0                 	# [3] 
	beq	$0,.BB17_test             	# [4] 
.Lt_1_4610:
.BB16_test:
	.loc	1	8	0
#   7    lwpf_enter();
#   8    lwpf_start(A);
	ldih	$5,lwpf_local_counter($31) 	!tprelhi	# [0] lwpf_local_counter
	vldd	$5,lwpf_local_counter($5) 	!tprello	# [1] lwpf_local_counter
	rcsr $1, 7
	vshff $1, $1, 0x24, $1
	rcsr $1, 6
	vshff $1, $1, 0xc4, $1
	rcsr $1, 5
	vshff $1, $1, 0xe0, $1
	rcsr $1, 4
	vsubl $5, $1, $5
	.loc	1	9	0
#   9    if (_MYID == 0){
	ldl	$0,80($sp)                	# [7] _temp_gra_spill0
	ldw	$0,0($0)                  	# [10] id:28 _MYID+0x0
.BB5_test:
	.loc	1	8	0
	ldl	$1,88($sp)                	# [0] _temp_gra_spill1
	vstd	$5,0($1)                 	# [3] id:27 lwpf_local_counter+0x0
	.loc	1	9	0
	beq	$0,.BB6_test              	# [3] 
.Lt_1_7938:
	rcsr $6, 7
	vshff $6, $6, 0x24, $6
	rcsr $6, 6
	vshff $6, $6, 0xc4, $6
	rcsr $6, 5
	vshff $6, $6, 0xe0, $6
	rcsr $6, 4
	
	.loc	1	14	0
#  10      int i;
#  11      for (i = 0; i < 10; i ++)
#  12        printf("%p\n", lwpf_local_counter);
#  13    }
#  14    lwpf_stop(A);
	vaddl	$6,$5,$6                	# [1] 
	.loc	1	57	0
	sextw	$0,$7                   	# [1] 
.BB11_test:
	.loc	1	14	0
	ldl	$17,88($sp)               	# [0] _temp_gra_spill1
	.loc	1	57	0
	ldl	$18,lwpf_global_counter_TEST($gp)	!literal	# [0] lwpf_global_counter_TEST
	sll	$7,5,$27                  	# [0] 
	mov	$31,$16                   	# [0] 
	.loc	1	56	0
	stw	$31,68($sp)               	# [1] anon14
	.loc	1	57	0
	s4subl	$27,$27,$27            	# [1] 
	mov	$31,$21                   	# [1] 
	mov	96,$19                    	# [1] 
	ldi	$20,68($sp)               	# [2] anon14
	.loc	1	14	0
	vstd	$6,0($17)                	# [3] id:27 lwpf_local_counter+0x0
	.loc	1	57	0
	addl	$18,$27,$18              	# [3] 
	stl	$31,0($sp)                	# [4] id:33
	bsr	$26,athread_put           	# [4] athread_put
.BB12_test:
	ldw	$0,68($sp)                	# [0] anon14
	ldih	$gp,0($26)               	!gpdisp!5	# [0] 0
	ldi	$gp,0+4($gp)              	!gpdisp!5	# [1] 0
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.Lt_1_6914             	# [4] 
.BB23_test:
.BB21_test:
#<loop> Loop body line 57
#<loop> unrolled 3 times
	ldw	$0,68($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB22_test             	# [4] 
.BB24_test:
#<loop> Part of loop body line 196608, head labeled .BB23_test
#<loop> unrolled 3 times
	ldw	$0,68($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB22_test             	# [4] 
.Lt_1_7426:
#<loop> Part of loop body line 196608, head labeled .BB23_test
	ldw	$0,68($sp)                	# [0] anon14
	cmpeq	$0,1,$0                 	# [3] 
	beq	$0,.BB23_test             	# [4] 
.Lt_1_6914:
.BB22_test:
	.loc	1	15	0
#  15    lwpf_exit();
	ldl	$26,32($sp)               	# [0] return_address
	ldi	$sp,96($sp)               	# [0] 
	ret	$31,($26),1               	# [3] 
.BB6_test:
	.loc	1	11	0
	stw	$31,72($sp)               	# [0] anon19
	.align	4
.Lt_1_6402:
.BB19_test:
#<loop> Loop body line 11, nesting depth: 1, iterations: 10
	.loc	1	12	0
	ldl	$16,.rodata($gp)          	!literal	# [0] .rodata
	ldl	$17,88($sp)               	# [0] _temp_gra_spill1
	ldi	$16,32($16)               	# [3] 
	.globl	printf
	bsr	$26,printf                	# [3] printf
.BB8_test:
#<loop> Part of loop body line 65536, head labeled .Lt_1_6402
	.loc	1	11	0
	ldw	$0,72($sp)                	# [0] anon19
	.loc	1	12	0
#	.body
#	.label_state 1
	ldih	$gp,0($26)               	!gpdisp!6	# [0] 0
#	.body
#	.restore $sp
	ldi	$gp,0+4($gp)              	!gpdisp!6	# [1] 0
	.loc	1	11	0
	addw	$0,1,$0                  	# [3] 
#	.body
#	.copy_state 1
	stw	$0,72($sp)                	# [4] anon19
#	.body
#	.copy_state 1
	cmpeq	$0,10,$0                	# [4] 
	beq	$0,.Lt_1_6402             	# [5] 
.BB9_test:
.BB20_test:
	ldl	$5,88($sp)                	# [0] _temp_gra_spill1
	ldl	$0,80($sp)                	# [0] _temp_gra_spill0
#	.body
#	.copy_state 1
	vldd	$5,0($5)                 	# [3] id:27 lwpf_local_counter+0x0
	ldw	$0,0($0)                  	# [3] id:28 _MYID+0x0
	br	$31,.Lt_1_7938             	# [3] 
.L_CC_test:
# PU cycle count: 0.000000
	.end	test

	.section .rodata
	.org 0x0
	.align	0
	# offset 0
	.ascii "A\0"
	.org 0x2
	.align	0
	# offset 2
	.ascii "B\0"
	.org 0x4
	.align	0
	# offset 4
	.ascii "C\0"
	.org 0x8
	.align	0
	# offset 8
	.ascii "LWPF_KERNELS_END\0"
	.org 0x20
	.align	0
	# offset 32
	.byte	0x25, 0x70, 0xa, 0x0	# %p\n\000
	.weak	_tdata_local_start
	.weak	_tdata_local_end
	.weak	_tdata_private_start
	.weak	_tdata_private_end
	.weak	_tdata_local_fix_end
	.section .tdata
	.align	0
	.section .tdata_local_fix
	.align	5
	.section .tdata_local
	.align	5
	.section .text1
	.align	4
	.section .data
	.align	5
	.section .rodata
	.align	3
#	.gpvalue 0

	.section .debug_info, "", "progbits"
	.align	0
	.byte	0x88, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_abbrev
	.long	0x65740108, 0x635f7473, 0x632e6570, 0x65706f00
	.long	0x2043436e, 0x32342e35, 0x77732d31, 0x3030352d
	.byte	0x00, 0x01, 0x00
	.long	.debug_line
	.long	0x6c620302, 0x5f667077, 0x74696e69, 0x5345545f
	.byte	0x54, 0x00, 0x01, 0x01, 0x04, 0x92, 0x1e, 0xe0
	.byte	0x00
	.quad	.BB1_lwpf_init_TEST
	.quad	.L_CC_lwpf_init_TEST
	.long	0x0000006b, 0x63620303, 0x00666e6f, 0x7fa09103
	.long	0x06010400, 0x74736574, 0x04010100, 0x00e01e92
	.quad	.BB1_test
	.quad	.L_CC_test
	.byte	0x00, 0x00

	.section .debug_frame, "", "progbits"
	.align	0

	.section .debug_aranges, "", "progbits"
	.align	0
	.byte	0x2c, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.byte	0x08, 0x00, 0x00, 0x00, 0x00, 0x00
	.quad	.BB1_lwpf_init_TEST
	.quad	.L_CC_test - .BB1_lwpf_init_TEST
	.long	0x00000000, 0x00000000, 0x00000000, 0x00000000

	.section .debug_pubnames, "", "progbits"
	.align	0
	.byte	0x2a, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.long	0x0000008c, 0x00000031, 0x6670776c, 0x696e695f
	.long	0x45545f74, 0x6b005453, 0x74000000, 0x00747365
	.byte	0x00, 0x00, 0x00, 0x00

	.section .eh_frame, "a", "progbits"
	.align	0

.LEHCIE:
	.long	.LEHCIE_end - .LEHCIE_begin
.LEHCIE_begin:
	.long 0x0
	.byte	0x01, 0x00, 0x01, 0x78, 0x1a, 0x0c, 0x1e, 0x00
	.align 3
.LEHCIE_end:

	.section .debug_abbrev, "", "progbits"
	.align	0
	.long	0x03011101, 0x13082508, 0x100b420b, 0x02000006
	.long	0x0b3a012e, 0x08030b3b, 0x408b0c3f, 0x110a400c
	.long	0x01011201, 0x03000013, 0x0b3a0005, 0x08030b3b
	.long	0x00000a02, 0x3a002e04, 0x030b3b0b, 0x8b0c3f08
	.byte	0x40, 0x0c, 0x40, 0x0a, 0x11, 0x01, 0x12, 0x01
	.byte	0x00, 0x00, 0x00, 0x00
	.section	.note.GNU-stack,"",@progbits
	.ident	"#SWCC Version 5.421-sw-500 : test_cpe.c compiled with : -O2 -msimd "

