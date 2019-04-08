	DEF_ASM_OP0(nop, SVMOPCODE_NOP) /* must be first OP0 */
	DEF_ASM_OP0(dbgbreak, SVMOPCODE_DBGBREAK)
	DEF_ASM_OP0(ret, SVMOPCODE_RET)
	DEF_ASM_OP0(leave, SVMOPCODE_LEAVE)
	DEF_ASM_OP0(vmexit, SVMOPCODE_VMEXIT)

ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_IM8S_R,	OPC_B, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_IM16S_R,	OPC_W, OPT_IM16S, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_IMM_R,		OPC_L, OPT_IM32, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_RIP_R,		0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_RBP8_R,	OPC_DISP8S, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_RBP_R,		0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_IND_R,		0, OPT_EA | OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_IMM_R,		OPC_Q, OPT_IM64, OPT_REG))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_RIP_R,		0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_RBP8_R,	OPC_DISP8S, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_RBP_R,		0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_IND_R,		0, OPT_EA | OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(movb, SVMOPCODE_MOVB_R_RIP,		0, OPT_REG, OPT_EA | OPT_RIP))
ALT(DEF_ASM_OP2(movb, SVMOPCODE_MOVB_R_RBP8,	OPC_DISP8S, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movb, SVMOPCODE_MOVB_R_RBP,		0, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movb, SVMOPCODE_MOVB_R_IND,		0, OPT_REG, OPT_EA | OPT_REG))
ALT(DEF_ASM_OP2(movsbl, SVMOPCODE_MOVSBL_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movsbl, SVMOPCODE_MOVSBL_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movsbl, SVMOPCODE_MOVSBL_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movsbl, SVMOPCODE_MOVSBL_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzbl, SVMOPCODE_MOVZBL_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzbl, SVMOPCODE_MOVZBL_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movzbl, SVMOPCODE_MOVZBL_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movzbl, SVMOPCODE_MOVZBL_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(movw, SVMOPCODE_MOVW_R_RIP,		0, OPT_REG, OPT_EA | OPT_RIP))
ALT(DEF_ASM_OP2(movw, SVMOPCODE_MOVW_R_RBP8,	OPC_DISP8S, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movw, SVMOPCODE_MOVW_R_RBP,		0, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movw, SVMOPCODE_MOVW_R_IND,		0, OPT_REG, OPT_EA | OPT_REG))
ALT(DEF_ASM_OP2(movswl, SVMOPCODE_MOVSWL_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movswl, SVMOPCODE_MOVSWL_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movswl, SVMOPCODE_MOVSWL_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movswl, SVMOPCODE_MOVSWL_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzwl, SVMOPCODE_MOVZWL_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzwl, SVMOPCODE_MOVZWL_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movzwl, SVMOPCODE_MOVZWL_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movzwl, SVMOPCODE_MOVZWL_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_R_RIP,		0, OPT_REG, OPT_EA | OPT_RIP))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_R_RBP8,	OPC_DISP8S, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_R_RBP,		0, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movl, SVMOPCODE_MOVL_R_IND,		0, OPT_REG, OPT_EA | OPT_REG))
ALT(DEF_ASM_OP2(movslq, SVMOPCODE_MOVSLQ_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movslq, SVMOPCODE_MOVSLQ_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movslq, SVMOPCODE_MOVSLQ_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movslq, SVMOPCODE_MOVSLQ_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzlq, SVMOPCODE_MOVZLQ_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(movzlq, SVMOPCODE_MOVZLQ_RIP_R,	0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(movzlq, SVMOPCODE_MOVZLQ_RBP_R,	0, OPT_EA | OPT_RBP, OPT_REG))
ALT(DEF_ASM_OP2(movzlq, SVMOPCODE_MOVZLQ_IND_R,	0, OPT_EA | OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_R_RIP,		0, OPT_REG, OPT_EA | OPT_RIP))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_R_RBP8,	OPC_DISP8S, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_R_RBP,		0, OPT_REG, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP2(movq, SVMOPCODE_MOVQ_R_IND,		0, OPT_REG, OPT_EA | OPT_REG))

ALT(DEF_ASM_OP2(lea, SVMOPCODE_LEA_RIP_R,		0, OPT_EA | OPT_RIP, OPT_REG))
ALT(DEF_ASM_OP2(lea, SVMOPCODE_LEA_RBP_R,		0, OPT_EA | OPT_RBP, OPT_REG))

ALT(DEF_ASM_OP1(setc, SVMOPCODE_SETC,			0, OPT_REG))
ALT(DEF_ASM_OP1(setb, SVMOPCODE_SETB,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnae, SVMOPCODE_SETNAE,		0, OPT_REG))
ALT(DEF_ASM_OP1(setnc, SVMOPCODE_SETNC,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnb, SVMOPCODE_SETNB,			0, OPT_REG))
ALT(DEF_ASM_OP1(setae, SVMOPCODE_SETAE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setz, SVMOPCODE_SETZ,			0, OPT_REG))
ALT(DEF_ASM_OP1(sete, SVMOPCODE_SETE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnz, SVMOPCODE_SETNZ,			0, OPT_REG))
ALT(DEF_ASM_OP1(setne, SVMOPCODE_SETNE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setbe, SVMOPCODE_SETBE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setna, SVMOPCODE_SETNA,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnbe, SVMOPCODE_SETNBE,		0, OPT_REG))
ALT(DEF_ASM_OP1(seta, SVMOPCODE_SETA,			0, OPT_REG))
ALT(DEF_ASM_OP1(sets, SVMOPCODE_SETS,			0, OPT_REG))
ALT(DEF_ASM_OP1(setns, SVMOPCODE_SETNS,			0, OPT_REG))
ALT(DEF_ASM_OP1(setl, SVMOPCODE_SETL,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnge, SVMOPCODE_SETNGE,		0, OPT_REG))
ALT(DEF_ASM_OP1(setnl, SVMOPCODE_SETNL,			0, OPT_REG))
ALT(DEF_ASM_OP1(setge, SVMOPCODE_SETGE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setng, SVMOPCODE_SETNG,			0, OPT_REG))
ALT(DEF_ASM_OP1(setle, SVMOPCODE_SETLE,			0, OPT_REG))
ALT(DEF_ASM_OP1(setg, SVMOPCODE_SETG,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnle, SVMOPCODE_SETNLE,		0, OPT_REG))
ALT(DEF_ASM_OP1(setp, SVMOPCODE_SETP,			0, OPT_REG))
ALT(DEF_ASM_OP1(setnp, SVMOPCODE_SETNP,			0, OPT_REG))
ALT(DEF_ASM_OP1(seto, SVMOPCODE_SETO,			0, OPT_REG))
ALT(DEF_ASM_OP1(setno, SVMOPCODE_SETNO,			0, OPT_REG))

ALT(DEF_ASM_OP1(ldsp, SVMOPCODE_LDSP_R,			0, OPT_REG))
ALT(DEF_ASM_OP1(ldsp, SVMOPCODE_LDSP_RBP,		0, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP1(ldsp, SVMOPCODE_LDSP_IND,		0, OPT_EA | OPT_REG))
ALT(DEF_ASM_OP1(stsp, SVMOPCODE_STSP_R,			0, OPT_REG))
ALT(DEF_ASM_OP1(stsp, SVMOPCODE_STSP_RBP,		0, OPT_EA | OPT_RBP))
ALT(DEF_ASM_OP1(stsp, SVMOPCODE_STSP_IND,		0, OPT_EA | OPT_REG))

ALT(DEF_ASM_OP1(ldbp, SVMOPCODE_LDBP_R,			0, OPT_REG))
ALT(DEF_ASM_OP1(stbp, SVMOPCODE_STBP_R,			0, OPT_REG))

ALT(DEF_ASM_OP2(add, SVMOPCODE_ADDSP_IM8S,		OPC_B, OPT_IM8S, OPT_RSP))
ALT(DEF_ASM_OP2(add, SVMOPCODE_ADDSP_IMM,		OPC_L, OPT_IM32, OPT_RSP))
ALT(DEF_ASM_OP2(sub, SVMOPCODE_SUBSP_R,			0, OPT_REG, OPT_RSP))
ALT(DEF_ASM_OP2(and, SVMOPCODE_ANDSP_IM8S,		OPC_B, OPT_IM8S, OPT_RSP))

ALT(DEF_ASM_OP1(enter, SVMOPCODE_ENTER,			OPC_L, OPT_IM32))
ALT(DEF_ASM_OP1(retn, SVMOPCODE_RETN,			OPC_W, OPT_IM16))

ALT(DEF_ASM_OP1(jc, SVMOPCODE_JC,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jb, SVMOPCODE_JB,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnae, SVMOPCODE_JNAE,			OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnc, SVMOPCODE_JNC,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnb, SVMOPCODE_JNB,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jae, SVMOPCODE_JAE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jz, SVMOPCODE_JZ,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(je, SVMOPCODE_JE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnz, SVMOPCODE_JNZ,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jne, SVMOPCODE_JNE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jbe, SVMOPCODE_JBE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jna, SVMOPCODE_JNA,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnbe, SVMOPCODE_JNBE,			OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(ja, SVMOPCODE_JA,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(js, SVMOPCODE_JS,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jns, SVMOPCODE_JNS,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jl, SVMOPCODE_JL,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnge, SVMOPCODE_JNGE,			OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnl, SVMOPCODE_JNL,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jge, SVMOPCODE_JGE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jng, SVMOPCODE_JNG,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jle, SVMOPCODE_JLE,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jg, SVMOPCODE_JG,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnle, SVMOPCODE_JNLE,			OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jo, SVMOPCODE_JO,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jno, SVMOPCODE_JNO,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jp, SVMOPCODE_JP,				OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jnp, SVMOPCODE_JNP,				OPC_JMP, OPT_ADDR))

ALT(DEF_ASM_OP1(call, SVMOPCODE_CALL_R,			0, OPT_REG))
ALT(DEF_ASM_OP1(call, SVMOPCODE_CALL,			OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(ncall0, SVMOPCODE_NCALL0_R,		0, OPT_REG))
ALT(DEF_ASM_OP2(ncall, SVMOPCODE_NCALL_IMM_R,	0, OPT_IM16, OPT_REG))
ALT(DEF_ASM_OP2(ntcall, SVMOPCODE_NTCALL_IMM_R,	0, OPT_IM16, OPT_REG))

ALT(DEF_ASM_OP1(pushl, SVMOPCODE_PUSHL_R,		0, OPT_REG))
ALT(DEF_ASM_OP1(pushq, SVMOPCODE_PUSHQ_R,		0, OPT_REG))
ALT(DEF_ASM_OP1(push32, SVMOPCODE_PUSH32_R,		0, OPT_REG))

ALT(DEF_ASM_OP2(imull, SVMOPCODE_IMULL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(imulq, SVMOPCODE_IMULQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(divl, SVMOPCODE_DIVL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(divq, SVMOPCODE_DIVQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(udivl, SVMOPCODE_UDIVL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(udivq, SVMOPCODE_UDIVQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(modl, SVMOPCODE_MODL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(modq, SVMOPCODE_MODQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(umodl, SVMOPCODE_UMODL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(umodq, SVMOPCODE_UMODQ_R_R,		0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(shll, SVMOPCODE_SHLL_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(shll, SVMOPCODE_SHLL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(shrl, SVMOPCODE_SHRL_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(shrl, SVMOPCODE_SHRL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(sarl, SVMOPCODE_SARL_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(sarl, SVMOPCODE_SARL_R_R,		0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(shlq, SVMOPCODE_SHLQ_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(shlq, SVMOPCODE_SHLQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(shrq, SVMOPCODE_SHRQ_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(shrq, SVMOPCODE_SHRQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(sarq, SVMOPCODE_SARQ_IMM_R,		0, OPT_IM8, OPT_REG))
ALT(DEF_ASM_OP2(sarq, SVMOPCODE_SARQ_R_R,		0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP1(jmp, SVMOPCODE_JMP_R,			0, OPT_REG))
ALT(DEF_ASM_OP1(jmp, SVMOPCODE_JMP_IND,			0, OPT_INDIR))
ALT(DEF_ASM_OP1(jmp, SVMOPCODE_JMPS,			OPC_SHORTJMP | OPC_JMP, OPT_ADDR))
ALT(DEF_ASM_OP1(jmps, SVMOPCODE_JMPS,			OPC_SHORTJMP, OPT_ADDR))

ALT(DEF_ASM_OP2(addl, SVMOPCODE_ADDL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(addl, SVMOPCODE_ADDL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(addl, SVMOPCODE_ADDL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(orl,  SVMOPCODE_ORL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(orl,  SVMOPCODE_ORL_IM8S_R,		OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(orl,  SVMOPCODE_ORL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(adcl, SVMOPCODE_ADCL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(adcl, SVMOPCODE_ADCL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(adcl, SVMOPCODE_ADCL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(sbbl, SVMOPCODE_SBBL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(sbbl, SVMOPCODE_SBBL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(sbbl, SVMOPCODE_SBBL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(andl, SVMOPCODE_ANDL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(andl, SVMOPCODE_ANDL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(andl, SVMOPCODE_ANDL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(subl, SVMOPCODE_SUBL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(subl, SVMOPCODE_SUBL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(subl, SVMOPCODE_SUBL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(xorl, SVMOPCODE_XORL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(xorl, SVMOPCODE_XORL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(xorl, SVMOPCODE_XORL_IMM_R,		OPC_L, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(cmpl, SVMOPCODE_CMPL_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(cmpl, SVMOPCODE_CMPL_IM8S_R,	OPC_L, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(cmpl, SVMOPCODE_CMPL_IMM_R,		OPC_L, OPT_IM, OPT_REG))

ALT(DEF_ASM_OP2(addq, SVMOPCODE_ADDQ_R_R,		OPC_Q, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(addq, SVMOPCODE_ADDQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(addq, SVMOPCODE_ADDQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(orq,  SVMOPCODE_ORQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(orq,  SVMOPCODE_ORQ_IM8S_R,		OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(orq,  SVMOPCODE_ORQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(adcq, SVMOPCODE_ADCQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(adcq, SVMOPCODE_ADCQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(adcq, SVMOPCODE_ADCQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(sbbq, SVMOPCODE_SBBQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(sbbq, SVMOPCODE_SBBQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(sbbq, SVMOPCODE_SBBQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(andq, SVMOPCODE_ANDQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(andq, SVMOPCODE_ANDQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(andq, SVMOPCODE_ANDQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(subq, SVMOPCODE_SUBQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(subq, SVMOPCODE_SUBQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(subq, SVMOPCODE_SUBQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(xorq, SVMOPCODE_XORQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(xorq, SVMOPCODE_XORQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(xorq, SVMOPCODE_XORQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))
ALT(DEF_ASM_OP2(cmpq, SVMOPCODE_CMPQ_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(cmpq, SVMOPCODE_CMPQ_IM8S_R,	OPC_Q, OPT_IM8S, OPT_REG))
ALT(DEF_ASM_OP2(cmpq, SVMOPCODE_CMPQ_IMM_R,		OPC_Q, OPT_IM, OPT_REG))

ALT(DEF_ASM_OP2(fcmps, SVMOPCODE_FCMPS_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fucmps, SVMOPCODE_FUCMPS_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fadds, SVMOPCODE_FADDS_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fsubs, SVMOPCODE_FSUBS_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fmuls, SVMOPCODE_FMULS_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fdivs, SVMOPCODE_FDIVS_R_R,		0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(fcmpd, SVMOPCODE_FCMPD_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fucmpd, SVMOPCODE_FUCMPD_R_R,	0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(faddd, SVMOPCODE_FADDD_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fsubd, SVMOPCODE_FSUBD_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fmuld, SVMOPCODE_FMULD_R_R,		0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fdivd, SVMOPCODE_FDIVD_R_R,		0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(fcvtl2s, SVMOPCODE_FCVTL2S_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtq2s, SVMOPCODE_FCVTQ2S_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvts2l, SVMOPCODE_FCVTS2L_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvts2q, SVMOPCODE_FCVTS2Q_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtuq2s, SVMOPCODE_FCVTUQ2S_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvts2uq, SVMOPCODE_FCVTS2UQ_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvts2d, SVMOPCODE_FCVTS2D_R_R, 0, OPT_REG, OPT_REG))

ALT(DEF_ASM_OP2(fcvtl2d, SVMOPCODE_FCVTL2D_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtq2d, SVMOPCODE_FCVTQ2D_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtd2l, SVMOPCODE_FCVTD2L_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtd2q, SVMOPCODE_FCVTD2Q_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtuq2d, SVMOPCODE_FCVTUQ2D_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtd2uq, SVMOPCODE_FCVTD2UQ_R_R, 0, OPT_REG, OPT_REG))
ALT(DEF_ASM_OP2(fcvtd2s, SVMOPCODE_FCVTD2S_R_R, 0, OPT_REG, OPT_REG))

#undef ALT
#undef DEF_ASM_OP0
#undef DEF_ASM_OP0L
#undef DEF_ASM_OP1
#undef DEF_ASM_OP2
#undef DEF_ASM_OP3
