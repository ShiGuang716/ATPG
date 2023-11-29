int g_TruthTable1[21][5] = 
{
	/*  0		1		x		d		dbar		*/
	{	ZERO,	ONE,	X,		D,		DBAR	},	/* and */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* nand */
	{	ZERO,	ONE,	X,		D,		DBAR	},	/* or */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* nor */
	{	ZERO,	ONE,	X,		D,		DBAR	},	/* pi */
	{	ZERO,	ONE,	X,		D,		DBAR	},	/* xor */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* xnor */
	{	X,  	X,  	X,  	X,  	X   	},	/* dff */
	{	X,  	X,  	X,  	X,  	X,  	},	/* dummy */
	{	ZERO,   ONE,	X,  	D,  	DBAR	},	/* buffer */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* inverter */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 11 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 12 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 13 */
	{	ZERO,  	ONE,  	X,  	D,  	DBAR,  	},	/* 14 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 15 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 16 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 17 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 18 */
	{	X,  	X,  	X,  	X,  	X,  	},	/* 19 */
	{	ZERO,   ONE,	X,  	D,  	DBAR	},	/* po */
};

int g_TruthTable2[21][5][5] = 
{
	/*	0		1		x		d		dbar 	*/
	/* and */
	{{	ZERO,	ZERO,	ZERO,	ZERO,	ZERO	},	/* 0 */
	{	ZERO,	ONE,	X,		D,		DBAR	}, 	/* 1 */
	{	ZERO,	X,		X,		X,		X		}, 	/* x */
	{	ZERO,	D,		X,		D,		ZERO	}, 	/* d */
	{	ZERO,	DBAR,	X,		ZERO,	DBAR	}}, /* dbar */
	/* nand */
	{{	ONE,	ONE,	ONE,	ONE,	ONE		},	/* 0 */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* 1 */
	{	ONE,	X,		X,		X,		X		},	/* x */
	{	ONE,	DBAR,	X,		DBAR,	ONE		},	/* d */
	{	ONE,	D,		X,		ONE,	D		}},	/* dbar */
	/* or */
	{{	ZERO,	ONE,	X,		D,		DBAR	},	/* 0 */
	{	ONE,	ONE,	ONE,	ONE,	ONE		},	/* 1 */
	{	X,		ONE,	X,		X,		X		},	/* x */
	{	D,		ONE,	X,		D,		ONE		},	/* d */
	{	DBAR,	ONE,	X,		ONE,	DBAR	}},	/* dbar */
	/* nor */
	{{	ONE,	ZERO,	X,		DBAR,	D		},	/* 0 */
	{	ZERO,	ZERO,	ZERO,	ZERO,	ZERO	},	/* 1 */
	{	X,		ZERO,	X,		X,		X		},	/* x */
	{	DBAR,	ZERO,	X,		DBAR,	ZERO	},	/* d */
	{	D,		ZERO,	X,		ZERO,	D		}},	/* dbar */
	/* pi */
	{{	X,		X,		X,		X,		X		},	/* 0 */
	{	X,		X,		X,		X,		X		},	/* 1 */
	{	X,		X,		X,		X,		X		},	/* x */
	{	X,		X,		X,		X,		X		},	/* d */
	{	X,		X,		X,		X,		X		}},	/* dbar */
	/*xor  */
	{{	ZERO,	ONE,	X,		D,		DBAR	},	/* 0 */
	{	ONE,	ZERO,	X,		DBAR,	D		},	/* 1 */
	{	X,		X,		X,		X,		X		},	/* x */
	{	D,		DBAR,	X,		ZERO,	ONE		},	/* d */
	{	DBAR,	D,		X,		ONE,	ZERO	}},	/* dbar */
	/*xnor */
	{{	ONE,	ZERO,	X,		DBAR,	D		},	/* 0 */
	{	ZERO,	ONE,	X,		D,		DBAR	},	/* 1 */
	{	X,		X,		X,		X,		X		},	/* x */
	{	DBAR,	D,		X,		ONE,	ZERO	},	/* d */
	{	D,		DBAR,	X,		ZERO,	ONE		}},	/* dbar */

	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* dff */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* dummy */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* buffer */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* not */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 11 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 12 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 13 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 14 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 15 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 16 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 17 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 18 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* 19 */
	{{0, }, {0, }, {0, }, {0, }, {0,}}, 				/* po */
};