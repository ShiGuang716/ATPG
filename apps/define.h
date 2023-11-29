// 逻辑门类型
#ifndef DEFINE_H_
#define DEFINE_H_
typedef unsigned char gtype;
#define 		AND 			0
#define 		NAND			1
#define 		OR  			2
#define 		NOR 			3
#define 		PI  			4 ///////////
#define 		XOR 			5
#define 		XNOR			6
#define 		DFF 			7
#define 		DUMMY   		8
#define 		BUFF			9
#define 		NOT 			10
#define			TIE_0			11
#define         MUX             12
#define         LATCH           13
#define         STEM            14
#define 		PO  			20 ///////////

typedef unsigned ftype;
#define         SA_0            0
#define         SA_1            1
#define 		OUTFAULT		7

// typedef struct FAULT
// {
// 	int gateId;	/* faulty gate */
// 	int faultLine;		/* faulty line, -1 if output fault */
// 	unsigned int  ftype;	/* fault type */
// 	unsigned int observe;
// } FAULT;

typedef struct MUXINFO
{
    long mux_dst_gate;
    long mux_s0_gate;
    long mux_a_gate;
    long mux_b_gate;
    int count;
} MUXINFO;

// static unsigned char add_SA0_Fault[4] = {0b00000001, 0b00000010, 0b00000100, 0b00001000};
// static unsigned char add_SA1_Fault[4] = {0b00010000, 0b00100000, 0b01000000, 0b10000000};
// static unsigned char del_SA0_Fault[4] = {0b11111110, 0b11111101, 0b11111011, 0b11110111};
// static unsigned char del_SA1_Fault[4] = {0b11101111, 0b11011111, 0b10111111, 0b01111111};

// static unsigned char add_SA_Fault[2][4] = {{0b00000001, 0b00000010, 0b00000100, 0b00001000},{0b00010000, 0b00100000, 0b01000000, 0b10000000}};
// static unsigned char del_SA_Fault[2][4] = {{0b11111110, 0b11111101, 0b11111011, 0b11110111},{0b11101111, 0b11011111, 0b10111111, 0b01111111}};

static unsigned short add_SA_Fault[2][8] = {{0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 0x0080},
											{0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000, 0x8000}};
static unsigned short del_SA_Fault[2][8] = {{0xfffe, 0xfffd, 0xfffb, 0xfff7, 0xffef, 0xffdf, 0xffbf, 0xff7f},
											{0xfeff, 0xfdff, 0xfbff, 0xf7ff, 0xefff, 0xdfff, 0xbfff, 0x7fff}};


// Yang Yun added
typedef struct STACK
{
	int last;
	int* list;
} STACK,*STACKPTR;

#define TEST_FOUND 0

#define 		UNDETECTED  	0
#define 		DETECTED		1
#define 		XDETECTED   	2   			/* Potentially detected */
#define 		REDUNDANT   	3
#define 	    PROCESSED	    4
/* typedef enum {LFREE,HEAD,BOUND} line_type; */
#define		LFREE		0
#define		HEAD		1
#define		BOUND		2


#define FREE    0
#define HEAD    1
#define BOUND   2

#define ZERO    0
#define ONE     1
#define X       2
#define D       3
#define DBAR    4

#define A_NOT(var)     ((var==ONE) ? ZERO :  \
			   (var==ZERO) ? ONE :  \
			   (var==D) ? DBAR : \
			   (var==DBAR) ? D : X)

#define copy(s1,s2) s2.last=s1.last;\
			  for(int i=s1.last;i>=0;i--) s2.list[i]=s1.list[i]

#define setline(obj, n0, n1) numzero[obj] = n0; \
                     numone[obj] = n1;

// #define setline(obj, n0, n1) obj.numzero = n0; \
//                      obj.numone = n1;					 
                    
#define delete(s,i) s.list[i]=s.list[(s.last)--]


// #define include_stem 

#define FAN
// #define ONEBYONE
// #define DEBUG 
//#define DEBUG_LEVEL_1
#ifdef DEBUG
#define Log(format, ...) printf("[%s][%s][%d]: " format "\n", __FILE__, \
								__FUNCTION__, __LINE__, ##__VA_ARGS__ ) 
#define Logerr(format, ...) printf("[%s][%s][%d] error: " format "\n", __FILE__, \
								__FUNCTION__, __LINE__, ##__VA_ARGS__ )
#else 
#define Log(format, ...)
#define Logerr(format, ...)
#endif


#define NUM  1

#define current_node g_tree.list[g_tree.last]

#define OVER_BACKTRACK  1
#define NO_TEST 2

/* Decision tree data structure */
typedef struct TREENODE
{
	long gate_id;		/* pointer to the gate */
	bool flag;		/* flag for backtracking */
	int pstack;		/* pointer for implication */
} TREETYPE, *TREEPTR;

struct ROOTTREE
{
	/* tree for backtracking */
	int last;
	TREEPTR list;
};

struct OBJ{
	long gate_id;
	unsigned int numzero;
	unsigned int numone;
};
#endif

/* mask bits for parallel pattern simulation */
#define BITSIZE	32		/* word size */
#define ALL0	0		/* (0,0,...,0) */
#define ALL1	(~0)		/* (1,1,...,1) */
#define MASK0	(~(ALL1<<1))	/* (1,1,...,1,0) */

#define setbit(x,y)  x|=(1<<y)
#define clrbit(x,y)  x&=~(1<<y)

static int BITMASK[32] =
{
	MASK0, MASK0 << 1, MASK0 << 2, MASK0 << 3, MASK0 << 4, MASK0 << 5, MASK0 << 6, MASK0 << 7, MASK0 << 8, MASK0 << 9, MASK0 << 10, MASK0 << 11, MASK0 << 12, MASK0 << 13, MASK0 << 14, MASK0 << 15, MASK0 << 16, MASK0 << 17, MASK0 << 18, MASK0 << 19, MASK0 << 20, MASK0 << 21, MASK0 << 22, MASK0 << 23, MASK0 << 24, MASK0 << 25, MASK0 << 26, MASK0 << 27, MASK0 << 28, MASK0 << 29, MASK0 << 30, MASK0 << 31
};

// #define 		ALL0			0
// #define 		ALL1			(~0)
// #define 		MASK0   		(~(ALL1<<1))