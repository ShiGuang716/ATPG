// #include <bits/stdc++.h>
// #include "define.h"
// #include <math.h>
#include <forward_list>
// #define	INFINITY	999999
#define	MAXOBJ		100000
#define	MAXTREE		40000



#define 	FORWARD		0
#define 	BACKWARD	1
#define 	BACKWARD	1
#define 	CONFLICT	2
#define SIZE_OF_FUT 	32 

#define set(var) (var=true)
#define reset(var) ((var)=false)
#define my_max(a,b)			( ((a) > (b)) ? (a) : (b))


#define inCnt(idx) (GA.V[idx].getInDegree())
#define outCnt(idx) (GA.V[idx].getOutDegree())
#define inNbr(i,j) (GA.V[i].getInNeighbor(j))
#define outNbr(i,j) (GA.V[i].getOutNeighbor(j))
#define inNbrs(i) (GA.V[i].getInNeighbors())
#define outNbrs(i) (GA.V[i].getOutNeighbors())

#define is_free(line) (ltype[line]==LFREE)
#define is_head(line) (ltype[line]==HEAD)
#define is_bound(line) (ltype[line]==BOUND)
#define is_fanout(line) (outCnt(line)>1)
#define is_reachable_from_fault(line) (freach[line])
#define is_conflict(line) (numzero[line]>0 && numone[line]>0)
#define is_justified(line) (changed[line])
#define is_unjustified(line) ((!changed[line])&&(FANoutput[line]!=X))
#define is_D_propagated(line) (FANoutput[line] ==D || FANoutput[line]==DBAR)
#define is_flagged(node) node.flag

#define MALLOC(_type,number) \
(_type *)malloc((unsigned)(sizeof(_type)*(number)))
#define ALLOCATE(pointer,_type,number) \
    pointer=MALLOC(_type,number)

#define newA(__E, __n) (__E *)malloc((__n) * sizeof(__E))
#define EMPTY (-1)
#define delete_last(s) --(s.last)
#define ptr_is_empty(s) (s->last<0)
#define stackclear(s) s.last=EMPTY
#define clearevent(depth) stackclear(g_pEventListStack[depth])
#define push(s,ele) s.list[++(s.last)]=ele
#define pushevent(gate) \
	if (!changed[gate]) \
	{ \
		push(g_pEventListStack[level[gate]], gate); \
		set(changed[gate]); \
	}
#define is_empty(s) (s.last<0)
#define pop(s) s.list[(s.last)--]
#define popevent(depth) pop(g_pEventListStack[depth])

#define schedule_output(pGate) \
	for (int mac_i = 0; mac_i < outCnt(pGate); mac_i++) \
	   pushevent(outNbr(pGate,mac_i))
#define schedule_input(gate,i) pushevent(inNbr(gate,i)) \
				 schedule_output(inNbr(gate,i))


#define outputs2EventList(pGate,iDFrontierCnt, pOutGate) \
	for (int i = 0; i <outCnt(pGate); i++) \
	{ \
		pOutGate = outNbr(pGate,i); \
        if (!changed[pOutGate]) \
		{ \
			push(g_pEventListStack[level[pOutGate]], pOutGate); \
			iDFrontierCnt++; \
			set(changed[pOutGate]); \
		} \
	}


#define setFreach(pGate, i) freach[pGate] = i
#define checkFreach(pGate, i) freach[pGate] == i
#define gate_eval1(g,v,f,i) \
	if(inCnt(g)==1) v=g_TruthTable1[type[g]][FANoutput[inNbr(g,0)]]; \
	else if(inCnt(g)==2) \
	   v=g_TruthTable2[type[g]][FANoutput[inNbr(g,0)]][FANoutput[inNbr(g,1)]];\
	else { \
	   f = (type[g]==NAND) ? AND : \
		   (type[g]==NOR) ? OR : type[g]; \
	   v=g_TruthTable2[f][FANoutput[inNbr(g,0)]][FANoutput[inNbr(g,1)]];\
	   for(i=2;i<inCnt(g);i++) \
		  v=g_TruthTable2[f][v][FANoutput[inNbr(g,i)]]; \
	   v=(type[g]==NAND||type[g]==NOR)? A_NOT(v) : v; \
	}


extern std::vector<long> index_head;
extern std::vector<long> index_PI;
extern std::vector<long> index_tie0;
extern std::vector<long> index_PO;
// extern std::vector<OBJ> curObj;
// extern std::vector<OBJ> initObj;
// extern std::vector<OBJ> fanObj;
// extern std::vector<OBJ> headObj;
// extern std::vector<OBJ> finalObj;
extern long n;
extern long noPI;
extern long noPO;
extern long noTIE_0;

extern int* ltype;
extern gtype *type;
extern long n;
extern int *changed;
extern int *freach;
extern int *freach1;
extern int *cobserve;
extern int *observe;
extern int *xpath;
extern unsigned int *level;
extern int * FANoutput;
extern int *numzero;
extern int *numone;
extern unsigned int* dpo;
extern int SELECTMODE;
extern int *cont0;
extern int *cont1;
extern unsigned short *fault_list;
extern int iNoStemGates;
extern int* pStemGates;
extern int* fosIndex;

int dy_id = INT_MAX;
// int *dyn_dom = new int[MAXGATE];
vector<int> dom_array;

extern long noRedundant_faults; 

typedef struct link
{
	int ngate;
	struct link *next;
} LINKTYPE;

extern LINKTYPE **u_path;

typedef struct FAULT
{
	int gateId;	/* faulty gate */
	int faultLine;		/* faulty line, -1 if output fault */
	unsigned int  ftype;	/* fault type */
	int observe;
} FAULT;


std::unordered_map<int,int> g_iHeadGateIndex;
/* static buffers for fan */
STACK g_unjustStack,		/* set of unjustified lines */
		  g_initObjStack,			/* set of initial objectives */
		  g_curObjStack,			/* set of current objectives */
		  g_fanObjStack,			/* set of fanout objectives */
		  g_headObjStack,			/* set of head objectives */
		  g_finalObjStack,		/* set of final objectives */
		  g_DfrontierStack,		/* set of Dfrotiers */
		  g_stack;			/* stack for backtracing */
struct ROOTTREE g_tree;
STACK *g_pEventListStack;	
STACK g_freeGatesStack,	/* fault free simulation */
		    g_faultyGatesStack,   	  /* list of faulty gates */
			g_evalGatesStack, 		  /* STEM_LIST to be simulated */
			g_activeStemStack;  	   /* list of active stems */
// STACK *levelStack;
std::vector<vector<int> > levelStack;
extern std::vector<vector<FAULT> > dfault; // save the faults that contribute for observe

int g_iMaxLevel = 0;
extern unsigned int max_level;
// char* ltype = NULL;
// gtype *type = NULL;


template <class vertex>
struct LTYPE
{
    graph<vertex> &GA;
    unsigned int current_level;
    std::mutex &mtx_ltype;
    LTYPE(graph<vertex> &_GA, unsigned int _current_level, std::mutex &_mtx_ltype) : GA(_GA), current_level(_current_level), mtx_ltype(_mtx_ltype) {}
    inline bool updateAtomic(long s, long d)
    {
        short inDegree = GA.V[d].getInDegree();
        int ltype_tmp = FREE;
        char oldV = 0;
        char newV = 0;
        bool done = false;
        for (int i = 0; i < inDegree; i++)
        {
            if (ltype[GA.V[d].getInNeighbor(i)] != FREE)
            {
                ltype_tmp = BOUND;
                break;
            }
        }
    #ifdef include_stem 
        if (ltype_tmp == FREE && GA.V[d].getOutDegree() == 1 && type[GA.V[d].getOutNeighbor(0)] == STEM)
        {
            ltype_tmp = HEAD;
        }
    #else
        if (ltype_tmp == FREE && GA.V[d].getOutDegree() != 1 ) 
		{
			ltype_tmp = HEAD;
		}
    #endif

        if (ltype_tmp == HEAD)
        {
            mtx_ltype.lock();
            index_head.emplace_back(d);
            mtx_ltype.unlock();
        }
        else if (ltype_tmp == BOUND)
        {
            for (int i = 0; i < inDegree; i++)
            {
                long pTempGate = GA.V[d].getInNeighbor(i);
                if (ltype[pTempGate] == FREE)
                {
                    writeUpdate(&ltype[pTempGate], (int)HEAD);
                    mtx_ltype.lock();
                    index_head.emplace_back(pTempGate);
                    mtx_ltype.unlock();
                }
            }
        }
        writeUpdate(&ltype[d], ltype_tmp);
        return true;
    }

    inline bool update(long s, long d)
    {
        return updateAtomic(s, d);
    }

    inline bool cond(long d)
    {
        if (level[d] == current_level)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

template <class vertex>
void initUniquePath(graph<vertex> &GA ,int iNoGate, int iMaxDPI) /* depth of event list, number of output */ //set_unique_path
{
	//OUTPUT:  freach
	int i, j;
	int pGateId, pInGateId, pOutGateId;
	LINKTYPE* pUPath;
	int iGateCount, k;

	for (i = iNoGate - 1; i >= 0; i--) //For every FANOUT
	{
		pGateId = i;
		if (u_path[pGateId] == NULL) //pGate->output <= 1
		{
			continue;
		}

		//pGate is a FANOUT has owns a pDominator !!!
		iGateCount = 0;
		setFreach(pGateId, i);

		for (j = 0; j< outCnt(pGateId); j++) //For every output of pGate
		{
			pOutGateId = outNbr(pGateId, j);
			if (!changed[pOutGateId])
			{
				push(g_pEventListStack[level[pOutGateId]], pOutGateId);
				set(changed[pOutGateId]);
				iGateCount++;
			}
		}
		
		for (j = level[pGateId] + 1; j < iMaxDPI; j++) //For every g_pEventListStack after index pGate->dpi
		{
			while (!is_empty(g_pEventListStack[j]))
			{
				pGateId = pop(g_pEventListStack[j]);
				reset(changed[pGateId]);
				setFreach(pGateId, i);
				iGateCount--;
				
				if (iGateCount == 0) //EXIT for iteration of pGate !!!
				{
					//pGate == g_net[i]->u_path->ngate !!!
					pUPath = u_path[i]; //g_net[i] == original pGate			
					for (k = 0; k < inCnt(pGateId); k++)
					{
						pInGateId = inNbr(pGateId, k); //pInGate ---> pGate
						if (checkFreach(pInGateId, i))
						{
							pUPath->next = (LINKTYPE *)malloc(sizeof(LINKTYPE));
							pUPath = pUPath->next;
							pUPath->ngate= pInGateId;
							pUPath->next = NULL;
						}
					}
					break;
				}
				
				for (k = 0; k< outCnt(pGateId); k++)
				{
					pOutGateId = outNbr(pGateId, k);
					if (!changed[pOutGateId])
					{
						push(g_pEventListStack[level[pOutGateId]], pOutGateId);
						set(changed[pOutGateId]);
						iGateCount++;
					}
				}
			}
		}
	}
}

template <class vertex>
int setCctParameters(graph<vertex> &GA ,int iNoGate, int iNoPI) //set_cct_parameters
{
	int iHeadCnt = 0;

	int current_level = 1;
	std::mutex mtx_ltype;
	bool *frontier_ltype = newA(bool, n);
	memset(frontier_ltype, (bool)0, sizeof(bool) * n);
	memset(ltype, 0, sizeof(char) * n);
	vertexSubset Frontier_ltype(n, iNoPI, frontier_ltype);
	unsigned int pTempGate = 0;
	for(int i = 0; i < index_PI.size(); i++)
	{
		pTempGate = index_PI[i];
		frontier_ltype[pTempGate] = (bool)1;
		ltype[pTempGate] = FREE;
    #ifdef include_stem 
		if (GA.V[pTempGate].getOutDegree() == 1 && type[GA.V[pTempGate].getOutNeighbor(0)] == STEM) 
		{
			ltype[pTempGate] = HEAD;
		}
    #else
        if (GA.V[pTempGate].getOutDegree() > 1 ) 
		{
			ltype[pTempGate] = HEAD;
		}
    #endif
			
	}
	while (!Frontier_ltype.isEmpty())
	{
		vertexSubset out = edgeMap(GA, Frontier_ltype, LTYPE<vertex>(GA, current_level, mtx_ltype));
		Frontier_ltype.del();
		Frontier_ltype = out;
		current_level++;
	}
	Frontier_ltype.del();

	for (int k = 0 ; k < iNoGate ; k++){
        if (is_head(k)) iHeadCnt++;
    }

	printf("iHeadCnt:%d\n",iHeadCnt);
	int j = iHeadCnt;
	for (int i = iNoGate - 1; i >= 0; i--)
	{
		if (is_head(i))
		{
			g_iHeadGateIndex[--j] = i;
			// printf("%d %d\n", j, i);
		}
	}

	/* alloacate space for sets (needed for the fan algorithm) */
	ALLOCATE(g_unjustStack.list, int, MAXOBJ);
	ALLOCATE(g_initObjStack.list, int, MAXOBJ);
	ALLOCATE(g_curObjStack.list, int, MAXOBJ);
	ALLOCATE(g_fanObjStack.list, int, MAXOBJ);
	ALLOCATE(g_headObjStack.list, int, iHeadCnt);
	ALLOCATE(g_finalObjStack.list, int, MAXOBJ);
	ALLOCATE(g_DfrontierStack.list, int, MAXOBJ);

	if (g_stack.list == NULL)
	{
		ALLOCATE(g_stack.list, int, iNoGate);
	}
	ALLOCATE(g_tree.list, TREETYPE, MAXTREE);

    return iHeadCnt;

}

bool allocateStacks(int iNoGate) //allocate_dynamic_buffers
{
	if ((g_faultyGatesStack.list = (int *)malloc((unsigned)(sizeof(int) * iNoGate))) == NULL)
	{
		return(false);
	}
	// if ((g_freeGatesStack.list = (int *)malloc((unsigned)(sizeof(int) * iNoGate))) == NULL)
	// {
	// 	return(FALSE);
	// }
	if ((g_evalGatesStack.list = (int *)malloc((unsigned)(sizeof(int) * iNoGate))) == NULL)
	{
		return(false);
	}
	if ((g_activeStemStack.list = (int *)malloc((unsigned)(sizeof(int) * iNoGate))) == NULL)
	{
		return(false);
	}
	return(true);
}

void allocateEventListStacks()
{
	int i;

    ALLOCATE(g_pEventListStack, STACK, g_iMaxLevel + 2);
	for (i = 0; i < g_iMaxLevel + 2; i++)
	{
		stackclear(g_pEventListStack[i]); //last = 0
	}
	
	for (i = 0; i < n; i++)
	{
		g_pEventListStack[level[i]].last++;
	}
	
	for (i = 0; i < g_iMaxLevel + 2; i++)
	{
		ALLOCATE(g_pEventListStack[i].list, int , g_pEventListStack[i].last + 1 + SIZE_OF_FUT);
		stackclear(g_pEventListStack[i]);
	}
}

template <class vertex>
void initNetAndFreach(graph<vertex> &GA,int iNoGate, int pFaultyGateId, int iMaxDPI) //init_net
{
	//INPUT:  pFaultyGateId
	//OUTPUT:  freach
	int i, j;
	int pGateId;

	/* clear changed ochange and set freach */
	for (i = 0; i < iNoGate; i++)
	{
		if (is_free(i))
		{
			set(changed[i]); //don't pay attention to LFREE gates
		}
		else
		{
			reset(changed[i]);
		}
		reset(freach[i]);
		FANoutput[i] = X;
		xpath[i] = 1;
	}

	/* clear all sets */
	for (i = 0; i < iMaxDPI; i++)
	{
		clearevent(i);
	}
	stackclear(g_DfrontierStack);
	stackclear(g_unjustStack);
	stackclear(g_initObjStack);
	stackclear(g_curObjStack);
	stackclear(g_fanObjStack);
	stackclear(g_headObjStack);
	stackclear(g_finalObjStack);
	stackclear(g_stack);
	stackclear(g_tree);

	/* flag all reachable gates from the faulty gate */
	pushevent(pFaultyGateId);
	for (i = level[pFaultyGateId]; i < iMaxDPI; i++)
		while (!is_empty(g_pEventListStack[i]))
		{
			pGateId = popevent(i);
			reset(changed[pGateId]); //changed == 0 ------------> Need to handle
			set(freach[pGateId]); //All outputs of pFaultyGateId --------------> freach = 1
			for (j = 0; j < outCnt(pGateId); j++)
			{
				pushevent(outNbr(pGateId,j));
			}
		}
}



template <class vertex>
struct DPO
{
    graph<vertex> &GA;
    DPO(graph<vertex> &_GA) : GA(_GA) {}
    inline bool updateAtomic(long s, long d)
    {
        bool done = false;
        unsigned int oldV, newV;
        while (!done)
        {
            oldV = dpo[d];
            newV = dpo[s] + 1;
            if (newV > oldV)
            {
                done = CAS(&dpo[d], oldV, newV);
            }
            else
            {
                return false;
            }
        }
        return true;
    }

    inline bool update(long s, long d)
    {
        return updateAtomic(s, d);
    }

    inline bool cond(long d)
    {
        return true;
    }
};

template <class vertex>  
long deleteRedundantFaults(graph<vertex> &GA, int iNoGate) //check_redundant_faults
{
	register int i, j;
	long noRedundant_faults;
	register int pFanoutGate;

	noRedundant_faults = 0;
	for (i = 0; i < iNoGate; i++)
	{
		if (outCnt(i) > 1) //for every FANOUT pFanoutGate !!!
		{
			pFanoutGate = i; //FANOUT
			for (j = 0; j< outCnt(pFanoutGate); j++)
			{
				//changed = 0
				(changed[outNbr(pFanoutGate, j)])++;
			}
			
			for (j = 0; j< outCnt(pFanoutGate); j++)
			{
				if (changed[outNbr(pFanoutGate, j)] > 1) //for every pFanoutGate->outList[j] && changed > 1
				{
					int temp = outNbr(pFanoutGate, j);
					if (fault_list[temp] != 0) {
						for (unsigned int j=0; j<8; j++) {
							if ((fault_list[temp] & add_SA_Fault[0][j]) != 0) {
								if (j < 7)
								{
									if (inNbr(temp,j) == pFanoutGate) //pFanoutGate is 2rd gate which faults!!
									{
										fault_list[temp] = fault_list[temp] & del_SA_Fault[0][j];
										noRedundant_faults++;
									}
								}
							}
						}
						for (unsigned int j=0; j<8; j++) {
							if ((fault_list[temp] & add_SA_Fault[1][j]) != 0) {
								if (j < 7)
								{
									if (inNbr(temp,j) == pFanoutGate) //pFanoutGate is 2rd gate which faults!!
									{
										fault_list[temp] = fault_list[temp] & del_SA_Fault[1][j];
										noRedundant_faults++;
									}
								}
							}
						}
					}
				}
				reset(changed[outNbr(pFanoutGate, j)]); //change to 0 instead of 1 !!
			}
		}
	}
	return(noRedundant_faults);
}

template<class vertex>
void initGateStackAndFreach(graph<vertex> &GA, int iNoGate, int iMaxLevelAdd2, int iNoPI) //init_simulation
{
	int i;
	int pGateId;
	int n = GA.n;

	/* clear changed ochange and set freach */
	for (i = 0; i < iNoGate; i++)
	{
		changed[i] = false; 
		freach[i] = false;
		cobserve[i] = ALL0;
		observe[i] = ALL0;
	}

	/* clear all sets */
	for (i = 0; i < iMaxLevelAdd2; i++)
	{
		stackclear(g_pEventListStack[i]);
	}
	stackclear(g_stack);
	stackclear(g_faultyGatesStack);
	stackclear(g_evalGatesStack);
	stackclear(g_activeStemStack);

	/* initialize all stacks */
	for (i = 0; i < iNoGate; i++)
	{
		if (fault_list[i] != 0) //Faulty
		{
			push(g_faultyGatesStack, i); //Core sentence
		}
		if (outCnt(i) != 1) //PO && FANOUT
		{
			push(g_evalGatesStack, i); //Core sentence
		}
	}

	/* schedule freach */
	for (i = 0; i <= g_faultyGatesStack.last; i++)
	{
		pGateId = g_faultyGatesStack.list[i];
		while (!freach[pGateId])
		{
			set(freach[pGateId]); //Core sentence
			if (outCnt(pGateId) == 1) //Spread freach within the FFR !!
			{
				pGateId = outNbr(pGateId, 0);
			}
		}
	}
	
}

/*  init fanout free region
	compute FFRindex for each gate
	for critical path tracing and dominator analysis
*/
template <class vertex>
void initStemGatesAndFOS(graph<vertex> &GA, int iNoGate, int* pStemGates, int iNoStemGates) //setfanoutstem
{
	//OUTPUT:  pStemGates & fos
	register int i, j;
	register int pGate;

	j = 0;
	for (i = 0; i < iNoGate; i++)
	{
		if (outCnt(i) != 1) //FANOUT or PO
		{
			pStemGates[j++] = i;
		}
	}

	stackclear(g_stack);
	//j == iNoStemGates
	for (i = 0; i < iNoStemGates; i++)
	{
		push(g_stack, pStemGates[i]);
		while (!is_empty(g_stack))
		{
			pGate = pop(g_stack);
			fosIndex[pGate] = pStemGates[i]; //Core sentence
			for (j = 0; j< inCnt(pGate); j++)
			{
				if (outCnt(inNbr(pGate, j)) == 1) //Spread fosIndex within the FFR !!
				{
					push(g_stack, inNbr(pGate, j));
				}
			}
		}
	}

	// for (i = 0; i < iNoGate; i++)
	// {
	// 	printf("pGate : %d, fosindex :%d\n", i, fosIndex[i]);
	// }
}

template <class vertex>
void preprocess(graph<vertex> &GA, int iNoGate, int iNoPI, int iMaxLevelAdd2) {
    
    ltype = newA(int, n);
    g_iMaxLevel = max_level;
    
    xpath = newA(int, n);
    memset(ltype, 0, sizeof(int) * n);
    FANoutput = newA(int, n); 
    memset(FANoutput,X,sizeof(int)*n);
	freach = newA(int, n);
    memset(freach, 0, sizeof(int)*n);
	freach1 = newA(int, n);
    memset(freach1, 0, sizeof(int)*n);
	changed = newA(int, n);
    memset(changed, false, sizeof(int)*n);
    u_path = newA(LINKTYPE* , n);
    // memset(u_path, NULL, sizeof(LINKTYPE*)*n);
    numone = newA(int, n);
    numzero = newA(int ,n);
    cont0 = newA(int, n);
    memset(cont0,0,sizeof(int)*n);
    cont1 = newA(int, n);
    memset(cont1,0,sizeof(int)*n);
	cobserve = newA(int, n);
	// memset(cobserve, 0, sizeof(int)*n);
	observe = newA(int, n);
	// memset(observe, 0, sizeof(int)*n);

	allocateStacks(iNoGate);
    ALLOCATE(g_stack.list, int, iNoGate + 10);
	stackclear(g_stack);
    allocateEventListStacks();
    Log("allocateEventListStacks end");
    setCctParameters(GA ,iNoGate, iNoPI);
    Log("setCctParameters end");
    
    fosIndex = newA(int, n);
	memset(fosIndex, 0, sizeof(int)*n);
	iNoStemGates = 0;
	for (int i = 0; i < iNoGate; i++)
	{
		if (is_fanout(i) || type[i] == PO) //FANOUT or PO
		{
			iNoStemGates++;
		}
	}
	printf("num of fanout + PO is %d\n",iNoStemGates);
	pStemGates = (int *) malloc((unsigned)(sizeof(int) * iNoStemGates));
	initStemGatesAndFOS(GA, iNoGate, pStemGates, iNoStemGates);
	dfault.resize(n);
	
	//compute dpo for each gate
    // bool *frontier_dpo = newA(bool, n);
    // memset(frontier_dpo, (bool)0, sizeof(bool) * n);
    // parallel_for(int i = 0; i < index_PO.size(); i++)
    // {
    //     frontier_dpo[index_PO[i]] = (bool)1;
    // }
    dpo = newA(unsigned int, n);
    memset(dpo, 0, sizeof(unsigned int) * n);
    // GA.transpose();
    // vertexSubset Frontier_dpo(n, noPO, frontier_dpo);
    // while (!Frontier_dpo.isEmpty())
    // {
    //     vertexSubset out = edgeMap(GA, Frontier_dpo, DPO<vertex>(GA));     
    //     Frontier_dpo.del();
    //     Frontier_dpo = out;
    // }
    // Frontier_dpo.del();
    // #ifdef DEBUG_LEVEL_1
    // Log("dpo:");
    // for(int i=0; i<n;i++) {
    //     printf("%d %d\n", i, dpo[i]);
    // }
    // #endif
    // GA.transpose();
    
    setGateTestability(GA, n);
    Log("setGateTestability end");
	
	initFanoutGateDominators(GA,iNoGate,iMaxLevelAdd2);
    Log("initFanoutGateDominators end");
    initUniquePath(GA ,iNoGate,iMaxLevelAdd2);
    Log("initUniquePath end");


    // levelStack.resize(max_level+2,vector<int>(0));
    // for(int i=0;i<n;i++) // for good simulation
    // {   
    //     // cout<< i<< ","<<level[i]<<endl;
    //     levelStack[level[i]].emplace_back(i);
    // }

	noRedundant_faults = deleteRedundantFaults(GA, n);
	printf("noRedundant_faults : %ld\n",noRedundant_faults);

	initGateStackAndFreach(GA, iNoGate, iMaxLevelAdd2, iNoPI);

	
	// check upath
	// for(int i = 0; i < n;i++){
	// 	if(u_path[i] != NULL){
	// 		printf("%d ",i);
	// 		for (LINKTYPE *pLink = u_path[i]; pLink != NULL; pLink = pLink->next){
	// 			printf("%d ",pLink->ngate);
	// 		}
	// 		printf("\n");
	// 	}
	// }
}

/*	propFault2Headline
	Propagates the fault into a head line.
	Changes the fault into equivalent fault in the head line.
*/
template <class vertex>
void propFault2Headline(graph<vertex> &GA,FAULT pFault) //faulty_line_is_free  
{
	//PRECONDITION:  is_free(pGate) == TRUE
	//OUTPUT:  pFault & pGate->output
	int i;
	int iValue;
	int pGateId;

	pGateId = pFault.gateId;
	
	if (pFault.faultLine != OUTFAULT)
	{
		//					  pFault -------> pGate
		pGateId = inNbr(pFault.gateId, pFault.faultLine); //pGate to upper side
		//pGate -------> pFault
	}
    // pGateId = inNbr(pFault.gateId, pFault.faultLine);

	FANoutput[pGateId] = (pFault.ftype == SA_0) ? D : DBAR;


	//Always: is_free(pGate) --------------> pGate->outCount == 1
	push(g_stack, outNbr(pGateId,0));

	while (!is_empty(g_stack))
	{
		pGateId = pop(g_stack);
		iValue = (type[pGateId] == AND || type[pGateId] == NAND) ? ONE : ZERO;
		for (i = 0; i< inCnt(pGateId); i++) //UP (Always:  pGate->inCount == 1)
		{
			//Execute only once for s9234tc.bench
            int temp_id = inNbr(pGateId,i);
			if (FANoutput[temp_id] == X)
			{
				FANoutput[temp_id]  = iValue;
			}
			else
			{
				FANoutput[pGateId] = g_TruthTable1[type[pGateId]][FANoutput[temp_id] ]; //Always:  i == 0
			}
            

		}	
		if (is_free(pGateId))
		{
			push(g_stack, outNbr(pGateId, 0)); //DOWN
		}
	}
	
	set(changed[pGateId]);
	push(g_stack, pGateId);
	schedule_output(pGateId);
	pFault.gateId = pGateId;
	pFault.faultLine = OUTFAULT; //pGate -------> pFault
	pFault.ftype = (FANoutput[pGateId] == D) ? SA_0 : SA_1; //Why this?

}

/*	refFaultyGateOutput
	box 1 of the FAN algorithm.
	Defines input and output values of the faulty gate.
	Returns the level of the highest level gate in which
	backward implication is required.
	If a CONFLICT condition (redundant faults) occurs,
	returns (-1).
*/
template <class vertex>
int refFaultyGateOutput(graph<vertex> &GA, FAULT pFault) //set_faulty_gate  
{
	//INPUT:  pFault
	//OUTPUT:  iLastDPI & g_stack & many gates's ouput
	int i, iLastDPI;
	int pGateId;
	int iOriginalLineValue, iLineValue;
	int pUpflowGateId;

	iLastDPI = 0;
	pGateId = pFault.gateId;
	Log("current gate :%d",pGateId);
	pUpflowGateId = (pFault.faultLine == OUTFAULT) ? pGateId : inNbr(pGateId,pFault.faultLine); //get the upflow gate  
	iLineValue = (pFault.ftype == SA_0) ? D : DBAR; //D: 1/0     DBAR: 0/1

	/* input stuck-at faults */
	if (pFault.faultLine != OUTFAULT) /////////////pUpflowGate -> pFault -> pGate
	{
		Log("input line pFault");
		/* input line pFault */
		FANoutput[pUpflowGateId] = (iLineValue == D) ? ONE : ZERO;  /* faulty line */ //get pUpflowGate's output!!
		push(g_stack, pUpflowGateId);
		//output change ------> push into g_stack !!
		switch (type[pGateId])
		{
		case AND:
		case NAND:
			set(changed[pGateId]);
			for (i = 0; i< inCnt(pGateId); i++)
			{
				if (i != pFault.faultLine) //For every other inputs
				{
					if (FANoutput[inNbr(pGateId,i)] == X)
					{
						FANoutput[inNbr(pGateId,i)] = ONE; //activate other inputs
						push(g_stack, inNbr(pGateId,i));
					}
					else if (FANoutput[inNbr(pGateId,i)] != ONE)
					{
						Logerr("Conflicted!");
                        return(-1); //conflict occurs!
					}
				}
			}
			FANoutput[pGateId] = (type[pGateId] == NAND) ? A_NOT(iLineValue) : iLineValue; //get pGate's output!!
			push(g_stack, pGateId);
			break;
		case OR:
		case NOR:
			set(changed[pGateId]);
			for (i = 0; i< inCnt(pGateId); i++)
			{
				if (i != pFault.faultLine) //For every other inputs
				{
					if (FANoutput[inNbr(pGateId,i)] == X)
					{
						FANoutput[inNbr(pGateId,i)] = ZERO; //activate other inputs
						push(g_stack, inNbr(pGateId,i));
					}
					else if (FANoutput[inNbr(pGateId,i)] != ZERO)
					{
						Logerr("Conflicted!");
                        return(-1); //conflict occurs!
					}
				}
			}
			FANoutput[pGateId] = (type[pGateId] == NOR) ? A_NOT(iLineValue) : iLineValue; //get pGate's output!!
			push(g_stack, pGateId);
			break;
		case NOT:
			FANoutput[pGateId] = A_NOT(iLineValue);
			push(g_stack, pGateId);
			break;
		case BUFF:
		case PO:
			FANoutput[pGateId] = iLineValue;
			push(g_stack, pGateId);
			break;
        // case STEM:
        //     changed[pGateId] = true;
        //     FANoutput[pGateId] = iLineValue;
		// 	push(g_stack, pGateId);
		case XOR:
		case XNOR:
			break; //where's code?? no support for XOR and XNOR?
		}

		/* schedule events */
		if (FANoutput[pGateId] != X)
		{
			schedule_output(pGateId); //lead to refFaultyGateOutput
		}
		for (i = 0; i< inCnt(pGateId); i++)
		{
			if (FANoutput[inNbr(pGateId,i)] != X) //Always != X
			{
				Log("iLastDPI : %d, level[inNbr(pGateId,i)] :%d, inNbr(pGateId,i) :%d",iLastDPI, level[inNbr(pGateId,i)], inNbr(pGateId,i));
				iLastDPI = max((int)level[inNbr(pGateId,i)], iLastDPI);
				Log("iLastDPI:%d",iLastDPI);
				schedule_input(pGateId, i); //lead to refFaultyGateOutput
			}
		}
	}
	else /////////////pGate -> pFault
	{
		/* output line pFault */
		//pGate->output = (iLineValue == D) ? ONE : ZERO;
		Log("out line fault");
		FANoutput[pGateId] = iLineValue;
		push(g_stack, pGateId);
		schedule_output(pGateId);
		Log("pGate %d inCount is : %d,ltype :%d",pGateId, inCnt(pGateId), ltype[pGateId]);
		if (is_head(pGateId)) //pGate->outCount >= 2 (FANOUT)
		{
			set(changed[pGateId]);
		}

		//pGate->inCount == 1
		else if (inCnt(pGateId) == 1) ////////////////////pGate->inList[0] ----> pGate ----> pFault 
		{
			set(changed[pGateId]);
			iOriginalLineValue = (iLineValue == D) ? ONE : ZERO;
			FANoutput[inNbr(pGateId,0)] = g_TruthTable1[type[pGateId]][iOriginalLineValue]; //get input gate's output
			push(g_stack, inNbr(pGateId,0));
			schedule_input(pGateId, 0);
			iLastDPI = level[inNbr(pGateId,0)];
		}
		else if ((iLineValue == D && (type[pGateId] == AND || type[pGateId] == NOR)) || //Why ???
			(iLineValue == DBAR && (type[pGateId] == NAND || type[pGateId] == OR))) //(pGate->inCount >= 2)
		{
			set(changed[pGateId]);
			iOriginalLineValue = (type[pGateId] == AND || type[pGateId] == NAND) ? ONE : ZERO;
			for (i = 0; i< inCnt(pGateId); i++)
			{
				if (FANoutput[inNbr(pGateId,i)] == X)
				{
					FANoutput[inNbr(pGateId,i)] = iOriginalLineValue;
					iLastDPI = max((int)level[inNbr(pGateId,i)], iLastDPI);
					Log("iLastDPI:%d",iLastDPI);
					push(g_stack, inNbr(pGateId,i));
					schedule_input(pGateId,i);
				}
			}		
		}
		else
		{
			pushevent(pGateId);
			iLastDPI = level[pGateId];
			Log("iLastDPI:%d",iLastDPI);
		}
	}
    Log("iLastDPI:%d",iLastDPI);
	return(iLastDPI);
}

template <class vertex>
int initFanoutGateDominators(graph<vertex> &GA, int g_iNoGate, int iMaxLevelAdd2) /* depth of event list, number of output */ //set_dominator
{
	//OUTPUT:  iNoDominator
	 int i, j;
	 int pGateId, pTempGateId;
	int pDomiGateId;
	int iGateCount;
	int iNoDominator = 0;

	for (i = g_iNoGate - 1; i >= 0; i--)
	{
    
		pGateId = i;
		if (outCnt(pGateId) <= 1)
		{
			u_path[pGateId] = NULL;
		}
        
		else //outCount >= 2 (FANOUT)
		{
			iGateCount = 0;
            
            //pTempGateId = outNbr(pGateId,j);       
			outputs2EventList(pGateId, iGateCount, pTempGateId); //pGate ----> g_pEventListStack & iGateCount

			//OUTPUT:  g_pEventListStack & iGateCount (changed)
           
			for (j = level[pGateId] + 1; j < iMaxLevelAdd2; j++)
			{
                
				while (!is_empty(g_pEventListStack[j])) //For every level stack
				{
					pDomiGateId = pop(g_pEventListStack[j]);
					reset(changed[pDomiGateId]);
					if (iGateCount <= 0)
					{
						continue;
					}
					iGateCount--;
					
					if (iGateCount == 0) //LAST ONE !!
					{

                        u_path[pGateId] = (LINKTYPE*)malloc(sizeof(LINKTYPE));
						u_path[pGateId]->ngate = pDomiGateId;
                        u_path[pGateId]->next = NULL;
                        // cout << u_path[pGateId]->ngate << "," << pGateId <<endl;
						iNoDominator++; //EXIT!!!

						break;
					}
					else //iGateCount != 0   NOT LAST ONE !!
					{
						//PO                    
						if (outCnt(pDomiGateId) == 0)
						{
							u_path[pGateId] = NULL; //blocked
							iGateCount = 0; //Empty g_pEventListStack[j] &&&&&&&& End while !!
						}
						
						//FANOUT
						else if (outCnt(pDomiGateId) > 1)
						{
							if (u_path[pDomiGateId] == NULL)
							{
								u_path[pGateId] = NULL; //blocked
								iGateCount = 0; //Empty g_pEventListStack[j] &&&&&&&& End while !!
							}
							else
							{   
                                // cout<< "pDomiGateId," << pDomiGateId << "," << type[pDomiGateId] <<endl;
								pTempGateId = u_path[pDomiGateId]->ngate; //next is u_path->ngate!!
								if (!changed[pTempGateId])
								{
									push(g_pEventListStack[level[pTempGateId]], pTempGateId);
									set(changed[pTempGateId]);
									iGateCount++;
								}
							}
						}
						
						//NORMAL GATE
						else //pDomiGate->outCount == 1
						{
							if (!changed[outNbr(pDomiGateId,0)])
							{
								pTempGateId = outNbr(pDomiGateId,0);
								push(g_pEventListStack[level[pTempGateId]], pTempGateId);
								set(changed[pTempGateId]);
								iGateCount++;
							}
						}
					}                
				}
			}
		}
	}
	return(iNoDominator);
}


template <class vertex>
int unique_sensitize(graph<vertex> &GA, long pGateId, long pFaultyGateId)
{
	int i;
	LINKTYPE *pLink;
	long pNextGateId;
	int iLastDPI;		/* the largest depth of sensitized lines */
	bool bFlag;
	int iValue;

	iLastDPI = (-1);
    // printf("pGateId:%d,pFaultyGateId:%d\n",pGateId, pFaultyGateId);
	/* sensitize the current pGate */
	if (pGateId != pFaultyGateId)
	{
		//pFaultyGate ----------> pGate
		//other Gates ---iValue---> pGate
		iValue = (type[pGateId] == AND || type[pGateId] == NAND) ? ONE :
			(type[pGateId] == OR || type[pGateId] == NOR) ? ZERO : X;

		// printf("iValue:%d\n",iValue);	
		if (iValue != X)
		{
			for (i = 0; i< inCnt(pGateId); i++)
            {
				int temp_gateId = inNbr(pGateId,i);
                if (FANoutput[temp_gateId] == X)
				{
					FANoutput[temp_gateId] = iValue;
					iLastDPI = max((int)level[temp_gateId], iLastDPI);
                    Log("iLastDPI:%d",iLastDPI);
					push(g_stack, temp_gateId);
					schedule_input(pGateId, i);
				}
            }
		}
	}

	/* sensitize pNextGate path */
	while (1)
	{
		// Get the outdegree of current gate
		Log("the outdegree of current gate(%ld):%d", pGateId, outCnt(pGateId));
		if (outCnt(pGateId) == 0)
		{// If the outdegree is zero, the current gate is primary output. So stop the unique_sensitize
			break;
		}
		else if (outCnt(pGateId)== 1) //pGate -> pNextGate
		{// If the outdegree is one, we can get the next gate and propagate the D or DBAR value.
			pNextGateId = outNbr(pGateId, 0);
		}
		else if (u_path[pGateId] == NULL)
		{
			Log("u_path[%ld] == NULL",pGateId);
			break;
		}
		else
		{
			// Otherwise, the output of the current gate is fanout, we can not continue propagate the D or DBAR value.
            // NOTICE: In fact, we can find the next gate that the D or DBAR value must propagate. But now we do not consider it.
            Log("The output of the gate(%ld) is fanout.", pGateId);
			pNextGateId = u_path[pGateId]->ngate;
		}
		//pGate -------------------> pNextGate
		//pNextGate->inList[i] ---iValue---> pNextGate
		iValue = (type[pNextGateId] == AND ||type[pNextGateId] == NAND) ? ONE :
			(type[pNextGateId] == OR || type[pNextGateId] == NOR) ? ZERO : X;
		
		Log("pNextGateId : %ld, type : %d, iValue : %d",pNextGateId, type[pNextGateId],iValue);

		if (iValue != X) //////////if iValue == X, then pNextGate can't have multiple outputs!!!
		{
			if (outCnt(pGateId) == 1)
			{
				/* one fanout */
				for (i = 0; i< inCnt(pNextGateId); i++)
                {
					int temp_gateId = inNbr(pNextGateId, i);
                    if (temp_gateId != pGateId && FANoutput[temp_gateId] == X)
					{
						FANoutput[temp_gateId]= iValue;
						iLastDPI = max((int)level[temp_gateId], iLastDPI);						
                        Log("iLastDPI:%d",iLastDPI);
						push(g_stack, temp_gateId);
						schedule_input(pNextGateId, i);
					}
                }
			}
			else //pGate->outCount > 1
			{
				/* multiple fanout */
				for (i = 0; i< inCnt(pNextGateId); i++)
                {
					int temp_gateId = inNbr(pNextGateId, i);
                    if (FANoutput[temp_gateId] == X)
					{
						set(bFlag);
						for (pLink = u_path[pGateId]->next; pLink != NULL; pLink = pLink->next)
                        {
							if (temp_gateId == pLink->ngate)
							{
								reset(bFlag);
								break;
							}
                        }
						if (bFlag)
						{
							FANoutput[temp_gateId] = iValue;
							iLastDPI = max((int)level[temp_gateId], iLastDPI);
                            Log("iLastDPI:%d",iLastDPI);
							push(g_stack, temp_gateId);
							schedule_input(pNextGateId, i);
						}
					}

                }
			}
		}
		pGateId = pNextGateId;
	}
    Log("iLastDPI:%d",iLastDPI);
	return(iLastDPI);
}


template<class vertex>
int getFaultyGateOutput(graph<vertex> &GA,int pGateId, FAULT pFault) //faulty_gate_eval
{
	//Precondition:  pGate == pFault->gate
	//OUTPUT:  pGate->output
	int i, j;
	int ivalue;//value
	//level fval;
	int iGateType;

	if (inCnt(pGateId) == 0) //pGate -> pFault
	{
		return(FANoutput[pGateId]);
	}

	if (pFault.faultLine == OUTFAULT) //(pGate->inList[0]) -> pGate -> pFault
	{
		j = 0;
		ivalue = FANoutput[inNbr(pGateId,0)];
	}
	else //(pGate->inList[j]) -> pFault -> pGate
	{
		j = pFault.faultLine;
		ivalue = FANoutput[inNbr(pGateId,0)];
		if (ivalue == ZERO && pFault.ftype == SA_1)
		{
			ivalue = DBAR; //ivalue = pGate->input
		}
		else if (ivalue == ONE && pFault.ftype == SA_0)
		{
			ivalue = D; //ivalue = pGate->input
		}
	}

	if (inCnt(pGateId) == 1)
	{
		ivalue = g_TruthTable1[type[pGateId]][ivalue]; //ivalue = pGate->output
	}
	else //pGate->inCount > 1
	{
		iGateType = (type[pGateId] == NAND) ? AND : (type[pGateId] == NOR) ? OR : type[pGateId];
		for (i = 0; i < j; i++)
        {
			ivalue = g_TruthTable2[iGateType][ivalue][FANoutput[inNbr(pGateId,i)]];
        }
		for (++i; i< inCnt(pGateId); i++)
        {
			ivalue = g_TruthTable2[iGateType][ivalue][FANoutput[inNbr(pGateId,i)]];
        }
		if (type[pGateId] == NAND || type[pGateId] == NOR)
		{
			ivalue = A_NOT(ivalue);
		}
	}

	if (pFault.faultLine == OUTFAULT) //pGate -> pFault
	{
		if (ivalue == ZERO && pFault.ftype == SA_1)
		{
			ivalue = DBAR;
		}
		else if (ivalue == ONE && pFault.ftype == SA_0)
		{
			ivalue = D;
		}
	}

	return(ivalue);
}

template<class vertex>
int eval(graph<vertex> &GA, int pGateId, FAULT pFault)
{
	int i, j;
	int iValue, iTempValue;
	int iNumX;
	int iTempValue2;
	auto *pGateInList = inNbrs(pGateId);

	reset(changed[pGateId]);
	//pGateInList = inNbrs(pGateId);

	/* if a line is a head line, stop */
	if (is_head(pGateId))
	{
		set(changed[pGateId]);
		return(FORWARD);
	}

	/* faulty pGate evaluation */
	if (pGateId == pFault.gateId)
	{
		Log("Current implication gate is faulty gate.");
        Log("Evaluate the faulty gate.");
		iValue = getFaultyGateOutput(GA, pGateId, pFault);

		// printf("getFaultyGateOutput!\n");
		// for (int i = 0; i < GA.n; i++)//lsj
    	// {
        // 	printf("%d ", FANoutput[i]);
    	// }
		// printf("\n");
	
		// Log("Faulty gate: %d iVal: %d output: %d.", pGateId, iValue, FANoutput[pGateId]);
		if (iValue == X)
		{
			if (FANoutput[pGateId] != X)
			{
				for (i = iNumX = 0; i< inCnt(pGateId); i++)
                {
					if (FANoutput[pGateInList[i]]== X)
					{
						iNumX++; j = i;
					}
                }
				if (iNumX == 1)
				{
					/* backward implication */
					// Log("Backward implication.");
					iValue = (FANoutput[pGateId] == D) ? ONE : (FANoutput[pGateId] == DBAR) ? ZERO : FANoutput[pGateId];
					iValue = g_TruthTable1[type[pGateId]][iValue];
					switch (type[pGateId])
					{
					case XOR:
					case XNOR:
						iTempValue = (j == 0) ? FANoutput[pGateInList[1]] : FANoutput[pGateInList[0]];
						if (iTempValue == ONE)
						{
							iValue = g_TruthTable1[NOT][iValue];
						}
						break;
					}
					FANoutput[pGateInList[j]] = iValue;
					set(changed[pGateId]);
					push(g_stack, pGateInList[j]);
					schedule_input(pGateId, j);
					return(BACKWARD);
				}
				else
				{
					Log("Current gate(%d) cannot backward imply, push it into the unjustified stack.", pGateId);
					push(g_unjustStack, pGateId);
					// unjustify.emplace_back(pGateId);
				}
			}
		}
		else if (iValue == FANoutput[pGateId])
		{
			Log("The ival is equal to output. Current gate cannot imply.");
			set(changed[pGateId]);
		}
		else if (FANoutput[pGateId] == X)
		{
			/* forward imp */
			// Log("Forward implication.");
			set(changed[pGateId]);
			FANoutput[pGateId] = iValue;
			push(g_stack, pGateId);
			schedule_output(pGateId);
		}
		else
		{
			Log("Conflicted!");
			return(CONFLICT);
		}
		return(FORWARD);
	}

	// printf("faulty pGate evaluation\n");
	// for (int i = 0; i < GA.n; i++)//lsj
    // {
    //     printf("%d ", FANoutput[i]);
    // }
	// printf("\n");

	/* fault free pGate evaluation */
	// Log("Evaluate the current gate(%d)", pGateId);
	gate_eval1(pGateId, iValue, iTempValue2, i);
	Log("The implicaiton gate: %d iVal: %d output: %d", pGateId, iValue, FANoutput[pGateId]);

	if (iValue == FANoutput[pGateId])
	{
		/* no event */
		Log("Current gate(%d) do not need imply.", pGateId);
		if (iValue != X)
		{
			set(changed[pGateId]);
		}
		return(FORWARD);
	}
	if (FANoutput[pGateId] == X)
	{
		/* forward implication */
		Log("Forward implication.");
		FANoutput[pGateId] = iValue;
		push(g_stack, pGateId);
		set(changed[pGateId]);
		schedule_output(pGateId);
		return(FORWARD);
	}

	if (iValue != X)
	{
		Log("Implication conflicted.");
		return(CONFLICT);
	}		/* conflict */

	/* backward implication */
	Log("Backward implication.");
	switch (type[pGateId])
	{
	case AND:
	case NAND:
	case OR:
	case NOR:
		iTempValue = (type[pGateId] == AND || type[pGateId] == NOR) ? ONE : ZERO;
		if (FANoutput[pGateId] == iTempValue)
		{
			Log("The output value is the non-control value, can imply directly.");
			set(changed[pGateId]);
			for (i = 0; i< inCnt(pGateId); i++)
            {
				if (FANoutput[pGateInList[i]] == X)
				{
					FANoutput[pGateInList[i]] = g_TruthTable1[type[pGateId]][iTempValue];
					push(g_stack, pGateInList[i]);
					schedule_input(pGateId, i);
				}
            }
			i = BACKWARD;
		}
		else
		{
			for (i = iNumX = 0; i< inCnt(pGateId); i++)
			{
				if (FANoutput[pGateInList[i]] == X)
				{
					iNumX++; j = i;
				}
			}
			Log("The number of unknown input values is %d.", iNumX);
			if (iNumX == 1)
			{
				FANoutput[pGateInList[j]] = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
				set(changed[pGateId]);
				push(g_stack, pGateInList[j]);
				schedule_input(pGateId, j);
				i = BACKWARD;
			}
			else
			{
				// unjustified line
				Log("Current gate(%d) cannot backward imply, push it into the unjustified stack.", pGateId);
				push(g_unjustStack, pGateId);
				// unjustify.emplace_back(pGateId);
				i = FORWARD;
			}
		}
		break;
	case BUFF:
	case NOT:
	case PO:
		FANoutput[pGateInList[0]] = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
		set(changed[pGateId]);
		push(g_stack, pGateInList[0]);
		schedule_input(pGateId, 0);
		i = BACKWARD;
		break;
	// case STEM: //lsj added
    //     FANoutput[pGateInList[0]] = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
    //     changed[pGateId] = true;
	// 	push(g_stack, pGateInList[0]);
	// 	schedule_input(pGateId, 0);
    //     break;
	case XOR:
	case XNOR:
		Log("Current implication gate is XOR or XNOR gate.");
		for (i = iNumX = 0; i< inCnt(pGateId); i++)
        {
			if (FANoutput[pGateInList[i]] == X)
			{
				iNumX++; j = i;
			}
        }
		Log("The number of unknown input values of gate(%d) is %d.", pGateId, iNumX);
		if (iNumX == 1)
		{
			iTempValue = (j == 0) ? FANoutput[pGateInList[1]] : FANoutput[pGateInList[0]];
			iValue = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
			if (iTempValue == ONE)
			{
				iValue = g_TruthTable1[NOT][iValue];
			}
			FANoutput[pGateInList[j]] = iValue;
			set(changed[pGateId]);
			push(g_stack, pGateInList[j]);
			schedule_input(pGateId, j);
			i = BACKWARD;
		}
		else
		{
			push(g_unjustStack, pGateId);
			// unjustify.emplace_back(pGateId);
			i = FORWARD;
		}
		break;
	}

	// printf("faulty free pGate evaluation\n");
	// for (int i = 0; i < GA.n; i++)//lsj
    // {
    //     printf("%d ", FANoutput[i]);
    // }
	// printf("\n");   

	return(i);
}


template <class vertex>
bool implyForwardAndBackward(graph<vertex> &GA,int iMaxDPI, bool bBackward, int iLastDPI, FAULT pFault) //imply
{
	int i, iStartDPI;
	int iStatus;
	int pGateId;

	if (bBackward)
	{
		iStartDPI = iLastDPI;
	}
	else
	{
		iStartDPI = 0;
	}
	Log("iStartDPI:%d",iStartDPI);

	while (true)
	{
		/* backward implication */
		Log("backward implication start.");
		if (bBackward)
		{
			for (i = iStartDPI; i >= 0; i--) //clear all stacks from 0 to iStartDPI
			{
				while (!is_empty(g_pEventListStack[i]))
				{
					pGateId = pop(g_pEventListStack[i]);
					Log("current level:%d, current gate:%d",i, pGateId);
					if ((iStatus = eval(GA,pGateId, pFault)) == CONFLICT)
					{
						return(false); //conflict!!
					}
				}
			}
		}
		/* forward implication */
		Log("forward implication start.");
		reset(bBackward);
		Log("iMaxDPI:%d",iMaxDPI);
		for (i = 0; i < iMaxDPI; i++)
		{
			while (!is_empty(g_pEventListStack[i]))
			{
				pGateId = pop(g_pEventListStack[i]);
				Log("current level:%d, current gate:%d",i, pGateId);
				if ((iStatus = eval(GA, pGateId, pFault)) == CONFLICT)
				{
					Log("conflict!");
					return(false); //conflict!!
				}
				else if (iStatus == BACKWARD)
				{
					Log("BACKWARD!");
					iStartDPI = i - 1;
					set(bBackward);
					break;
				}
			}
			if (bBackward)
			{
				break;
			}
		}
		if (!bBackward)
		{
			break;
		}
	}
	return(true);
}

template<class vertex>
void update_Dfrontier(graph<vertex> &GA)
{

	int pGateId;
	//int iFirst; //No use !!

	//iFirst = INFINITY;
	for (int i = 0; i <= g_DfrontierStack.last;)
	{
		pGateId = g_DfrontierStack.list[i];
		switch (FANoutput[pGateId])
		{
		case D:
		case DBAR:
			//lsj added
            // if(outCnt(pGateId)==1 && type[outNbr(pGateId,0)]==STEM){
            //     FANoutput[outNbr(pGateId,0)] = FANoutput[pGateId];
            //     pGateId = outNbr(pGateId,0);
            // }
			for (int j = 0; j< outCnt(pGateId); j++)
			{
				push(g_DfrontierStack, outNbr(pGateId, j));
			}
			delete(g_DfrontierStack, i);
			break;

			
		case X: //Keep X gates
			//if (pGate->index < iFirst)
			//{
			//	iFirst = pGate->index;
			//}
			i++;
			break;

			
		default:
			/* 1 or 0 */
			delete(g_DfrontierStack, i);
		}
	}
}

//检查是否存在从 D-frontier 中的门到主输出的任何路径，使得沿路径的所有线都在X
template<class vertex>
bool Xpath(graph<vertex> &GA,int pGateId)
{
	int i;

	/* base step --- if no X-path exist, return FALSE  */
	if (FANoutput[pGateId] != X || xpath[pGateId] == 0)
	{
		xpath[pGateId] = 0;
		return(false);
	}

	/* base step --- if an X-path exist, return TRUE */
	if ((type[pGateId] == PO) || (xpath[pGateId] == 2))
	{
		xpath[pGateId] = 2;
		return(true);
	}

	/* induction step --- else, go to next step */
	for (i = 0; i < outCnt(pGateId); i++)
    {
		if (Xpath(GA,outNbr(pGateId,i)))
		{
			xpath[pGateId] = 2;
			return(true);
		}
    }

	xpath[pGateId] = 0;
	return(false);
} 


int closest_po(STACKPTR pObjectiveStack, int *piClose)
{
	int i, iDPO;
	int pGate;

	if (ptr_is_empty(pObjectiveStack))
	{
		return 0;
	}
	*piClose = pObjectiveStack->last;
	iDPO = dpo[pObjectiveStack->list[*piClose]];
	for (i = (pObjectiveStack->last) - 1; i >= 0; i--)
		if (dpo[pObjectiveStack->list[i]] < iDPO)
		{
			iDPO = dpo[pObjectiveStack->list[i]];
			*piClose = i;
		}
	pGate = pObjectiveStack->list[*piClose];
	return pGate;
}



template <class vertex>
void setGateTestability(graph<vertex> &GA,int iNoGate) //set_testability
{
	int i, j, iDepth;

	/* cont0 and cont1 */
	//cont0 & cont1 small is better!!
	for (i = 0; i < iNoGate; i++)
	{
		if (is_free(i) || is_head(i))
		{
			cont0[i] = 0;
		}
		else
		{
			iDepth = (-1);
			for (j = 0; j< inCnt(i); j++)
			{
				iDepth = max(iDepth, cont0[inNbr(i,j)]);
			}
			cont0[i] = iDepth + 1; //cont0 = max{pGate->inList[j]->cont0} + 1
		}
		cont1[i] = cont0[i];
	}

	/* iDepth from output */
	//dpo small is better!!
	for (i = iNoGate - 1; i >= 0; i--)
	{
		if (type[i] == PO)
		{
			dpo[i] = 0;
		}
		else
		{
			iDepth = (-1);
			for (j = 0; j<outCnt(i); j++)
			{
				iDepth = max(iDepth, (int)dpo[outNbr(i,j)]);
			}
			dpo[i] = iDepth + 1;  //dpo = max{pGate->outList[j]->dpo} + 1
		}
	}
}

void printStackInfo()
{
    Log("init obj: ");
    for (int i = 0; i <= g_initObjStack.last; i++)
    {
    Log(" (%d, %d, %d) ", g_initObjStack.list[i], numzero[g_initObjStack.list[i]], numone[g_initObjStack.list[i]]);
    }
    Log("\n");
    Log("cur obj: ");
    for (int i = 0; i <= g_curObjStack.last; i++)
    {
        Log(" (%d, %d, %d) ", g_curObjStack.list[i], numzero[g_curObjStack.list[i]], numone[g_curObjStack.list[i]]);
    }
    Log("\n");
    Log("head obj: ");
    for (int i = 0; i <= g_headObjStack.last; i++)
    {
        Log(" (%d, %d, %d) ", g_headObjStack.list[i], numzero[g_headObjStack.list[i]], numone[g_headObjStack.list[i]]);
    }
    Log("\n");
    Log("fan obj: ");
    for (int i = 0; i <= g_fanObjStack.last; i++)
    {
    Log(" (%d, %d, %d) ", g_fanObjStack.list[i], numzero[g_fanObjStack.list[i]], numone[g_fanObjStack.list[i]]);
    }
    Log("\n");
    Log("final obj: ");
    for (int i = 0; i <= g_finalObjStack.last; i++)
    {
        Log(" (%d, %d, %d) ", g_finalObjStack.list[i], numzero[g_finalObjStack.list[i]], numone[g_finalObjStack.list[i]]);
    }
    Log("\n");
}

template <class vertex>
int backtrace(graph<vertex> &GA ,int iState)
{
	int i;
    // int k;
	int j;
	int v1;
	int a_curr_obj;
	int n0, n1, nn0, nn1;
	int easiest=0;
    int easy_cont;
    Log("istate = 81");
	/* box 1: Initialization of objective and its logic level */
	if (iState == 81)
	{
		// printStackInfo();
        Log("Multiple backtrace from the set of initial objectives.");
        copy(g_initObjStack, g_curObjStack);
        Log("Let the set of initial objectives be the set of current objectives.");
		for (i = 0; i <= g_initObjStack.last; i++)
		{
			a_curr_obj = g_initObjStack.list[i];
			switch (FANoutput[a_curr_obj])
			{
			case ZERO:
			case DBAR:
				setline(a_curr_obj, 1, 0); break;	/* unjustified lines */
			case ONE:
			case D:
				setline(a_curr_obj, 0, 1); break;
			default:
				/* Dfrontier */
                Log("Current objective(%d) is a D-frontier.", a_curr_obj);
				switch (type[a_curr_obj])
				{
				case AND:
				case NOR:
					setline(a_curr_obj, 0, 1); break;
				case NAND:
				case OR:
					setline(a_curr_obj, 1, 0); break;
				case XOR:
				case XNOR:
					setline(a_curr_obj, 1, 0); break;
				}
			}

		}
		iState = 82;
	}

	while (true)
	{
        switch (iState)
		{
			/* Box 2,3,4 of figure 8 */
		case 82:
            Log("istate = 82");
            // printStackInfo();
			if (is_empty(g_curObjStack)) 		/* box 2 */
			{
				Log("The set of current objectives is empty.");
                if (is_empty(g_fanObjStack))
				{
					Log("The set of fanout-point objectives is empty.");
                    iState = 103;
                    Log("istate = 103");
				}	/* box 4 */
				else
				{
					Log("The set of fanout-point objectives is not empty.");
                    iState = 86;
                    Log("istate = 86");
                    // printStackInfo();
				}
			}
			else
			{
				Log("The set of  current objectives is not empty.");
                a_curr_obj = pop(g_curObjStack);	/* box 3 */
				iState = 85;
                Log("istate = 85");
                // printStackInfo();
			}
			break;

			/* Box 5,9,10,11,12 of figure 8 */
		case 85:
            // Log("istate = 85");
            Log("whether the current objective(%d) is the head line", a_curr_obj);
			if (is_head(a_curr_obj))
			{
				/* box 5 */
                Log("The objective line(%d) is head line.", a_curr_obj);
				push(g_headObjStack, a_curr_obj);	/* box 12 */
 
			}
			else
			{
				/* box 9,10,11 */
                Log("The objective line(%d) is not head line. type is %d", a_curr_obj,type[a_curr_obj]);
				switch (type[a_curr_obj])
				{
				case AND:
					n0 = numzero[a_curr_obj];
					n1 = numone[a_curr_obj];
					v1 = ZERO;
					break;
				case OR:
		            n0 = numzero[a_curr_obj];
					n1 = numone[a_curr_obj];
					v1 = ONE;
					break;
				case NAND:
		            n1 = numzero[a_curr_obj];
					n0 = numone[a_curr_obj];
					v1 = ZERO;
					break;
				case NOR:
		            n1 = numzero[a_curr_obj];
					n0 = numone[a_curr_obj];
					v1 = ONE;
					break;
				case NOT:
		            n1 = numzero[a_curr_obj];
					n0 = numone[a_curr_obj];
					v1 = X;
					break;
				case XOR:
					j = 0;
					if ((v1 = FANoutput[inNbr(a_curr_obj, j)]) == X)
					{
						v1 =  FANoutput[inNbr(a_curr_obj, ++j)];
					}
					if (v1 == ONE)
					{
                        n1 = numzero[a_curr_obj];
                        n0 = numone[a_curr_obj];
					}
					else
					{
                        n0 = numzero[a_curr_obj];
                        n1 = numone[a_curr_obj];
					}
					v1 = X;
					break;
				case XNOR:
					j = 0;
					if ((v1 = FANoutput[inNbr(a_curr_obj, j)]) == X)
					{
						v1 =  FANoutput[inNbr(a_curr_obj, ++j)];
					}
					if (v1 == ZERO)
					{
                        n1 = numzero[a_curr_obj];
                        n0 = numone[a_curr_obj];
					}
					else
					{
                        n0 = numzero[a_curr_obj];
                        n1 = numone[a_curr_obj];
					}
					v1 = X;
					break;
				default:
					/* BUFF, PO, PI */
                    n0 = numzero[a_curr_obj];
                    n1 = numone[a_curr_obj];
					v1 = X;
					break;
				}
                // Log("v1:%d",v1);
				/* Find the easiest input. */
				auto *input = inNbrs(a_curr_obj);
				easy_cont = INT_MAX;
				easiest = 0;
				if (v1 == ZERO)  
				{
					/* and, nand */
					for (i = 0; i< inCnt(a_curr_obj); i++)
                    {
						if (FANoutput[input[i]] == X)
						{
							if (easy_cont > cont0[input[i]])
							{
								easy_cont = cont0[input[i]];
								easiest = i;
							}
						}
                    }
				}
				else
				{
					/* or, nor, xor,xnor */
					for (i = 0; i< inCnt(a_curr_obj); i++)
                    {
						if (FANoutput[input[i]] == X)
						{
							if (easy_cont > cont1[input[i]])
							{
								easy_cont = cont1[input[i]];
								easiest = i;
							}
						}
                    }
				}
               

				for (i = 0; i<inCnt(a_curr_obj); i++)
				{
				                  
                    if (FANoutput[input[i]] == X)
					{

                        if (i == easiest)
						{
							nn0 = n0; nn1 = n1;                           
						}
						else if (v1 == ZERO)
						{
                            nn0 = 0; nn1 = n1;
						}
						else if (v1 == ONE)
						{
                            nn0 = n0; nn1 = 0;
						}
						else
						{
							/* xor,xnor */
                            if (n0 > n1)
							{
								nn0 = n0; nn1 = n1;
							}
							else
							{
								nn0 = n1; nn1 = n0;
							}
						}
						if (nn0 > 0 || nn1 > 0)
						{
							if (outCnt(input[i]) > 1 || type[input[i]] == STEM)
							{
								if (numzero[input[i]] == 0 && numone[input[i]] == 0)
								{
									push(g_fanObjStack, input[i]);
									Log("The input gate is fanout point.");
								}
								numzero[input[i]] += nn0;
								numone[input[i]] += nn1;
							}                                   
							else
							{
								setline(input[i], nn0, nn1);
								push(g_curObjStack, input[i]);
							}

						}
					}
				} /* for */
			}
            // printStackInfo();
			iState = 82;
			break;

			/* Box 6,7,8 of figure 8 */
		case 86:
            // printStackInfo();
			a_curr_obj = closest_po(&g_fanObjStack, &i);
            	/* box 6 */
            Log("Take out a fanout-point objective p(%d) closest to primary output.", a_curr_obj);
			delete(g_fanObjStack, i);
			if (FANoutput[a_curr_obj] != X)
			{
				iState = 82; break;
			}
			if (is_reachable_from_fault(a_curr_obj))
			{
				/* box 7 */
                Log("is_reachable_from_fault");
                Log("The fanout-point objective p(%d) is reachable from the fault line.", a_curr_obj);
                // printStackInfo();
				iState = 85; break;
			}
			if (!is_conflict(a_curr_obj))
			{
				/* box 8 */
                Log("is_conflict");
                Log("n0(p) and n1(p) are not both non-zero.");
                // printStackInfo();
				iState = 85; break;

			}
			Log("Both n0(p) and n1(p) are both non-zero.");
            // Let the fanout point objective be final objective to assign a value
            Log("Let the fanout-point objective be final objective to assign a value.");
            push(g_finalObjStack, a_curr_obj);		/* box 12 in figure 10 */
			iState = 93;
            Log("istate = 93");
            // printStackInfo();
			break;

		default:
			return(iState);
		}
	}
}


template <class vertex>
void find_final_objective(graph<vertex> &GA,bool *backtrace_flag, bool fault_propagated_to_po, int nog, int *last_Dfrontier)
{
	int i;
	int p;
	int state;


	if (*backtrace_flag)
	{   
        Log("The backtrace flag is on. Go to reset.");
		state = 107;
	}	/* box 1 */
	else if (is_empty(g_fanObjStack))
	{
		Log("The set of fanout-point objectives is empty.");
        state = 103;
	}	/* box 2 */
	else
	{
		Log("Multiple backtrace from a fanout-point objective.");
        state = 86;
	}
	while (true)
		switch (state)
		{
		case 103:
            Log("istate = 103");
            // printStackInfo();
			/* box 3,4,5,6 */
			if (is_empty(g_headObjStack))
			{
				Log("The set of head objectives is empty. Go to reset.");
                state = 107;
			}	/* box 3 */
			else
			{
				Log("Take out a head objective.");
                p = g_headObjStack.list[0];
				for (i = 1; i <= g_headObjStack.last; i++)
                {
					g_headObjStack.list[i - 1] = g_headObjStack.list[i];
                }
				delete_last(g_headObjStack);
				if (FANoutput[p] == X)
				{
					/* box 4,5 */
                    Log("The head line(%d) is unspecified.", p);
					push(g_finalObjStack, p);		/* box 6 */
                    Log("Let the head objective be final objective.");
                    // printStackInfo();
					state = 93;
				}
				else
				{   
                    // printStackInfo();
                    Log("The head line (%d) is specified.", p);
					state = 103;
				}
			}
			break;
		case 107:
            Log("istate = 107");
            Log("Reset backtrace flag and let all the sets of objectives be empty.");
            // printStackInfo();
			/* box 7,8,9,10,11 */
            reset(*backtrace_flag);	/* box 7 */
			for (i = 0; i < nog; i++)
			{
				/* initialization */
				numzero[i] = 0;
				numone[i] = 0;
			}
			stackclear(g_initObjStack);
			stackclear(g_curObjStack);
			stackclear(g_fanObjStack);
			stackclear(g_headObjStack);
			stackclear(g_finalObjStack);
            Log("Is there any unjustified lines?");
			if (!is_empty(g_unjustStack))
			{
				/* box 8 */
                Log("There is some unjustified lines.");
				copy(g_unjustStack, g_initObjStack);	/* box 9 */
                Log("Let all the unjustified lines be the set of initial objectives.");

				if (fault_propagated_to_po)
				{
					/* box 10 */
                    Log("Fault signal propagated to a primary output.");
					(*last_Dfrontier) = 0; 
					state = 81;
					break;
				}
			}
            Log("Dfrontier:");
			for(int i = g_DfrontierStack.last; i >= 0; i--){
				Log("%d",g_DfrontierStack.list[i]);
			}
			
            Log("Fault signal does not propagate to a primary output.");
			(*last_Dfrontier) = closest_po(&g_DfrontierStack, &i); 
            // numone[*last_Dfrontier] = 0;
            // numzero[*last_Dfrontier] = 0;
			push(g_initObjStack, (*last_Dfrontier));

            Log("Add a gate in D-frontier(%d) to the set of initial objectives.", *last_Dfrontier);

			state = 81;
			break;
		case 86:
		case 81:
            Log("istate = 81");
            Log("Multiple backtrace.");
			state = backtrace(GA,state);
			break;
		default:
            Log("return");
			return;				/* exit */
		}
}


template <class vertex>
bool backtrack(graph<vertex> &GA ,int faulty_gateId, int *last, int nog)
{
	int p;
	int i, j;
	int value;

	while (!is_empty(g_tree))
		if (is_flagged(current_node))
		{
			delete_last(g_tree);
		}
		else
		{
			/* update & remove duplicate unjustified lines */
			for (i = g_unjustStack.last; i >= 0; i--)
			{
				p = g_unjustStack.list[i];
				if (is_justified(p))
				{
					delete(g_unjustStack, i);
				}
				else
				{
					set(changed[p]);
					push(g_finalObjStack, p);
				}
			}
			while (!is_empty(g_finalObjStack))
				reset(changed[pop(g_finalObjStack)]);

			/* restore and schedule events */
			value = A_NOT(FANoutput[current_node.gate_id]);
			set(current_node.flag);
			for (i = current_node.pstack; i <= g_stack.last; i++)
			{
				p = g_stack.list[i];
				FANoutput[p] = X;
				reset(changed[p]);
				for (j = 0; j< outCnt(p); j++)
                {
					reset(changed[outNbr(p,j)]);
                }
			}
			*last = 0;
			for (i = current_node.pstack + 1; i <= g_stack.last; i++)
			{
				p = g_stack.list[i];
				for (j = 0; j< outCnt(p); j++)
				{
					if (FANoutput[outNbr(p,j)] != X)
					{
						if (!is_justified(outNbr(p,j)))
						{
							push(g_unjustStack, outNbr(p,j));
						}
						pushevent(outNbr(p,j));
					}
					if ( (*last) < level[outNbr(p,j)])
					{
						(*last) = level[outNbr(p,j)];
					}
				}
			}
			g_stack.last = current_node.pstack;
			FANoutput[current_node.gate_id] = value;
			if (is_head(current_node.gate_id))
			{
				set(changed[current_node.gate_id]);
			}
			else
			{
				pushevent(current_node.gate_id);
			}
			schedule_output(current_node.gate_id);

			/* update Dfrontier */
			stackclear(g_DfrontierStack);
			push(g_DfrontierStack, faulty_gateId);
			update_Dfrontier(GA);

			/* update unjustified set */
			for (i = g_unjustStack.last; i >= 0; i--)
            {
				if (FANoutput[g_unjustStack.list[i]] == X)
				{
					delete(g_unjustStack, i);
				}
            }
			/* reset xpath */
			for (i = faulty_gateId; i < nog; i++)
            {
				xpath[i] = 1;
            }
			return(true);
		}
	return(false);
}

template <class vertex>
void restore_faults(graph<vertex> &GA,FAULT pFault)
{
	int i, j;
	int pGateId;
	int k;
	int value, gtype;

	stackclear(g_fanObjStack);
	pGateId = pFault.gateId;
	if (pFault.faultLine != OUTFAULT)
	{
		pGateId =  inNbr(pGateId, pFault.faultLine); //inList[pFault->line] -> pFault -> pGate
	}
	push(g_fanObjStack, pGateId);

	while (true)
	{
		pGateId = outNbr(pGateId,0);
		for (i = 0; i< inCnt(pGateId); i++)
        {
			if (FANoutput[inNbr(pGateId,i)] == ZERO || FANoutput[inNbr(pGateId,i)] == ONE)
			{
				push(g_fanObjStack, inNbr(pGateId,i));
			}
        }
		if (is_head(pGateId))
		{
			break;
		}
	}

	while (!is_empty(g_fanObjStack))
	{
		pGateId = pop(g_fanObjStack);    
		if (FANoutput[pGateId] == D)
		{
			FANoutput[pGateId] = ONE;
		}
		else if (FANoutput[pGateId]  == DBAR)
		{
			FANoutput[pGateId]  = ZERO;
		}
		if (!(type[pGateId] == PI || FANoutput[pGateId]  == X))
		{
			stackclear(g_curObjStack);
			push(g_curObjStack, pGateId);
			while (!is_empty(g_curObjStack))
			{
				pGateId = pop(g_curObjStack);
				switch (type[pGateId])
				{
				case PI:
					break;
				case XOR:
					FANoutput[inNbr(pGateId,0)] = ZERO;
					push(g_curObjStack, inNbr(pGateId,0));
					for (j = 1; j< inCnt(pGateId) ; j++)
					{
						FANoutput[inNbr(pGateId,j)] = FANoutput[pGateId];
						push(g_curObjStack, inNbr(pGateId,j));
					}
					break;
				case XNOR:
					FANoutput[inNbr(pGateId,0)]= ONE;
					push(g_curObjStack, inNbr(pGateId,0));
					for (j = 1; j< inCnt(pGateId); j++)
					{
						FANoutput[inNbr(pGateId,j)] = FANoutput[pGateId];
						push(g_curObjStack, inNbr(pGateId,j));
					}
					break;
				case PO:
				case BUFF:
				case NOT:
					FANoutput[inNbr(pGateId,0)] = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
					push(g_curObjStack, inNbr(pGateId,0));
					break;
				default:
					/* and,or,nor,nand */
					value = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
					gtype = (type[pGateId] == AND || type[pGateId] == NAND) ? ONE : ZERO;
					if (value == gtype)
					{
						for (j = 0; j< inCnt(pGateId); j++)
						{
							FANoutput[inNbr(pGateId,j)] = value;
							push(g_curObjStack, inNbr(pGateId,j));
						}
					}
					else
					{
						k = 0;
						for (j = 1; j< inCnt(pGateId); j++)
                        {
							if (level[inNbr(pGateId,j)] < level[inNbr(pGateId,k)])
							{
								k = j;
							}
                        }
						FANoutput[inNbr(pGateId,k)] = value;
						push(g_curObjStack, inNbr(pGateId,k));
					} 
					break;
				} /* switch */
			} /* while */
		} /* if */
	} /* while */
}

template <class vertex>
void justify_free_lines(graph<vertex> &GA, int iNoPI, FAULT pFault, FAULT pCurFault)
{
	int i, j, k;
	int pGateId;
	int iValue, gtype;

	for (i = 0; i < iNoPI; i++)
	{
		// if (g_iHeadGateIndex[i] < 0) //never happen!!
		// {
		// 	break;
		// }
		pGateId = g_iHeadGateIndex[i];
		if (pGateId == pCurFault.gateId && pFault.gateId != -1)
		{
			restore_faults(GA, pFault);
			continue;
		}
		if (FANoutput[pGateId] == D)
		{
			FANoutput[pGateId] = ONE;
		}
		else if (FANoutput[pGateId] == DBAR)
		{
			FANoutput[pGateId] = ZERO;
		}
		if (!(type[pGateId] == PI || FANoutput[pGateId] == X)) //pGate->type != PI && pGate->output != X
		{
			//	  	headlines1[i]=p->output;
			stackclear(g_curObjStack);
			push(g_curObjStack, pGateId);
			while (!is_empty(g_curObjStack))
			{
				pGateId = pop(g_curObjStack);
				switch (type[pGateId])
				{
				case PI:
					break;
				case XOR:
                    FANoutput[inNbr(pGateId,0)] = ZERO;
					push(g_curObjStack, inNbr(pGateId,0));
					for (j = 1; j< inCnt(pGateId); j++)
					{
						FANoutput[inNbr(pGateId,j)] = FANoutput[pGateId];
						push(g_curObjStack, inNbr(pGateId,j));
					}
					break;
				case XNOR:
					FANoutput[inNbr(pGateId,0)] = ONE;
					push(g_curObjStack, inNbr(pGateId,0));
					for (j = 1; j< inCnt(pGateId); j++)
					{
						FANoutput[inNbr(pGateId,j)] = FANoutput[pGateId];
						push(g_curObjStack, inNbr(pGateId,j));
                    }
					break;
				case PO:
				case BUFF:
				case NOT:
					FANoutput[inNbr(pGateId,0)] = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
					push(g_curObjStack, inNbr(pGateId,0));
					break;
				default:
					/* and,or,nor,nand */
					iValue = g_TruthTable1[type[pGateId]][FANoutput[pGateId]];
					gtype = (type[pGateId] == AND || type[pGateId] == NAND) ? ONE : ZERO;
					if (iValue == gtype)
					{
						for (j = 0; j< inCnt(pGateId); j++)
						{
							FANoutput[inNbr(pGateId,j)] = iValue;
							push(g_curObjStack,inNbr(pGateId,j));
						}
					}
					else
					{
						k = 0;
						for (j = 1; j< inCnt(pGateId); j++){
							if (level[inNbr(pGateId,j)]<level[inNbr(pGateId,k)])
							{
								k = j;
							}
                        }
						FANoutput[inNbr(pGateId,k)] = iValue;
						push(g_curObjStack, inNbr(pGateId,k));
					} 
					break;
				} /* switch */
			} /* while */
		} /* if */
	} /* for */
}


template <class vertex>
int fan1(graph<vertex>& GA, int iNoGate, int iMaxLevelAdd2, int iNoPI, int iNoPO, FAULT pCurrentFault, int iMaxBackTrack, int *piNoBackTrack)
{
	//INPUT:  pCurrentFault
	//OUTPUT:  iState & g_iPatternsForOneTime & g_net[i] (PI)
    *piNoBackTrack = 0;
	int i;
	int pGateId, pTempGateId;
	int iLastDPI;
	int pLastDFtrGateId;
	bool bBackwardFlag, bBacktraceFlag, bFaultPropagatedToPO;
	bool bDfrontierChanged, bDone;
	int iState;
	FAULT pOriginFault = {-1,-1,0,0};
    // memset(changed, (bool)0, sizeof(bool)*n);
    // memset(FANoutput,X,sizeof(char)*n);
    // memset(u_path, NULL, sizeof(LINKTYPE*)*n);
    // memset(freach, (bool)0, sizeof(bool)*n);
    memset(numzero, 0, sizeof(int)*n);
    memset(numone, 0, sizeof(int)*n);
	reset(bDone);
	reset(bBackwardFlag);
	set(bBacktraceFlag);
	reset(bFaultPropagatedToPO);
	pLastDFtrGateId = 0;

	pGateId = pCurrentFault.gateId;

    Log("init");
	initNetAndFreach(GA,iNoGate, pGateId, iMaxLevelAdd2);
    Log("init end");

 	/* initializaiton */
	if (pCurrentFault.faultLine != OUTFAULT) //pGate -----> pCurrentFault -----> pGate
	{
		//					  pCurrentFault -------> pGate
		pGateId = inNbr(pGateId,pCurrentFault.faultLine);//pGate to upper side
        //pGate -------> pCurrentFault
	}
	if (is_free(pGateId)) //come here only once for s9234tc.bench,       pGate->outCount == 1
	{
		/* box 1 */
        Log("The faulty line is free line.");
		pOriginFault.gateId = pCurrentFault.gateId;
		pOriginFault.faultLine= pCurrentFault.faultLine;
		pOriginFault.ftype = pCurrentFault.ftype;
        Log("Propagate the fault to a head line.");
        // Propagate the fault to a head line
        Log("PropFault2Headline start!");
		propFault2Headline(GA,pCurrentFault);
        
        Log("PropFault2Headline end!");
		iLastDPI = 0;
	}
	else
	{
		iLastDPI = refFaultyGateOutput(GA,pCurrentFault);
	}
	Log("iLastDPI:%d",iLastDPI);
	if (iLastDPI == (-1)) //Never come here !!
	{
		//STOP***********************STOP
		return(NO_TEST);
	}
    
	pGateId = pCurrentFault.gateId; //pGate comes back to pCurrentFault->gate
	Log("pGateId:%d",pGateId);
    // printf("pGateId:%d,faultline:%d,faulttype:%d\n",pGateId,pCurrentFault.faultLine, pCurrentFault.ftype);

	push(g_DfrontierStack, pGateId);
    Log("Unique sensitize start!");
	i = unique_sensitize(GA,pGateId, pGateId);//唯一敏化
	if ((iLastDPI = max(i, iLastDPI)) > 0)
	{
		set(bBackwardFlag);
	}
	iState = 93;
    Log("Unique sensitize end!");
    Log("iLastDPI:%d",iLastDPI);


// 	/* main loop of fan algorithm */
	while (bDone == false)
	{
		switch (iState)
		{
		case 93:
			/* box 3,4,5,6 */
            Log("Implication start!");
			if (!implyForwardAndBackward(GA, iMaxLevelAdd2, bBackwardFlag, iLastDPI, pCurrentFault))
			{
				/* box 3 */
                Logerr("Implication is conflicted! Backtrack.");
				iState = 98;
				break;
			}
            Log("Implication end!");
			if (FANoutput[pGateId] == ZERO || FANoutput[pGateId] == ONE) //???
			{
				iState = 98; break;
			}
            // After implication, update the unjustified lines.
            Log("Update the unjustified lines.");
// 			/* update unjustified lines and delete duplicated lines 
// 																	final_obj should be empty */
			for (i = g_unjustStack.last; i >= 0; i--)
			{
				pTempGateId = g_unjustStack.list[i];
				if (is_justified(pTempGateId))
				{
					delete(g_unjustStack, i);
				}
				else
				{
					set(changed[pTempGateId]);
					push(g_finalObjStack, pTempGateId);
				}
			}
            Log("Update the unjustified lines finished.");
            // Box 4: is continuation of backtrace meaningful.
            Log("Make sure whether the continuation of backtrace meaningful.");
			while (!is_empty(g_finalObjStack))
				reset(changed[pop(g_finalObjStack)]);

			/* check for backtrace */
			for (i = g_initObjStack.last; i >= 0; i--)
            {
				if (is_justified(g_initObjStack.list[i]))
				{
					delete(g_initObjStack, i);
				}
            }
            Log("Fault signal propagated to a primary output?");
			reset(bFaultPropagatedToPO);
			for (i = 0; i < iNoPO; i++)
            {
				if (FANoutput[index_PO[i]] == D || FANoutput[index_PO[i]] == DBAR)
				{
					Log("Yes, fault signal propagated to a primary output.");
                    set(bFaultPropagatedToPO);
					break;
				}
            }
			if (pLastDFtrGateId != 0) 
			{
				Log("In the last process of backtrace, the initial objectives have a D-frontier.");
                if (FANoutput[pLastDFtrGateId] == X)
				{
					Log("The D-frontier has not changed.");
                    reset(bDfrontierChanged);
				}
				else
				{
					Log("The D-frontier has changed.");
                    set(bDfrontierChanged);
				}
			}
			else
			{
				set(bDfrontierChanged);
			}

			if (is_empty(g_initObjStack) && bDfrontierChanged)	/* box 4, 4-1 */
			{
				Log("It is not meaningful to continue the last backtrace.");
                Log("Set backtrace flag.");
                set(bBacktraceFlag);
			}

			if (bFaultPropagatedToPO)
			{
				iState = 99;				/* box 4-3 */
				for (i = g_unjustStack.last; i >= 0; i--)
                {
					if (is_unjustified(g_unjustStack.list[i]) && is_bound(g_unjustStack.list[i]))
					{
						Log("There is some unjustified bound lines. Need backtrace.");
                        iState = 97;
						break;
					}
                }
			}
			else
			{
				/* box 5 */
				/* update Dfrontier */
                Log("Fault signal does not propagate the primary output.");
				Log("Dfrontier count is %d",g_DfrontierStack.last + 1);
				if (!is_empty(g_DfrontierStack))
				{
					Log("Update D-frontier start!");
                    update_Dfrontier(GA);
                    Log("Update D-frontier end!");
				}
				Log("Dfrontier count is %d",g_DfrontierStack.last + 1);

				for (i = pGateId; i < iNoGate; i++)  //???
                {
					if (xpath[i] == 2)
					{
						xpath[i] = 1;
					}
                }

				for (i = g_DfrontierStack.last; i >= 0; i--)
                {
					if (!Xpath(GA,g_DfrontierStack.list[i]))
					{
						delete(g_DfrontierStack, i);
					}
                }

				if (is_empty(g_DfrontierStack))
				{
					Log("The number of D-frontier is zero. Backtrack.");
                    iState = 98;
				}
				else if (g_DfrontierStack.last == 0)
				{
					/* box 6 */
					Log("The number of D-frontier is one. Unique sensitize.");
                    Log("Unique sensitize start!");
                    Log("Dfrontier: %d", g_DfrontierStack.list[0]);
                    iLastDPI = unique_sensitize(GA, g_DfrontierStack.list[0], pGateId);
                    Log("iLastDPI:%d",iLastDPI);                    
                    if (iLastDPI  > 0)
					{
						Log("Unique sensitize end!");
                        set(bBackwardFlag);
						iState = 93;
					}
					else if (iLastDPI == 0)  
					{
						iState = 93;
					}
					else
					{
						iState = 97; 
					}

				}
				else
				{   
                    Log("The number of D-frontier is many. Backtrace.");
					iState = 97;
				}
			}
			break;

		case 97:
			/* box 7 */
            Log("Find final objective start!");
			Log("bBacktraceFlag :%d",bBacktraceFlag);
			find_final_objective(GA, &bBacktraceFlag, bFaultPropagatedToPO, iNoGate, &pLastDFtrGateId);
            Log("Find final objective end!");
			while (!is_empty(g_finalObjStack))
			{
				pTempGateId = pop(g_finalObjStack);
				if (numzero[pTempGateId] > numone[pTempGateId])
				{
					FANoutput[pTempGateId] = ZERO;
				}
				else
				{
					FANoutput[pTempGateId] = ONE;
				}
				(g_tree.last)++;
				current_node.gate_id = pTempGateId;
				reset(current_node.flag);
				push(g_stack, pTempGateId);
				if (is_head(pTempGateId))
				{
					set(changed[pTempGateId]);
				}
				else
				{
					pushevent(pTempGateId);
				}		/**** study ****/
				schedule_output(pTempGateId);
				current_node.pstack = g_stack.last;
                Log("find a final obj");
			}

			reset(bBackwardFlag);
			iState = 93;
			break;

		case 98:
			/* box 8 */
            Log("backtrack start");
			if (is_flagged(current_node))
			{
				(*piNoBackTrack)++;
			}
			Log("piNoBackTrack : %d",(*piNoBackTrack));
			for (i = 0; i < iMaxLevelAdd2; i++)
            {
				while (!is_empty(g_pEventListStack[i]))
					reset(changed[pop(g_pEventListStack[i])]);
            }
			// printf("num of backtrack: %d\n",(*piNoBackTrack));
			if (*piNoBackTrack > iMaxBackTrack)
			{   
                Log("The number of back track is up to the maximum number.");
				// printf("num of backtrack: %d\n",(*piNoBackTrack));
				iState = OVER_BACKTRACK;
				//printf("over backtrack! current fault id :%d, faultline:%d, fault type: sa%d\n", pCurrentFault.gateId, pCurrentFault.faultLine, pCurrentFault.ftype);
			}
			else if (backtrack(GA, pGateId, &iLastDPI, iNoGate))
			{
				Log("Backtrack end!");
                Log("Set the backtrace flag.");
                if (iLastDPI > 0)
				{
					set(bBackwardFlag);
				}
				else
				{
					reset(bBackwardFlag);
				}
				set(bBacktraceFlag);
				pLastDFtrGateId = 0;
                Log("Go to implication.");
				iState = 93;
			}
			else
			{
				iState = NO_TEST;
				// printf("NO_TEST\n");
				//printf("no test! current fault id :%d, faultline:%d, fault type: sa%d\n", pCurrentFault.gateId, pCurrentFault.faultLine, pCurrentFault.ftype);
			}
			break;

		case 99:
            /* box 9 */
            Log("Justify free lines start!");
			justify_free_lines(GA, iNoPI, pOriginFault, pCurrentFault);
            Log("Justify free lines end!");
			iState = TEST_FOUND;
			
			break;
		default:
			set(bDone);
		}
    }

	if (pOriginFault.gateId != -1)
	{
		pCurrentFault.gateId = pOriginFault.gateId;
		pCurrentFault.faultLine = pOriginFault.faultLine;
		pCurrentFault.ftype = pOriginFault.ftype;
		// free((FAULT)pOriginFault);
	}

	// for(int i = 0; i < iNoGate; i++){ //lsj added
	// 	if(type[i] == STEM && type[inNbr(i,0)] == PI){
	// 		FANoutput[inNbr(i,0)] = FANoutput[i];
	// 	}
	// }
    // free(g_unjustStack.list);
    // free(g_initObjStack.list);
	// free(g_curObjStack.list);
	// free(g_fanObjStack.list);
	// free(g_headObjStack.list);
	// free(g_finalObjStack.list);
	// free(g_DfrontierStack.list);

    // free(g_stack.list);
	// free(g_tree.list);
	return(iState);
    Log("end!");
}


template<class vertex>
int dynamic_unique_sensitize(graph<vertex> &GA, int *Dfront, int nod, int maxdpi, vector<int> dom_array, int faulty_gate)
{
	int ndom = 0, ngate;
	int dy_id2;
	int i, j, k; //,l,m,n;
	int gut, nextgate; //gate, g not used
	//int Dominator;
	int ndominator = 0,no_dom = 0; //new_dom not used
	int flag = false;
	int debug = true;
	int v1,send = -1;
	vector<int> xpo; 
	int nxpo;

	/* pass 1: D-frontier propagation */
	++dy_id;
	for (i = nxpo = 0; i <= nod; i++)
	{
		gut = Dfront[i];
		if (freach1[gut] < dy_id)
		{
			push(g_pEventListStack[level[gut]], gut);
			freach1[gut] = dy_id;
		}
	}
	for (i = 0; i < maxdpi; i++)
	{
		while (!is_empty(g_pEventListStack[i]))
		{
			gut = pop(g_pEventListStack[i]);
			if (type[gut] == PO)
			{
				// xpo[nxpo++] = gut;
				xpo.emplace_back(gut);
				nxpo++;
			}
			for (j = 0; j< outCnt(gut); j++)
			{
				nextgate = outNbr(gut, j);
				if ((FANoutput[nextgate] == X) && (freach1[nextgate] < dy_id))
				{
					push(g_pEventListStack[level[nextgate]], nextgate);
					freach1[nextgate] = dy_id;
				}
			}
		}
	}

	/* pass 2: Backward propagation --- X-path */
	dy_id2 = dy_id + 1;
	for (i = 0; i < nxpo; i++)
	{
		gut = xpo[i];
		freach1[gut] = dy_id2;
		for (j = 0; j< inCnt(gut); j++)
		{
			nextgate = inNbr(gut, j);
			if (freach1[nextgate] == dy_id)
			{
				push(g_pEventListStack[level[nextgate]], nextgate);
				freach1[nextgate] = dy_id2;
			}
		}
	}
	for (i = maxdpi - 1; i >= 0; i--)
		while (!is_empty(g_pEventListStack[i]))
		{
			gut = pop(g_pEventListStack[i]);
			for (j = 0; j< inCnt(gut); j++)
			{
				nextgate = inNbr(gut, j);
				if (freach1[gut] == dy_id)
				{
					push(g_pEventListStack[level[nextgate]], nextgate);
					freach1[gut] = dy_id2;
				}
			}
		}

	/* pass 3: Compute dominators */
	dy_id = dy_id2 + 1;
	for (i = ngate = 0; i <= nod; i++)
	{
		gut = Dfront[i];
		push(g_pEventListStack[level[gut]], gut);
		ngate++;
		freach1[gut] = dy_id;
	}
	for (i = k = 0; i < maxdpi; i++)
	{
		if (ngate == 1 && g_pEventListStack[i].last == 0)
		{
			// dom_array[k++] = g_pEventListStack[i].list[0];
			dom_array.emplace_back(g_pEventListStack[i].list[0]);
			k++;
		}
		while (!is_empty(g_pEventListStack[i]))
		{
			gut = pop(g_pEventListStack[i]);
			ngate--;
			if (type[gut] == PO)
			{
				ngate = INT_MAX; break;
			}
			for (j = 0; j< outCnt(gut); j++)
			{
				nextgate = outNbr(gut, j);
				if (freach1[nextgate] == dy_id2)
				{
					push(g_pEventListStack[level[nextgate]], nextgate);
					freach1[nextgate] = dy_id;
					ngate++;
				}
			}
		}
		if (ngate == INT_MAX)
		{
			break;
		}
	}

	/* Assign non-controlling values to dominators */
	send = (-1);
	while (--k >= 0)
	{
		gut = dom_array[k];
		/*
			printf("dominator: gut=%d #Dfrontier=%d",gut->index,nod+1);
		*/
		if (gut == faulty_gate)
		{
			continue;
		}
		v1 = (type[gut] == AND || type[gut] == NAND) ? ONE : (type[gut] == OR || type[gut] == NOR) ? ZERO : X;
		if (v1 != X)
		{
			for (i = 0; i<inCnt(gut); i++)
			{
				nextgate = inNbr(gut, i);
				if (freach1[nextgate] < dy_id && FANoutput[nextgate] == X)
				{
					FANoutput[nextgate] = v1;
					/*
					printf("\tmandatory signal assignment: gut=%d val=%d", nextgate->index, v1);
					*/
					send = max((int)level[nextgate], send);
					push(g_stack, nextgate);
					schedule_input(gut, i);
				}
			}
		}
	}


	return(send);
}


template <class vertex>
int fan2(graph<vertex>& GA, int nog, int maxdpi, int npi, int npo, FAULT cf, int maxbacktrack, int *nbacktrack)
{
	register int i;
	register int gut, g;
	int last_dpi;
	int last_Dfrontier;
	bool backward_flag, backtrace_flag, fault_propagated_to_po;
	bool Dfrontier_changed, done;
	int state;
	FAULT pOriginFault = {-1,-1,0};

	*nbacktrack = 0;
	reset(done);
	reset(backward_flag);
	set(backtrace_flag);
	reset(fault_propagated_to_po);
	last_Dfrontier = 0;

	gut = cf.gateId;

	initNetAndFreach(GA, nog, gut, maxdpi);		/* initializaiton */

	if (cf.faultLine != OUTFAULT)
	{
		gut = inNbr(gut, cf.faultLine);
	}
	if (is_free(gut))
	{
		/* box 1 */
		// pOriginFault = (FAULT)malloc(sizeof(FAULT));
		pOriginFault.gateId = cf.gateId;
		pOriginFault.faultLine= cf.faultLine;
		pOriginFault.ftype = cf.ftype;
		propFault2Headline(GA, cf);
		last_dpi = 0;
	}
	else
	{
		last_dpi = refFaultyGateOutput(GA, cf);
	}
	if (last_dpi == (-1))
	{
		return(NO_TEST);
	}

	gut = cf.gateId;

	push(g_DfrontierStack, gut);
	i = unique_sensitize(GA, gut, gut);
	if ((last_dpi = max(i, last_dpi)) > 0)
	{
		set(backward_flag);
	}
	state = 93;

	/* main loop of fan algorithm */
	while (done == false)
	{
		switch (state)
		{
		case 93:
			/* box 3,4,5,6 */

			if (!implyForwardAndBackward(GA, maxdpi, backward_flag, last_dpi, cf))
			{
				/* box 3 */
				state = 98;
				break;
			}
			if (FANoutput[gut] == ZERO || FANoutput[gut] == ONE)
			{
				state = 98; break;
			}

			/* update unjustified lines and delete duplicated lines 
																	final_obj should be empty */
			for (i = g_unjustStack.last; i >= 0; i--)
			{
				g = g_unjustStack.list[i];
				if (is_justified(g))
				{
					delete(g_unjustStack, i);
				}
				else
				{
					set(changed[g]);
					push(g_finalObjStack, g);
				}
			}
			while (!is_empty(g_finalObjStack))
				reset(changed[pop(g_finalObjStack)]);

			/* check for backtrace */
			for (i = g_initObjStack.last; i >= 0; i--)
				if (is_justified(g_initObjStack.list[i]))
				{
					delete(g_initObjStack, i);
				}

			reset(fault_propagated_to_po);
			for (i = 0; i < npo; i++)
				if (FANoutput[index_PO[i]] == D || FANoutput[index_PO[i]] == DBAR)
				{
					set(fault_propagated_to_po);
					break;
				}

			if (last_Dfrontier != 0)
			{
				if (FANoutput[last_Dfrontier] == X)
				{
					reset(Dfrontier_changed);
				}
				else
				{
					set(Dfrontier_changed);
				}
			}
			else
			{
				set(Dfrontier_changed);
			}

			if (is_empty(g_initObjStack) && Dfrontier_changed)	/* box 4, 4-1 */
			{
				set(backtrace_flag);
			}

			if (fault_propagated_to_po)
			{
				state = 99;				/* box 4-3 */
				for (i = g_unjustStack.last; i >= 0; i--)
					if (is_unjustified(g_unjustStack.list[i]) && is_bound(g_unjustStack.list[i]))
					{
						state = 97;
						break;
					}
			}
			else
			{
				/* box 5 */
				/* update Dfrontier */
				if (!is_empty(g_DfrontierStack))
				{
					update_Dfrontier(GA);
				}
				for (i = gut; i < nog; i++)
					if (xpath[i] == 2)
					{
						xpath[i] = 1;
					}
				for (i = g_DfrontierStack.last; i >= 0; i--)
					if (!Xpath(GA, g_DfrontierStack.list[i]))
					{
						delete(g_DfrontierStack, i);
					}

				if (is_empty(g_DfrontierStack))
				{
					state = 98;
				}
				/*when dfrontier is not zero*/
				else
				{
					/* box 6 */
					if (dy_id >= INT_MAX - 3)
					{
						for (i = 0; i < nog; i++)
							freach1[i] = 0;
						dy_id = 0;
					}
					if ((last_dpi = dynamic_unique_sensitize(GA, g_DfrontierStack.list, g_DfrontierStack.last, maxdpi, dom_array, gut)) > 0)
					{
						set(backward_flag);
						state = 93;
					}
					else if (last_dpi == 0)
					{
						state = 93;
					}
					else
					{
						state = 97;
					}
				}
			}
			break;

		case 97:
			/* box 7 */

			find_final_objective(GA, &backtrace_flag, fault_propagated_to_po, nog, &last_Dfrontier);

			while (!is_empty(g_finalObjStack))
			{
				g = pop(g_finalObjStack);
				if (numzero[g] > numone[g])
				{
					FANoutput[g] = ZERO;
				}
				else
				{
					FANoutput[g] = ONE;
				}
				(g_tree.last)++;
				current_node.gate_id = g;
				reset(current_node.flag);
				push(g_stack, g);
				if (is_head(g))
				{
					set(changed[g]);
				}
				else
				{
					pushevent(g);
				}		/**** study ****/
				schedule_output(g);
				current_node.pstack = g_stack.last;
			}

			reset(backward_flag);
			state = 93;
			break;

		case 98:
			/* box 8 */

			if (is_flagged(current_node))
			{
				(*nbacktrack)++;
			}
			for (i = 0; i < maxdpi; i++)
				while (!is_empty(g_pEventListStack[i]))
					reset(changed[pop(g_pEventListStack[i])]);
			if (*nbacktrack > maxbacktrack)
			{
				state = OVER_BACKTRACK;
			}
			else if (backtrack(GA, gut, &last_dpi, nog))
			{
				if (last_dpi > 0)
				{
					set(backward_flag);
				}
				else
				{
					reset(backward_flag);
				}
				set(backtrace_flag);
				last_Dfrontier = 0;
				state = 93;
			}
			else
			{
				state = NO_TEST;
			}
			break;

		case 99:
			justify_free_lines(GA, npi, pOriginFault, cf);
			state = TEST_FOUND;
		default:
			set(done);
		}
	}

	if (pOriginFault.gateId != 0)
	{
		cf.gateId = pOriginFault.gateId;
		cf.faultLine = pOriginFault.faultLine;
		cf.ftype = pOriginFault.ftype;
		// free((char*)pOriginFault);
	}
	return(state);
}
