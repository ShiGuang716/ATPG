 

#define inCnt(idx) (GA.V[idx].getInDegree())
#define outCnt(idx) (GA.V[idx].getOutDegree())
#define inNbr(i,j) (GA.V[i].getInNeighbor(j))
#define outNbr(i,j) (GA.V[i].getOutNeighbor(j))
#define is_empty(s) (s.last<0)
#define pop(s) s.list[(s.last)--]

/* macros for parallel gate evaluation */
/* one input gate */
//pgate_eval1
#define evalGateWithInputs1_1(pGateId, iValue) \
	iValue = (type[pGateId] == STEM)? output[inNbr(pGateId,0)] : \
    (type[pGateId] == NOT || type[pGateId] == NAND || type[pGateId] == NOR)? ~output[inNbr(pGateId,0)] : \
    output[inNbr(pGateId,0)]


/* two input gates */
//pgate_eval2
#define evalGateWithInputs2_1(pGateId, iValue) \
	switch (type[pGateId]) \
	{ \
	   case AND: iValue = output[inNbr(pGateId,0)] & output[inNbr(pGateId,1)]; break; \
	   case NAND: iValue = ~(output[inNbr(pGateId,0)] & output[inNbr(pGateId,1)]); break; \
	   case OR: iValue = output[inNbr(pGateId,0)] | output[inNbr(pGateId,1)]; break; \
	   case NOR: iValue = ~(output[inNbr(pGateId,0)] | output[inNbr(pGateId,1)]); break; \
	   case XOR: iValue = output[inNbr(pGateId,0)] ^ output[inNbr(pGateId,1)]; break; \
	   case XNOR: iValue = ~(output[inNbr(pGateId,0)] ^ output[inNbr(pGateId,1)]); break; \
	}

/* 3-input gate */
//pgate_eval3
#define evalGateWithInputs3_1(pGateId,iValue) \
	switch (type[pGateId]) \
	{ \
	   case AND: iValue = output[inNbr(pGateId,0)] & output[inNbr(pGateId,1)] & output[inNbr(pGateId,2)]; break; \
	   case NAND: iValue = ~(output[inNbr(pGateId,0)] & output[inNbr(pGateId,1)] & output[inNbr(pGateId,2)]); break; \
	   case OR: iValue = output[inNbr(pGateId,0)] | output[inNbr(pGateId,1)] | output[inNbr(pGateId,2)]; break; \
	   default: iValue = ~(output[inNbr(pGateId,0)] | output[inNbr(pGateId,1)] | output[inNbr(pGateId,2)]); \
	}


/* more than 4-inputs */
//pgate_eval4
#define evalGateWithInputs4_1(pGateId, iValue) \
switch(type[pGateId]) { \
   case AND: \
	  iValue=output[inNbr(pGateId,0)]&output[inNbr(pGateId,1)]&output[inNbr(pGateId,2)]&output[inNbr(pGateId,3)]; \
	  for(int i=4;i<inCnt(pGateId);i++) iValue&=output[inNbr(pGateId,i)]; break; \
   case NAND: \
	  iValue=output[inNbr(pGateId,0)]&output[inNbr(pGateId,1)]&output[inNbr(pGateId,2)]&output[inNbr(pGateId,3)]; \
	  for(int i=4;i<inCnt(pGateId);i++) iValue&=output[inNbr(pGateId,i)]; \
	  iValue=~iValue; break; \
   case OR: \
	  iValue=output[inNbr(pGateId,0)]|output[inNbr(pGateId,1)]|output[inNbr(pGateId,2)]|output[inNbr(pGateId,3)]; \
	  for(int i=4;i<inCnt(pGateId);i++) iValue|=output[inNbr(pGateId,i)]; break; \
   default: \
	  iValue=output[inNbr(pGateId,0)]|output[inNbr(pGateId,1)]|output[inNbr(pGateId,2)]|output[inNbr(pGateId,3)]; \
	  for(int i=4;i<inCnt(pGateId);i++) iValue|=output[inNbr(pGateId,i)]; \
	  iValue=~iValue; \
}

extern int *output;
extern int *output1;
extern long noPI;
// extern bool *changed;
// extern bool *freach;
extern gtype *type;
// extern STACK *levelStack;
extern std::vector<vector<int>> levelStack;



template <class vertex>
void GoodSim(graph<vertex> &GA, unsigned int max_level) //GoodSim
{    
    int n = GA.n;
	int iValue = 0;
    // cout<< max_level <<endl;
    for(int i = noPI; i < n; i++)
    {
            int pGateId = i;
            int indegree = inCnt(pGateId);
		    switch (indegree)
		    {
		    case 1:
			    evalGateWithInputs1_1(pGateId, output[pGateId]); 
                break;
		    case 2:
			    evalGateWithInputs2_1(pGateId, output[pGateId]); 
                break;
		    case 3:
			    evalGateWithInputs3_1(pGateId, output[pGateId]); 
                break;
		    default:
			    evalGateWithInputs4_1(pGateId, iValue);
			    output[pGateId] = iValue;
		    }

		    freach[pGateId] = false;
		    changed[pGateId] = false;
        
    }

	// for(int i = 1; i < max_level +2; i++)
    // {
    //     for(int j=0;j< levelStack[i].size();j++)
    //     {
    //         int pGateId = levelStack[i][j];
    //         int indegree = inCnt(pGateId);
	// 	    switch (indegree)
	// 	    {
	// 	    case 1:
	// 		    evalGateWithInputs1(pGateId, output[pGateId]); 
    //             break;
	// 	    case 2:
	// 		    evalGateWithInputs2(pGateId, output[pGateId]); 
    //             break;
	// 	    case 3:
	// 		    evalGateWithInputs3(pGateId, output[pGateId]); 
    //             break;
	// 	    default:
	// 		    evalGateWithInputs4(pGateId, iValue);
	// 		    output[pGateId] = iValue;
	// 	    }

	// 	    freach[pGateId] = false;
	// 	    changed[pGateId] = false;
    //     }
    // }
}

