#include <string.h>
#include <cstring> 
#include <vector>
#include <stack>
#include <mutex>
#include <forward_list>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include "ligra.h"
#include "define.h"
#include "truthtable.h"
// #include "preprocess.h"
// #include "implication.h"
// #include "backtrace.h"
// #include "justification.h"
// #include "backtrack.h"
// #include "FAN.h"
// #include "FAN1.h"
#include "FAN2.h"
#include "GOODSIM2.h"
#include "FAULTSIM.h"

using namespace std;

std::vector<char *> patterns_pool; 
std::vector<long> index_PI;
std::vector<long> index_PO;
std::vector<long> index_tie0;
std::vector<long> index_head;
std::vector<int> curObj;
std::vector<int> initObj;
std::vector<int> fanObj;
std::vector<int> headObj;
std::vector<int> finalObj;
// std::vector<long> g_stack;
// std::stack<TREETYPE> g_tree;
std::vector<int> Dfrontier;
std::vector<int> unjustify;

int* FANoutput;
bool* active_set;
int* changed;
int* ltype = NULL;
unsigned int* dpo = NULL;
int* fosIndex = NULL;
int *freach = NULL;
int *freach1 = NULL;
int *observe = NULL;
int* cobserve = NULL;
long n;
long noPI = 0, noPO = 0, noTIE_0 = 0, noSTEM = 0;
long num_of_faults;
long detected_faults;
long noRedundant_faults = 0;
int generated_patterns;
int eff_patterns;
int test_patterns;
gtype *type = NULL;
unsigned int *level;

int *output;
int *output1;
unsigned char *cnt_in;
unsigned short *fault_list;
int *xpath;
int iState;
timer atpg_timer;
timer fan_timer;
timer faultsim_timer;
timer goodsim_timer;
double fan_time = 0;
double faultsim_time = 0;
double goodsim_time = 0;
double atpg_time = 0;
unsigned int max_level = 0;
LINKTYPE **u_path = NULL;
int *numzero;
int *numone;
int *cont0;
int *cont1;
int iNoBackTrack = 0;
int g_iMaxBackTrack1 = 10;		/* maximum backtracking of FAN1 */
int g_iMaxBackTrack2 = 0;		/* maximum backtracking of FAN2 */
int iNoStemGates;
int *pStemGates;
std::vector<vector<FAULT> > dfault; // save the faults that contribute for observe
std::vector<int> gate_index;



//读取gate_type文件用到（不重要）
void split(std::string s, std::string delimiter, std::vector<std::string> &res)
{
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	res.clear();
   
	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
	{
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	} 
     
	res.push_back(s.substr(pos_start));
}


//多路选择器（不重要）
long binarySearch(MUXINFO* muxinfo, long n, long target) {
	int low = 0, high = n, middle = 0;
	while(low < high) {
		middle = (low + high)/2;
		if(target == muxinfo[middle].mux_dst_gate) {
			return middle;
		} else if(target < muxinfo[middle].mux_dst_gate) {
			high = middle;
		} else if(target > muxinfo[middle].mux_dst_gate) {
			low = middle + 1;
		}
	}
	printf("binarySearch Error!");
	exit(-1);
}
//（不重要）
void insert_active_gate(forward_list<long> &fl,const long &insert_gate_id){
    if (fl.empty()) {
		fl.push_front(insert_gate_id);
		return;
	}
	auto prev=fl.before_begin();
    auto curr=fl.begin();
    while(1){
		long curr_gate_id = *curr;
        if (level[insert_gate_id] < level[curr_gate_id]){
            fl.insert_after(prev, insert_gate_id);
            return;
        } else {
            prev=curr;
            ++curr;
            if(curr==fl.end()) break;
        }
    }
    fl.insert_after(prev, insert_gate_id);
}


//确定各门的逻辑级
template <class vertex>
struct LEVELIZE
{
	unsigned int *level;
	graph<vertex> &GA;
	LEVELIZE(unsigned int *_level, graph<vertex> &_GA) : level(_level), GA(_GA) {}    
  
	inline bool updateAtomic(long s, long d)
	{ 
		bool done = false;
		volatile unsigned int newV, oldV;
		while (!done)
		{ 
			oldV = level[d];
			newV = level[s] + 1;
			if (newV > oldV)
			{
				done = CAS(&level[d], oldV, newV);
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




// 逻辑仿真
template <class vertex>
struct GOOD_SIM {
	unsigned int* logic;
	unsigned char* cnt_in;
	gtype* type;
	MUXINFO* muxinfo;
	long mux_gate_num;
	graph<vertex>& GA;
	GOOD_SIM(unsigned int* _logic, unsigned char* _cnt_in, gtype* _type, graph<vertex>& _GA, long _mux_gate_num, MUXINFO* _muxinfo) : 
	logic(_logic), cnt_in(_cnt_in), type(_type), GA(_GA), muxinfo(_muxinfo), mux_gate_num(_mux_gate_num) { }
	inline bool updateAtomic (long s, long d) {
	
		switch(type[d]) {
			case AND:
			case NAND:
			writeAnd(&logic[d], logic[s]);break;
			case OR:
			case NOR:
			writeOr(&logic[d], logic[s]);break;
			case XOR:
			case XNOR:
			writeXor(&logic[d], logic[s]);break;
			case MUX:
			{
				long mux_index = binarySearch(muxinfo, mux_gate_num, d);
					
				long s0, a, b;
				s0	= muxinfo[mux_index].mux_s0_gate;
				a	= muxinfo[mux_index].mux_a_gate;
				b	= muxinfo[mux_index].mux_b_gate;
				if (s != s0 && s != a && s != b) {
					printf("muxinfo error!");
					exit(-1);
				}

				int mux_count = writeAdd(&muxinfo[mux_index].count, (int)1);
				if (mux_count == 2) {
					logic[d] = logic[a] & ~logic[s0] | logic[b] & logic[s0];
				}
				break;
			}
			default:
			logic[d] = logic[s];break;
		}

		volatile unsigned char newV, oldV;
		do {oldV = cnt_in[d]; newV = oldV + 1;}
		while (!CAS(&cnt_in[d], oldV, newV));
		
		if (newV == (unsigned char)GA.V[d].getInDegree()) {
			
			
			
			if (type[d] == NAND || type[d] == NOR || type[d] == XNOR || type[d] == NOT) writeNotAll(&logic[d]);
			
			
			return true;
		} else {
			return false;
		}
	}
	inline bool update (long s, long d) {
		return updateAtomic(s,d);
	}
	inline bool cond (long d) {
		
		return (cnt_in[d] < GA.V[d].getInDegree());
	}

	inline long binarySearch(MUXINFO* muxinfo, long n, long target) {
		int low = 0, high = n, middle = 0;
		while(low < high) {
			middle = (low + high)/2;
			if(target == muxinfo[middle].mux_dst_gate) {
				return middle;
			} else if(target < muxinfo[middle].mux_dst_gate) {
				high = middle;
			} else if(target > muxinfo[middle].mux_dst_gate) {
				low = middle + 1;
			}
		}
		printf("binarySearch Error!");
		exit(-1);
	}
};


template <class vertex>
void fault_sim(graph<vertex> &GA)
{
    // good sim
    Log("Good Sim Start!\n");
	goodsim_timer.next();
    GoodSim(GA,max_level);
	goodsim_time += goodsim_timer.next();

	Log("Good Sim Finished!\n");
	faultsim_timer.next();
	Log("Fault Sim Start!\n");
	for (long i=0; i<n; i++) {
		if (fault_list[i] != 0) {
			// 先处理SA0 fault
			for (unsigned int j=0; j<8; j++) {
				if ((fault_list[i] & add_SA_Fault[0][j]) != 0) {
					memcpy(output1, output, sizeof(int) * n);
					long faulty_gate;
					if(j < 7){
						faulty_gate = GA.V[i].getInNeighbor(j);
						output1[faulty_gate] = ALL0;
					}
					else{
						faulty_gate = i; //output fault
						output1[faulty_gate] = ALL0;
					}

					forward_list<long> active_gates;
					unordered_set<long> visited_gates;
                    if(j < 7){
                        active_gates.push_front(i);
					    visited_gates.insert(i);
                    }
					else{
                        for(long k=0; k<GA.V[i].getOutDegree(); k++) {
							long temp = GA.V[i].getOutNeighbor(k);
							if (visited_gates.find(temp) == visited_gates.end()) {
								insert_active_gate(active_gates, temp);
								visited_gates.insert(temp);
							}
						}
                    }
					
					while (!active_gates.empty()) {
						auto curr_iter = active_gates.begin();
						long curr_gate_id = *curr_iter;
						active_gates.pop_front();

						output1[curr_gate_id] = output1[GA.V[curr_gate_id].getInNeighbor(0)];
						for (int k=1; k<GA.V[curr_gate_id].getInDegree(); k++) {
							switch(type[curr_gate_id]) {
								case AND:
								case NAND:
								output1[curr_gate_id] = output1[curr_gate_id] & output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case OR:
								case NOR:
								output1[curr_gate_id] = output1[curr_gate_id] | output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case XOR:
								case XNOR:
								output1[curr_gate_id] = output1[curr_gate_id] ^ output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case MUX:
								{
									// i = GA.V[curr_gate_id].getInDegree();
									// long mux_index = binarySearch(muxinfo, mux_gate_num, curr_gate_id);

									// long s0, a, b;
									// s0	= muxinfo[mux_index].mux_s0_gate;
									// a	= muxinfo[mux_index].mux_a_gate;
									// b	= muxinfo[mux_index].mux_b_gate;

									// output1[curr_gate_id] = output1[a] & !output1[s0] | output1[b] & output1[s0];
                    				break;
								}
								default:break;
							}
						}
						if (type[curr_gate_id] == NAND || type[curr_gate_id] == NOR || 
							type[curr_gate_id] == XNOR || type[curr_gate_id] == NOT) 
						{
							output1[curr_gate_id] = ~output1[curr_gate_id];
						}
						if ( (j < 7) && (curr_gate_id == i) && (GA.V[faulty_gate].getOutDegree() > 1) ){
							output1[faulty_gate] = output[faulty_gate];
						}
						if (output1[curr_gate_id] != output[curr_gate_id]) {
							if (type[curr_gate_id] == PO) {
								detected_faults++;
								if (detected_faults%100==0) cout << "detected_faults: " << detected_faults << endl;
								fault_list[i] = fault_list[i] & del_SA_Fault[0][j];
								break;
							} else {
								for(long k=0; k<GA.V[curr_gate_id].getOutDegree(); k++) {
									long temp = GA.V[curr_gate_id].getOutNeighbor(k);
									if (visited_gates.find(temp) == visited_gates.end()) {
										insert_active_gate(active_gates, GA.V[curr_gate_id].getOutNeighbor(k));
										visited_gates.insert(temp);
									}
								}
							}
						}
					// while end
					}
				}
			}
			// 再处理SA1 fault
			for (unsigned int j=0; j<8; j++) {
				if ((fault_list[i] & add_SA_Fault[1][j]) != 0) {
					// bool debug = false;
					// if (FANoutput[index_PI[0]] == 1 && FANoutput[index_PI[1]] == 2
					// 	&& FANoutput[index_PI[2]] == 0 && FANoutput[index_PI[3]] == 1
					// 	&& FANoutput[index_PI[4]] == 0) {
					// 		if (i==7 && j==1) {
					// 			debug = true;
					// 		}
					// 	}										
					memcpy(output1, output, sizeof(int) * n);
					long faulty_gate;
					if(j < 7){
						faulty_gate = GA.V[i].getInNeighbor(j);
						output1[faulty_gate] = ALL1;
					}
					else{
						faulty_gate = i; //output fault
						output1[faulty_gate] = ALL1;
					}

					forward_list<long> active_gates;
					unordered_set<long> visited_gates;
					if(j < 7){
                        active_gates.push_front(i);
					    visited_gates.insert(i);
                    }
					else{
                        for(long k=0; k<GA.V[i].getOutDegree(); k++) {
							long temp = GA.V[i].getOutNeighbor(k);
							if (visited_gates.find(temp) == visited_gates.end()) {
								insert_active_gate(active_gates, temp);
								visited_gates.insert(temp);
							}
						}
                    }

					while (!active_gates.empty()) {
						auto curr_iter = active_gates.begin();
						long curr_gate_id = *curr_iter;
						active_gates.pop_front();

						output1[curr_gate_id] = output1[GA.V[curr_gate_id].getInNeighbor(0)];
						for (int k=1; k<GA.V[curr_gate_id].getInDegree(); k++) {
							switch(type[curr_gate_id]) {
								case AND:
								case NAND:
								output1[curr_gate_id] = output1[curr_gate_id] & output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case OR:
								case NOR:
								output1[curr_gate_id] = output1[curr_gate_id] | output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case XOR:
								case XNOR:
								output1[curr_gate_id] = output1[curr_gate_id] ^ output1[GA.V[curr_gate_id].getInNeighbor(k)];
								break;
								case MUX:
								{
									// i = GA.V[curr_gate_id].getInDegree();
									// long mux_index = binarySearch(muxinfo, mux_gate_num, curr_gate_id);

									// long s0, a, b;
									// s0	= muxinfo[mux_index].mux_s0_gate;
									// a	= muxinfo[mux_index].mux_a_gate;
									// b	= muxinfo[mux_index].mux_b_gate;

									// output1[curr_gate_id] = output1[a] & !output1[s0] | output1[b] & output1[s0];
                    				break;
								}
								default:break;
							}
						}
						if (type[curr_gate_id] == NAND || type[curr_gate_id] == NOR || 
							type[curr_gate_id] == XNOR || type[curr_gate_id] == NOT) 
						{
							output1[curr_gate_id] = ~output1[curr_gate_id];
						}
						// if (debug) {
						// 	printf("debug: output[%ld]==%d, output1[%ld]==%d\n",curr_gate_id,output[curr_gate_id],curr_gate_id,output1[curr_gate_id]);
						// 	if(curr_gate_id == 11) {
						// 		printf("gate %ld:\n", curr_gate_id);
						// 		for(int k=0; k<GA.V[curr_gate_id].getInDegree(); k++) {
						// 			long temp = GA.V[curr_gate_id].getInNeighbor(k);
						// 			printf("incoming edge: %ld, output[%ld]=%d, output1[%ld]=%d\n",temp,temp,output[temp],temp,output1[temp]);
						// 		}
						// 	}
						// }
						if ( (j < 7) && (curr_gate_id == i) && (GA.V[faulty_gate].getOutDegree() > 1) ){
							output1[faulty_gate] = output[faulty_gate];
						}
						if (output1[curr_gate_id] != output[curr_gate_id]) {
							if (type[curr_gate_id] == PO) {
								detected_faults++;
								if (detected_faults%100==0) cout << "detected_faults: " << detected_faults << endl;
								fault_list[i] = fault_list[i] & del_SA_Fault[1][j];
								break;
							} else {
								for(int k=0; k<GA.V[curr_gate_id].getOutDegree(); k++) {
									long temp = GA.V[curr_gate_id].getOutNeighbor(k);
									if (visited_gates.find(temp) == visited_gates.end()) {
										insert_active_gate(active_gates, GA.V[curr_gate_id].getOutNeighbor(k));
										visited_gates.insert(temp);
									}
								}
							}
						}
					// while end
					}
				}
			}
		}
	}
	faultsim_time += faultsim_timer.next();
	Log("Fault Sim Finished!\n");

}


//读取gate_type文件和csr文件
template <class vertex>
void Compute(graph<vertex> &GA, commandLine P)
{ 	
	// srand((unsigned)time(NULL));
	// default_random_engine e(static_cast<unsigned int>(time(nullptr)));
	// std::uniform_int_distribution<> u(0, 1);
	char *typeFilename = P.getOptionValue("-t");
	char *muxinfoFilename = P.getOptionValue("-m");
    n = GA.n;
	long mux_gate_num = 0;
	std::ifstream in(typeFilename);
	std::ofstream out("output.txt", std::ios::binary);
	std::string line;
	std::vector<std::string> container;
	long count = 0;

    type = new gtype[n];
	FILE *f = fopen("log.txt", "w");
	MUXINFO *muxinfo;
    // 读取gate type文件
	if (in.is_open())
	{
		while (std::getline(in, line))
		{
			split(line, " ", container);
			long gate_id = (long)std::stol(container[0]);
			int gate_type = (int)std::stoi(container[1]);
			type[gate_id] = gate_type;

			if (gate_type == PI)
			{
				index_PI.push_back(gate_id);
				noPI++;
			}
			else if (gate_type == PO)
			{
				index_PO.push_back(gate_id);
				noPO++;
			}  
			else if (gate_type == TIE_0)
			{
				index_tie0.push_back(gate_id);
				noTIE_0++;
			}
			else if (gate_type == STEM)
			{       
				noSTEM++;
			}
			count++;
		}
	}
	in.close();
    printf("number of PI: %ld\n", noPI);
	printf("number of TIE_0: %ld\n", noTIE_0);
	printf("number of PO: %ld\n", noPO);
	printf("number of STEM: %ld\n", noSTEM);
	printf("count:%ld n:%ld\n", count, n);
    if (count != n)
	{
		printf("type file error!\n");
		exit(-1);
	}
	printf("type file loaded!\n");

    
	// // 读取mux info文件
	// std::ifstream mux_in(muxinfoFilename);
	// long mux_count = 0;   
	// if (mux_in.is_open())
	// {
	// 	std::getline(mux_in, line);
	// 	mux_gate_num = (long)std::stol(line);
	// 	muxinfo = new MUXINFO[mux_gate_num];

	// 	while (std::getline(mux_in, line))
	// 	{
	// 		split(line, " ", container);
	// 		muxinfo[mux_count].mux_dst_gate = (long)std::stol(container[0]);
	// 		muxinfo[mux_count].mux_s0_gate = (long)std::stol(container[1]);
	// 		muxinfo[mux_count].mux_a_gate = (long)std::stol(container[2]);
	// 		muxinfo[mux_count].mux_b_gate = (long)std::stol(container[3]);
         
	// 		if (mux_count > 0 && muxinfo[mux_count].mux_dst_gate < muxinfo[mux_count - 1].mux_dst_gate)
	// 		{
	// 			printf("muxinfo should be sorted!");
	// 			exit(-1);
	// 		}
	// 		mux_count++;
	// 	}
	// 	if (mux_count != mux_gate_num)
	// 	{
	// 		printf("muxinfo file error!\n");
	// 		exit(-1);
	// 	}
	// }
	// mux_in.close();
	// if (mux_count != mux_gate_num)
	// {
	// 	printf("mux gate num error!\n");
	// 	exit(-1);
	// }
	// printf("muxinfo file loaded!\n");

    // levelize
	//unsigned int *level = new unsigned int[n];
	level = new unsigned int[n];

	bool *frontier = newA(bool, n);

	printf("levelize start!\n");

	parallel_for(long i = 0; i < n; i++)
	{
		level[i] = 0;
		if (type[i] == PI || type[i] == TIE_0)
		{
			frontier[i] = 1;
		}
		else 
		{
			frontier[i] = 0;
		} 
	} 
	vertexSubset Frontier(n, noPI + noTIE_0, frontier);
	while (!Frontier.isEmpty())
	{
		vertexSubset out = edgeMap(GA, Frontier, LEVELIZE<vertex>(level, GA));
		Frontier.del();
		Frontier = out;
	}
	Frontier.del();
	for (long i = 0; i < n; i++)
	{
		if (type[i] == PO)
		{
			if (level[i] > max_level)
			{
				max_level = level[i];
			}
		} 
	}  
	printf("max level: %d\n", max_level);
	// for (long i = 0; i < n; i++)
	// {
	// 	if (type[i] == PO)
	// 	{
	// 		level[i] = max_level;
	// 	} 
	// }
	printf("levelize finished!\n");

    // 生成故障列表
	fault_list = new unsigned short[n];
    output = new int[n];
    output1 = new int[n];
    cnt_in = new unsigned char[n];

	num_of_faults = 0;
	detected_faults = 0;

	printf("fault list generating start!\n");
	memset(fault_list, 0, sizeof(unsigned short) * n);

	parallel_for(long i = 0; i < n; i++)
	{
		//add input SA0 & SA1 fault
		int inputs = GA.V[i].getInDegree();
		if (inputs > 1)
		{ 
			if (type[i] == AND || type[i] == NAND || type[i] > PI)
			{
				for (int j = 0; j < inputs; j++)
				{
					fault_list[i] = fault_list[i] | add_SA_Fault[1][j];
					writeAdd(&num_of_faults, (long)1);
				}
			} 
			if (type[i] == OR || type[i] == NOR || type[i] > PI)
			{
				for (int j = 0; j < inputs; j++)
				{
					fault_list[i] = fault_list[i] | add_SA_Fault[0][j];
					writeAdd(&num_of_faults, (long)1);
				}
			}
		}
		if (type[i] == PO)
		{
			fault_list[i] = fault_list[i] | add_SA_Fault[0][0] | add_SA_Fault[1][0];
			writeAdd(&num_of_faults, (long)2);
		}
		//add output SA0 & SA1 fault for fanout gate 
		int outputs = GA.V[i].getOutDegree();
		if(outputs > 1)
		{
			fault_list[i] = fault_list[i] | add_SA_Fault[0][7] | add_SA_Fault[1][7];
			writeAdd(&num_of_faults, (long)2);
		}
	}

	printf("num of faults: %ld\n", num_of_faults);
	printf("fault list generating finished!\n");
	
	gate_index.resize(n,0);
	for(int i = 0; i < n; i++){
		gate_index[i] = i;
	}
	sort(gate_index.begin(), gate_index.end(), [GA](int a, int b){
		return GA.V[a].getInDegree() > GA.V[b].getInDegree();
	});

    // ATPG开始
    // FAN算法预处理（概要视频08：00）
    Log("Preprocess start!");
	preprocess(GA, n, noPI, max_level + 2);
	Log("Preprocess end.");

// FAN算法开始
#ifdef FAN
	Log("FAN start!");
	// FANoutput = newA(char, n);
	// freach = newA(bool, n);
	// active_set = newA(bool, n);
	// changed = newA(bool, n);
#ifdef ONEBYONE
	int num = 0;
	for (long i = 0; i < n; i++)
	{
		if (fault_list[i] != 0)
		{
			// 先处理SA0
			for (unsigned int j = 0; j < 8; j++)
			{      
				if ((fault_list[i] & add_SA_Fault[0][j]) == 0)
					continue;
				num++;
				if (num == NUM) 
				{
					Log("SA0 fault.");
					Log("Fault gate: %ld, Fault line: %d, gate type: %d.", i, j, type[i]);
					printf("Fault gate: %ld, Fault line: %d, gate type: %d, sa0.\n", i, j, type[i]);
                    FAULT pCurrentFault;
                    pCurrentFault.gateId=i;
                    pCurrentFault.faultLine = j;
                    pCurrentFault.ftype = SA_0;
					int rt = fan1(GA, n, max_level+2,  noPI, noPO, pCurrentFault, g_iMaxBackTrack1, &iNoBackTrack);
					printf("fan end\n");

					for(int k = 0; k < noPI; k++){
						printf("%d ", FANoutput[k]);
					}
					printf("\n");					
				}
			}  
			// 再处理SA1
			for (unsigned int j = 7; j < 8; j++)
			{     
				if ((fault_list[i] & add_SA_Fault[1][j]) == 0)
					continue;
				num++; 
				if (num == NUM)             
				{  
					Log("SA1 fault.");
					Log("Fault gate: %ld, Fault line: %d, gate type: %d.", i, j, type[i]);
					printf("Fault gate: %ld, Fault line: %d, gate type: %d, sa1.\n", i, j, type[i]);
                    FAULT pCurrentFault;
                    pCurrentFault.gateId=i;
                    pCurrentFault.faultLine = j;
                    pCurrentFault.ftype = SA_1;
					int rt = fan1(GA, n, max_level+2,  noPI, noPO, pCurrentFault, g_iMaxBackTrack1, &iNoBackTrack);
					printf("fan end\n");

					for(int k = 0; k < noPI; k++){
						printf("%d ", FANoutput[k]);
					}
					printf("\n");						
				}
			}   
		}
	}   
#else 
	generated_patterns = 0;
	int rt;
	long fault_in_list = num_of_faults;
	noRedundant_faults = 0;
	atpg_timer.start();
	goodsim_timer.start();
	faultsim_timer.start();
	fan_timer.start();
	split(typeFilename, "_gatetype", container);
	string cctFilename = container[0];
	string resultfile = cctFilename + "_result.txt";
	ofstream ofs;
    ofs.open(resultfile ,ios::out);
	printf("----------\t----------\t--------\t-----------\t--------\n");
	printf("#生成向量数\t#故障覆盖率\t#剩余故障\t本轮检测故障\t运行时间\n");
	for (long i = 0; i < n; i++)
	{ 
		if (fault_list[gate_index[i]] != 0)
		{
			for(int sa = 0; sa < 2; sa++)//每次取一个故障
			{
				for (unsigned int j = 0; j < 8; j++)
				{
					if ((fault_list[gate_index[i]] & add_SA_Fault[sa][j]) == 0)
						continue;
					Log("SA0 fault.");
					Log("Fault gate: %ld, Fault line: %d, gate type: %d.", i, j, type[i]);
					// printf("Fault gate: %ld, outcnt : %d ,Fault line: %d, FFR index : %d, sa0.\n", i, outCnt(i), j, fosIndex[i]);
					//if(generated_patterns>5000)
					FAULT pCurrentFault;
					pCurrentFault.gateId = gate_index[i];
					pCurrentFault.faultLine = j;//有一位是输出线其他都是输入
					pCurrentFault.ftype = sa;
					Log("fan start");
					
					//故障取出来先fan算法	测试生成每生成32个向量就合到一起，每一位代表一个向量，放到一个一维int数组，每一位都利用起来相当于32个向量压缩成了一个向量。				
					fan_timer.next();
					rt = fan1(GA, n, max_level+2,  noPI, noPO, pCurrentFault, g_iMaxBackTrack1, &iNoBackTrack);
					fan_time += fan_timer.next();
					
					Log("fan end");
					// rt = fan(GA, i, j, 0);
					//return;
					if (rt == TEST_FOUND) {
						detected_faults++;
						// if (detected_faults%100==0) cout << "detected_faults: " << detected_faults << endl;
						fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[sa][j];
						generated_patterns++;
						int bit = (generated_patterns % 32);
						// srand((unsigned)time(NULL));
						// default_random_engine e(static_cast<unsigned int>(time(nullptr)));
						// std::uniform_int_distribution<> u(0, 1);
						parallel_for(long k=0; k<noPI; k++) {
							if (type[k] == PI) {
								int temp = rand() % 2;
								// int temp = u(e);
								if (FANoutput[k] == 0 || (FANoutput[k] == X && temp ==0)) {
									clrbit(output1[k], bit);
								} else {
									setbit(output1[k], bit);
								}
							} else if (type[k] == TIE_0) {
								clrbit(output1[k], bit);
							} else if (type[k]==AND || type[k]==NAND) {
								setbit(output1[k], bit);
							}
						}

						// for(long k = 0; k < noPI; k++){
						// 	output1[k] & BITMASK[bit] != ALL0 ? ofs << 1 : ofs << 0;
						// }
						// ofs << endl;
						
						// memcpy(output, output1, sizeof(int) * n);

						if (bit == 0) {
							// fault_sim(GA);
							faultsim_timer.next();
							Fault0_Simulation(GA, n, max_level + 2, noPI, noPO);
							faultsim_time += faultsim_timer.next();
							// if((fault_list[i] & add_SA_Fault[0][j]) != 0) {
							// 	fault_list[i] = fault_list[i] & del_SA_Fault[0][j];						
							// 	printf("test pattern error!\n");
							// 	printf("SA1 Fault gate: %ld, Fault line: %d, gate type: %d.\n", i, j, type[i]);
							// 	exit(-1);
							// }
							long curr_detected_faults = fault_in_list - num_of_faults + detected_faults;
							fault_in_list = num_of_faults - detected_faults;
							double coverage=(double)detected_faults/(double)num_of_faults*100.00;
							atpg_time = goodsim_time + faultsim_time + fan_time;
							ofs << coverage << " " << atpg_time << endl;
							printf("----------\t----------\t--------\t-----------\t--------\n");
							printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns, coverage, fault_in_list, curr_detected_faults, atpg_time);
						}
						Log("FAN success!%d", generated_patterns);
					}
					else if(rt == NO_TEST){
						noRedundant_faults++;
						fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[sa][j];
						Log("Fan end, no pattern generated!");					
					}				
				} 
			}
		}
	}

    if (generated_patterns % 32 != 0) {
		faultsim_timer.next();
		Fault0_Simulation(GA, n, max_level + 2, noPI, noPO);
		faultsim_time += faultsim_timer.next();
		long curr_detected_faults = fault_in_list - num_of_faults + detected_faults;
		fault_in_list = num_of_faults - detected_faults;
		double coverage=(double)detected_faults/(double)num_of_faults*100.00;
		// printf("----------\t----------\t--------\t-----------\t--------\n");
		// printf("#生成向量数\t#故障覆盖率\t#剩余故障\t本轮检测故障\t运行时间\n");
		atpg_time = goodsim_time + faultsim_time + fan_time;
		ofs << coverage << " " << atpg_time << endl;
		printf("----------\t----------\t--------\t-----------\t--------\n");
		printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns, coverage, fault_in_list, curr_detected_faults, atpg_time);
	}

#endif 
	
	printf("1. Test pattern generation results\n");
	printf("   Number of test patterns                     : %d\n", generated_patterns);
	printf("   num of faults                               : %ld\n", num_of_faults);
	printf("   Number of detected faults                   : %ld\n", detected_faults);
	printf("   Number of identified redundant faults       : %ld\n", noRedundant_faults);
	printf("   Number of aborted faults                    : %ld\n", num_of_faults - detected_faults - noRedundant_faults);
	printf("   Fault coverage                              : %.3lf %%\n", (double) detected_faults / (double) num_of_faults * 100.0);
	printf("\n");

	printf("2. CPU time\n");
	printf("   good simulation time                        : %.3lf sec\n", goodsim_time);
	printf("   fault simulation time                       : %.2lf sec\n", faultsim_time - goodsim_time);
	printf("   fan time                                    : %.2lf sec\n", fan_time);
	printf("   total time                                  : %.2lf sec\n", atpg_time );
	printf("\n");


	free(FANoutput);   
	free(freach);
	free(active_set);
	free(changed);
#endif
	free(ltype);
	free(dpo);
	free(observe);
	free(cobserve);
	delete []type;
	delete []level;
	delete []fault_list;



}