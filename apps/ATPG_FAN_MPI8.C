#include <string.h>
#include <cstring> 
#include <vector>
#include <stack>
#include <mutex>
#include <forward_list>
#include <unordered_set>
#include <unordered_map>
#include "ligra.h"
#include "define.h"
#include "truthtable.h"
#include "FAN2.h"
#include "GOODSIM2.h"

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
int iNoStemGates;
int *pStemGates;
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
timer sychronize_timer;
double fan_time = 0;
double atpg_time = 0;
double faultsim_time = 0;
double goodsim_time = 0;
double sychronize_time = 0;
unsigned int max_level = 0;
LINKTYPE **u_path = NULL;
int *numzero;
int *numone;
int *cont0;
int *cont1;
int iNoBackTrack = 0;
int g_iMaxBackTrack1 = 100;		/* maximum backtracking of FAN */
int g_iMaxBackTrack2 = 0;		/* maximum backtracking of FAN */
// std::vector<long> gnet;
std::vector<int> gate_index;

int process_id;
int processes;
int local_bitsize;
int generated_patterns_global;
long detected_faults_global;
unsigned int finish_flag;
unsigned int finish_flag_global;

std::vector<vector<FAULT> > dfault; // save the faults that contribute for observe
// unsigned int *local_output;


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
	/**
	 * TODO()输入门一定在前几个吗？门id是拓扑排序好的吗？
	 * 
	 * 根据PI门的输入值模拟出其他门的正常值
	 * 输入：
	 * output[]
	 * 输出
	 * output[]
	*/
    GoodSim(GA,max_level);
	goodsim_time += goodsim_timer.next();

	Log("Good Sim Finished!\n");
	faultsim_timer.next();
	Log("Fault Sim Start!\n");
	//int under = ceil((float)(n/processes));
	//long i = process_id * under; i < under + process_id * under ; i++
	for (long i = process_id; i < n ; i+=processes) {
		if (fault_list[gate_index[i]] != 0) {
			// 先处理SA0 fault
			for (unsigned int j=0; j<8; j++) {
				/**
				 * TODO()
				 * active_gates与visited_gates需要每次重新初始化吗？
				 * 每次进行的是32个向量的故障模拟还是1个向量的故障模拟？
				 * 
				 * 对一个门的一个故障：
				 * 1、	若为输入导致的故障，将输入门的值置0
				 * 		若为输出导致的故障，将门的值置0
				 * 2、	若为输入导致的故障，将门id压入active_gates与visited_gates
				 * 		若为输出导致的故障，对门的每个出邻居，如果其不在visited_gates中，则压入active_gates与visited_gates
				 * 3、	模拟故障的传播：如果传播到PO门，则故障可检测。如果传播到其他门，则继续传播。
				*/
				if ((fault_list[gate_index[i]] & add_SA_Fault[0][j]) != 0) {
					memcpy(output1, output, sizeof(unsigned int) * n);
					long faulty_gate;
					if(j < 7){
						faulty_gate = GA.V[gate_index[i]].getInNeighbor(j);
						output1[faulty_gate] = ALL0;
					}
					else{
						faulty_gate = gate_index[i]; //output fault
						output1[faulty_gate] = ALL0;
					}

					forward_list<long> active_gates;
					unordered_set<long> visited_gates;
                    if(j < 7){
                        active_gates.push_front(gate_index[i]);
					    visited_gates.insert(gate_index[i]);
                    }
					else{
                        for(long k=0; k<GA.V[gate_index[i]].getOutDegree(); k++) {
							long temp = GA.V[gate_index[i]].getOutNeighbor(k);
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
						/**
						 * TODO()
						 * 这一步是复原吗？为什么有出度>1的判断？
						 * 有必要复原output1吗？
						*/
						if ( (j < 7) && (curr_gate_id == gate_index[i]) && (GA.V[faulty_gate].getOutDegree() > 1) ){
							output1[faulty_gate] = output[faulty_gate];
						}
						if (output1[curr_gate_id] != output[curr_gate_id]) {
							if (type[curr_gate_id] == PO) {
								detected_faults++;
								fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[0][j];
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
				if ((fault_list[gate_index[i]] & add_SA_Fault[1][j]) != 0) {
					// bool debug = false;
					// if (FANoutput[index_PI[0]] == 1 && FANoutput[index_PI[1]] == 2
					// 	&& FANoutput[index_PI[2]] == 0 && FANoutput[index_PI[3]] == 1
					// 	&& FANoutput[index_PI[4]] == 0) {
					// 		if (i==7 && j==1) {
					// 			debug = true;
					// 		}
					// 	}										
					memcpy(output1, output, sizeof(unsigned int) * n);
					long faulty_gate;
					if(j < 7){
						faulty_gate = GA.V[gate_index[i]].getInNeighbor(j);
						output1[faulty_gate] = ALL1;
					}
					else{
						faulty_gate = gate_index[i]; //output fault
						output1[faulty_gate] = ALL1;
					}

					forward_list<long> active_gates;
					unordered_set<long> visited_gates;
					if(j < 7){
                        active_gates.push_front(gate_index[i]);
					    visited_gates.insert(gate_index[i]);
                    }
					else{
                        for(long k=0; k<GA.V[gate_index[i]].getOutDegree(); k++) {
							long temp = GA.V[gate_index[i]].getOutNeighbor(k);
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
						if ( (j < 7) && (curr_gate_id == gate_index[i]) && (GA.V[faulty_gate].getOutDegree() > 1) ){
							output1[faulty_gate] = output[faulty_gate];
						}
						if (output1[curr_gate_id] != output[curr_gate_id]) {
							if (type[curr_gate_id] == PO) {
								detected_faults++;
								fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[1][j];
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



template <class vertex>
void Compute(graph<vertex> &GA, commandLine P)
{ 
    // int current_pid = GetCurrentPid(); // or you can set a outside program pid
    // float cpu_usage_ratio = GetCpuUsageRatio(current_pid);
    // float memory_usage = GetMemoryUsage(current_pid);
	srand((unsigned)time(NULL));
	char *typeFilename = P.getOptionValue("-t");
	char *muxinfoFilename = P.getOptionValue("-m");
    n = GA.n;
	long mux_gate_num = 0;
	std::ifstream in(typeFilename);
	std::ofstream out("output.txt", std::ios::binary);
	std::string line;
	std::vector<std::string> container;
	long count = 0;

    process_id = GA.getProcessId();
	processes = GA.getProcessSize();
    // local_bitsize = BITSIZE / processes;
    local_bitsize = 32;
    for (int i=0; i<processes; i++){
        setbit(finish_flag, i);
    }
    finish_flag_global = finish_flag;

	/**
	 * 读取输入文件，按行读取
	 * 每行表示一个gate，包含 gate_id、gate_type数据
	 * 
	 * 输入：
	 * 文件名 typeFilename
	 * 输出：
	 * type[i]		gate_id==i的gate的type
	 * index_PI		gate_type==PI的gate的id
	 * noPI			gate_type==PI的gate的数量
	*/

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
    if (process_id == 0) {
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
    }

	/**
	 * 层级化，给每个gate赋予一个level,使用层次遍历
	 * 
	 * 输入：
	 * GA			图数据
	 * type[i]		gate_id==i的gate的type
	 * noPI			gate_type==PI的gate数量
	 * noTIE_0		gate_type==TIE_0的gate数量
	 * 
	 * 输出：
	 * level[i]		gate_id==i的gate的level
	 * maxlevel		图中level最高的gate的level
	*/
	level = new unsigned int[n];

	bool *frontier = newA(bool, n);

	if (process_id == 0) {
	    printf("levelize start!\n");
    }

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
	if (process_id == 0) {  
	    printf("max level: %d\n", max_level);
    }

	if (process_id == 0) {
	    printf("levelize finished!\n");
    }

	/**
	 * 生成故障列表
	 * 
	 * 输入：
	 * GA				图数据
	 * add_SA_Fault		故障注入信息
	 * 
	 * 输出：
	 * fault_list[i]	gate_id==i的gate的故障信息
	 * 					0000 0001 0000 0010 
	 * 					低八位表示SA0故障，高八位表示SA1故障，
	 * 					7位和15位表示输出故障，其他表示输入故障
	 * num_of_faults	生成的故障总数
	*/
	fault_list = new unsigned short[n];
    // local_output = new unsigned int[n];
    output = new int[n];
    output1 = new int[n];
    cnt_in = new unsigned char[n];

	num_of_faults = 0;
	detected_faults = 0;

	if (process_id == 0) {
        printf("fault list generating start!\n");
    }
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

	if (process_id == 0) { 
	    printf("num of faults: %ld\n", num_of_faults);
	    printf("fault list generating finished!\n");
    }

	gate_index.resize(n,0);
	for(int i = 0; i < n; i++){
		gate_index[i] = i;
	}
	
	//int under = ceil((float)(n/processes));
	// for(int i = 0 ; i < processes ; i++){
	// 	sort(gate_index.begin() + under * i, gate_index.begin() + under * (i + 1) - 1 , [GA](int a, int b){
	// 		return GA.V[a].getInDegree() > GA.V[b].getInDegree();
	// 	});
	// }
#ifdef SORT
	sort(gate_index.begin(), gate_index.end(), [GA](int a, int b){
			return GA.V[a].getInDegree() > GA.V[b].getInDegree();
		});
#endif
    // ATPG开始
    // FAN算法预处理
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
	for (long i = 546; i < n; i++)
	{
		if (fault_list[i] != 0)
		{
			// 先处理SA0
			for (unsigned int j = 7; j < 8; j++)
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
	sychronize_timer.start();
	if(process_id == 0) {
		printf("----------\t----------\t--------\t-----------\t--------\n");
		printf("#生成向量数\t#故障覆盖率\t#剩余故障\t本轮检测故障\t运行时间\n");
	}

	//准备写文件
	// ofstream ofs;
	// ofs.open("test_pattern_"+to_string(process_id),ios::out);


	//long i = process_id * under; i < under + process_id * under ; i++
	for (long i = process_id; i < n ; i+=processes)
	{ 
		if (fault_list[gate_index[i]] != 0)
		{
			// 先处理SA0
			for (unsigned int j = 0; j < 8; j++)
			{
				if ((fault_list[gate_index[i]] & add_SA_Fault[0][j]) == 0)
					continue;
				Log("SA0 fault.");
				Log("Fault gate: %ld, Fault line: %d, gate type: %d.", i, j, type[gate_index[i]]);
				// printf("Fault gate: %ld, Fault line: %d, gate type: %d, sa0.\n", i, j, type[i]);
				//if(generated_patterns>5000)
                FAULT pCurrentFault;
                pCurrentFault.gateId=gate_index[i];
                pCurrentFault.faultLine = j;
                pCurrentFault.ftype = SA_0;
                Log("fan start");

				fan_timer.next();
				/**
				 * fan1函数
				 * 输入：
				 * GA				图数据
				 * iNoGate			图中顶点数量
				 * iMaxLevelAdd2	图中顶点的最高层级+2
				 * inoPI			PI门数量
				 * iNoPO			PO门数量
				 * pCurrentFault	故障（门id，故障序号，故障类型）
				 * iMaxBackTrack	最高回溯次数
				 * piNoBackTrack	指向整数的指针，用于记录回溯次数
				 * 
				 * 输出：
				 * rt				函数运行结果
				 * FANoutput[]		门输出（向量）
				*/
                rt = fan1(GA, n, max_level+2,  noPI, noPO, pCurrentFault, g_iMaxBackTrack1, &iNoBackTrack);
				fan_time += fan_timer.next();
				
                Log("fan end");
				// rt = fan(GA, i, j, 0);
                //return;
                if(rt == NO_TEST){
					noRedundant_faults++;
					fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[0][j];
					Log("Fan end, no pattern generated!");					
				} 
				if (rt == TEST_FOUND) {
					detected_faults++;
					fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[0][j];
					generated_patterns++;
					int bit = (generated_patterns % local_bitsize);
					
					parallel_for(long k=0; k<n; k++) {
						if (type[k] == PI) {
							int temp = rand() % 2;
							if (FANoutput[k] == 0 || (FANoutput[k] == X && temp ==0)) {
								clrbit(output[k], bit);
							} else {
								setbit(output[k], bit);
							}
						} else if (type[k] == TIE_0) {
                            clrbit(output[k], bit);
                        } else if (type[k]==AND || type[k]==NAND) {
                            setbit(output[k], bit);
                        }
					}

					// for(long k = 0 ; k < n ; k++){
					// 	 output[k] & BITMASK[bit] != ALL0 ? ofs << 1 : ofs << 0;
					// }
					// ofs << endl;

					
					//结束标志看每个机器是不是跑完了
					//generated_patterns同步每个机器生成的测试向量数量汇总到一起
					//detected_faults同步故障数量
					sychronize_timer.next();
					if (bit == 0) {
                        // parallel_for(long k=0; k<n; k++) {
                        //     output1[k] = local_output[k] << (local_bitsize * process_id);
                        // }
                        // MPI_Allreduce(output1, output, n, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
                        MPI_Allreduce(&finish_flag,&finish_flag_global,1,MPI_UNSIGNED,MPI_BAND,MPI_COMM_WORLD);
						fault_sim(GA);
                        MPI_Reduce(&generated_patterns,&generated_patterns_global,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
						MPI_Reduce(&detected_faults,&detected_faults_global,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
						if (process_id == 0) {
                            long curr_detected_faults = fault_in_list - num_of_faults + detected_faults_global;
							fault_in_list = num_of_faults - detected_faults_global;
                            double coverage=(double)detected_faults_global/(double)num_of_faults*100.00;
							printf("----------\t----------\t--------\t-----------\t--------\n");
							atpg_time += atpg_timer.next();
							printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns_global, coverage, fault_in_list, curr_detected_faults, atpg_time);
                        }
					}
					sychronize_time += sychronize_timer.next();
					Log("FAN success!%d", generated_patterns);
				}			
		 	} 
			// 再处理SA1 
			for (unsigned int j = 0; j < 8; j++)
			{
				if ((fault_list[gate_index[i]] & add_SA_Fault[1][j]) == 0)
					continue;
				Log("SA1 fault.");
				Log("Fault gate: %ld, Fault line: %d, gate type: %d.", i, j, type[gate_index[i]]);
				// printf("Fault gate: %ld, Fault line: %d, gate type: %d, sa1.\n", i, j, type[i]);
				//if(generated_patterns>5000)
                FAULT pCurrentFault;
                pCurrentFault.gateId=gate_index[i];
                pCurrentFault.faultLine = j;
                pCurrentFault.ftype = SA_1;
                Log("fan start");

                fan_timer.next();
                rt = fan1(GA, n, max_level+2,  noPI, noPO, pCurrentFault, g_iMaxBackTrack1, &iNoBackTrack);
				fan_time += fan_timer.next();

                Log("fan end");
                if(rt == NO_TEST){
					//noRedundant_faults++;
					fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[1][j];
					Log("Fan end, no pattern generated!");					
				} 
				if (rt == TEST_FOUND) {
					detected_faults++;
					fault_list[gate_index[i]] = fault_list[gate_index[i]] & del_SA_Fault[1][j];
					generated_patterns++;
					int bit = (generated_patterns % local_bitsize);
					//srand((unsigned)time(NULL));
					parallel_for(long k=0; k<n; k++) {
						if (type[k] == PI) {
							int temp = rand() % 2;
							if (FANoutput[k] == 0 || (FANoutput[k] == X && temp ==0)) {
                                clrbit(output[k], bit);
							} else {
								setbit(output[k], bit);
							}
						} else if (type[k] == TIE_0) {
                            clrbit(output[k], bit);
                        } else if (type[k]==AND || type[k]==NAND) {
                            setbit(output[k], bit);
                        }
					}
					
					// for(long k = 0 ; k < n ; k++){
					// 	 output[k] & BITMASK[bit] != ALL0 ? ofs << 1 : ofs << 0;
					// }
					// ofs << endl;
					
					sychronize_timer.next();
					if (bit == 0) {
                        // parallel_for(long k=0; k<n; k++) {
                        //     output1[k] = local_output[k] << (local_bitsize * process_id);
                        // }
                        // MPI_Allreduce(output1, output, n, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
						MPI_Allreduce(&finish_flag,&finish_flag_global,1,MPI_UNSIGNED,MPI_BAND,MPI_COMM_WORLD);
						fault_sim(GA);
                        MPI_Reduce(&generated_patterns,&generated_patterns_global,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
						MPI_Reduce(&detected_faults,&detected_faults_global,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
						if (process_id == 0) {
                            long curr_detected_faults = fault_in_list - num_of_faults + detected_faults_global;
							fault_in_list = num_of_faults - detected_faults_global;
                            double coverage=(double)detected_faults_global/(double)num_of_faults*100.00;
							printf("----------\t----------\t--------\t-----------\t--------\n");
							atpg_time += atpg_timer.next();
							printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns_global, coverage, fault_in_list, curr_detected_faults, atpg_time);
                        }
					}
					sychronize_time += sychronize_timer.next();
					Log("FAN success!%d", generated_patterns);
				}
			}
		}
	}
	sychronize_timer.next();

   if (generated_patterns % local_bitsize != 0) {
        // parallel_for(long k=0; k<n; k++) {
        //     output1[k] = local_output[k] << (local_bitsize * process_id);
        // }
        // MPI_Allreduce(output1, output, n, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
		MPI_Allreduce(&finish_flag,&finish_flag_global,1,MPI_UNSIGNED,MPI_BAND,MPI_COMM_WORLD);
        fault_sim(GA);
        MPI_Reduce(&generated_patterns,&generated_patterns_global,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&detected_faults,&detected_faults_global,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		if (process_id == 0) {
            long curr_detected_faults = fault_in_list - num_of_faults + detected_faults_global;
			fault_in_list = num_of_faults - detected_faults_global;
            double coverage=(double)detected_faults_global/(double)num_of_faults*100.00;
			printf("----------\t----------\t--------\t-----------\t--------\n");
			atpg_timer.next();
			printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns_global, coverage, fault_in_list, curr_detected_faults, atpg_time);
        }

    }
	
    clrbit(finish_flag, process_id);
    while (finish_flag_global != 0) {
        // MPI_Allreduce(output1, output, n, MPI_UNSIGNED, MPI_BOR, MPI_COMM_WORLD);
        MPI_Allreduce(&finish_flag,&finish_flag_global,1,MPI_UNSIGNED,MPI_BAND,MPI_COMM_WORLD);
		MPI_Reduce(&generated_patterns,&generated_patterns_global,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(&detected_faults,&detected_faults_global,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		//cout << "process " << process_id << ":" << finish_flag << " " << finish_flag_global << endl;
		// for(int i = 0 ; i < processes ; i++){
		// 	if(i == process_id){
		// 		printf("   process_id:								   : %d\n",process_id);
		// 		printf("   good simulation time                        : %.3lf sec\n", goodsim_time);
		// 		printf("   fault simulation time                       : %.2lf sec\n", faultsim_time);
		// 		printf("   fan time                                    : %.2lf sec\n", fan_time);
		// 		printf("   total time								   : %.2lf sec\n", fan_time+faultsim_time+goodsim_time);
		// 		printf("\n");
		// 	}
		// }
	} 
	sychronize_time += sychronize_timer.next();
	sychronize_time -= faultsim_time;
    cout << "process " << process_id << ": goodsim " << goodsim_time << " faultsim " << faultsim_time <<" fan "<< fan_time <<" sychronize " << sychronize_time<<endl;
	
	// int RankID;
	// MPI_Comm_rank(MPI_COMM_WORLD, &RankID);
	// cout << " RankID : "<< RankID << endl;
	if(process_id == 0){
		printf("   process_id                                   : %d\n",process_id);
		printf("   good simulation time                         : %.3lf sec\n", goodsim_time);
		printf("   fault simulation time                        : %.2lf sec\n", faultsim_time);
		printf("   fan time                                     : %.2lf sec\n", fan_time);
		printf("   sychronize time                               : %.2lf sec\n", sychronize_time);
		printf("   total time                                   : %.2lf sec\n", atpg_time);
		printf("   generated_patterns                           : %d\n", generated_patterns_global);
		printf("   Fault coverage                               : %.3lf %%\n", (double)detected_faults_global/(double)num_of_faults*100.00);
		printf("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// if (process_id == 0) {
    //         long curr_detected_faults = fault_in_list - num_of_faults + detected_faults_global;
	// 		fault_in_list = num_of_faults - detected_faults_global;
    //         double coverage=(double)detected_faults_global/(double)num_of_faults*100.00;
	// 		printf("----------\t----------\t--------\t-----------\t--------\n");
	// 		atpg_timer.next();
	// 		printf("%d\t\t%lf %%\t%ld\t\t%ld\t\t%.2lf sec\n", generated_patterns_global, coverage, fault_in_list, curr_detected_faults, atpg_time);
    //     }
	
#endif 
	
	// printf("1. Test pattern generation results\n");
	// printf("   Number of test patterns                     : %d\n", generated_patterns);
	// printf("   num of faults                               : %ld\n", num_of_faults);
	// printf("   Number of detected faults                   : %ld\n", detected_faults);
	// printf("   Number of identified redundant faults       : %ld\n", noRedundant_faults);
	// printf("   Number of aborted faults                    : %ld\n", num_of_faults - detected_faults - noRedundant_faults);
	// printf("   Fault coverage                              : %.3lf %%\n", (double) detected_faults / (double) num_of_faults * 100.0);
	// printf("\n");

	// printf("2. CPU time\n");

	

	// std::cout << "current pid: " << current_pid << std::endl;
    // std::cout << "cpu usage ratio: " << cpu_usage_ratio * 100 << "%" << std::endl;
    // std::cout << "memory usage: " << memory_usage << "MB" << std::endl;
	free(FANoutput);   
	free(freach);
	free(active_set);
	free(changed);
#endif
	free(ltype);
	free(dpo);
	delete []type;
	delete []level;
	delete []fault_list;


}