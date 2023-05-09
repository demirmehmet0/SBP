/*
This code is adapted from Boyar and Peralta's algorithm which can be found here: https://github.com/thomaspeyrin/XORreduce

In order to run this program, you need to compile it with C++11 or newer. The following three arguments can be added:
TIME_LIMIT - The amount of time (in seconds) allowed to run the program (default 1000)
XOR2C - Cost of implementing XOR2 gate (default 1.0)
XOR3C - Cost of implementing XOR3 gate (default 1.625)
An example command will be:
        g++ -std=c++11 -o main_globalopt.out -D XOR2C=1.0 -D XOR3C=1.625 -D TIME_LIMIT=1000 main_globalopt.cpp

To run the program, we feed the matrix into the program:
        ./main_globalopt.out < ./test/testmat.txt

An example of an input matrix will be:
15 15
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 1 0 0 0 1 0 1 0 0
0 0 0 0 0 0 0 0 0 1 1 0 1 0 0
1 0 0 0 0 0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 1 0 1
0 0 0 0 1 0 0 0 0 0 0 0 1 0 0
0 0 0 0 1 0 0 0 1 0 0 0 1 0 0
0 0 0 0 0 1 0 1 0 0 0 0 0 0 0
*/


#include <math.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

const int MaxBaseSize=1000;
const bool PRINTROWS=true;
// const float XOR2 = 1.0, XOR3 = 1.625;
        
#define SIZE 32
#define LOOP_SIZE 100
#define LARGE 100000
#define computeSum(table,table_size,s) s=0; for(int sum_i=0; sum_i<table_size; sum_i++) s+=table[sum_i];

// #ifndef OPTION
// #define OPTION -1
// #endif

#ifndef XOR4C
#ifndef XOR3C
#ifndef XOR2C
#define XOR2C 1.0
#endif
#endif
#endif

#ifndef TIME_LIMIT
#define TIME_LIMIT 1000
#endif

enum class Gate { XOR2, XOR3, XOR4 };

int NumInputs;
int NumTargets;
int XorCount, Xor2Count, Xor3Count, Xor4Count, depth;
float XorCost1, XorCost2, XorCost3, XorCost4;
long long int Target[MaxBaseSize];
int Dist[MaxBaseSize]; //distance from current base to Target[i]
int NDist[MaxBaseSize]; //what Dist would be if NewBase was added
long long int Base[MaxBaseSize];
string Program[MaxBaseSize];
int BaseSize;
int TargetsFound;
int InitDist[MaxBaseSize]; // storing the initial distance
long long int NewBase; //global variable containing a candidate new base
mt19937 rand_generator;

struct Element
{
    int parent_i;
    int parent_j;
    int parent_k;
    int parent_l;
    int newDist[SIZE];
    Gate gate;
};

void InitBase(); // refresh the base for subsequent rounds
void ReadTargetMatrix();
bool is_base(long long int x);
int NewDistance(int u); //calculates the distance from the base to Target[u]
void TotalDistance(); //returns the sum of distances to targets
void TotalDistanceXOR3(Gate gate);
int NewDistanceXOR3(int u, Gate gate);
bool reachable(long long int T, int K, int S);
bool EasyMove(); //if any two bases add up to a target, pick them
int EasyMoveXOR3();
void PickNewBaseElement();
void PickNewBaseElementXOR3();
void refreshDist(); // refresh the distance to the targets for subsequent rounds
void refreshDistAndTarget();
int RNBP(Element AllElements[], int counter); // RNBP
int A1(Element AllElements[], int counter); // A1
int A2(Element AllElements[], int counter); // A2
int calculateDist(int A[],int length);
int calculateNorm(int A[],int length);
bool filtering(int tempDist[], vector<int> filter_indices);
vector<string> direct;
unordered_map<string, int> depth_map;
int depths[1000];
unordered_map<string, vector<string> > iwsec_graph;
vector<string> insert_order;

vector<vector<int> > TargetMatrix;
vector<int> row, col;

float XOR2_ASIC1 = 2.0;
float XOR3_ASIC1 = 3.25;
float XOR4_ASIC1 = 5.0;

float XOR2_ASIC2 = 1.981;
float XOR3_ASIC2 = 3.715;
float XOR4_ASIC2 = 5.5;

float XOR2_ASIC3 = 2.5;
float XOR3_ASIC3 = 4.2;
float XOR4_ASIC3 = 6.25;

float XOR2_ASIC4 = 3.33;
float XOR3_ASIC4 = 4.66;
float XOR4_ASIC4 = 5.99;

template<class T>
void ShuffleFromIndices(vector<T> &vec, const vector<int> &indices) {
    vector<T> temp = vec;
    for(int i=0;i<temp.size();i++) {
        temp[i] = vec[indices[i]];
    }
    vec = temp;
}

int iwsec_compute_depth() {

    XorCount = 0;
    Xor2Count = 0;
    Xor3Count = 0;
    Xor4Count = 0;
    XorCost1 = 0;
    XorCost2 = 0;
    XorCost3 = 0;
    XorCost4 = 0;
    
    unordered_map<string, int> depth;

    for(auto &i:insert_order) {
        if (iwsec_graph[i].size() != 0) {
            
            if (iwsec_graph[i].size() == 2) {
                depth[i] = max(depth[iwsec_graph[i][0]], depth[iwsec_graph[i][1]]) + 1;
            } else {
                depth[i] = max({depth[iwsec_graph[i][0]], depth[iwsec_graph[i][1]], depth[iwsec_graph[i][2]]}) + 1;
            }

            // Count XOR2 and XOR3 cost
            if (iwsec_graph[i].size() == 2) {
                XorCost1 += XOR2_ASIC1;
                XorCost2 += XOR2_ASIC2;
                XorCost3 += XOR2_ASIC3;
                XorCost4 += XOR2_ASIC4;
                Xor2Count += 1;
            } else {
                XorCost1 += XOR3_ASIC1;
                XorCost2 += XOR3_ASIC2;
                XorCost3 += XOR3_ASIC3;
                XorCost4 += XOR3_ASIC4;
                Xor3Count += 1;
            }
            XorCount += 1;
        }
    }

    return max_element(depth.begin(), depth.end(),
    [](const pair<string, int>& p1, const pair<string, int>& p2) {
        return p1.second < p2.second; })->second;
}

bool iwsec_xor3_postprocessing() {

    unordered_map<string, int> fan_out;
    vector<string> variables;

    bool changed = false;

    // Compute fan out of temporary variables
    for (auto &i:iwsec_graph) {
        for (auto &j:i.second) {
            if (j[0] == 't') {
                fan_out[j] += 1;
                variables.push_back(j);
            }
        }
    }

    unordered_map<string, int> local_depth;

    // Random shuffle to arbitrarily pick in case of multiple candidates.
    // In case of equal depth, this ensures random tie break.
    random_shuffle(variables.begin(), variables.end());

    for (auto &var:variables) {
        auto i = fan_out[var];
        if (i == 1 && iwsec_graph[var].size() == 2) {
            for (auto &j:iwsec_graph) {
                // Ensure that it is still XOR2
                if (j.second.size() == 2) {
                    int idx = -1;
                    if (j.second[0] == var) {
                        idx = 0;
                    } else if (j.second[1] == var) {
                        idx = 1;
                    }
                    
                    // Found the temporary variable.
                    if (idx != -1) {
                        // Check if other index is also temporary variable with fan out 1.
                        if (j.second[(idx + 1) % 2][0] == 't') {
                            // If current variable depth is greater than the other variable depth, then continue without reducing. 
                            if (fan_out[j.second[(idx + 1) % 2]] == 1 && local_depth[j.second[(idx + 1) % 2]] < local_depth[j.second[idx]]) {
                                continue;
                            }
                        }

                        // Delete the temporary variable and add its inputs
                        j.second.erase(j.second.begin() + idx);
                        j.second.push_back(iwsec_graph[var][0]);
                        j.second.push_back(iwsec_graph[var][1]);
                        // Delete the temporary variable
                        iwsec_graph[var] = {};
                        changed = true;
                    }
                }
            }
        }
        if (iwsec_graph[var].size() == 2) {
            local_depth[var] = max(local_depth[iwsec_graph[var][0]], local_depth[iwsec_graph[var][1]]) + 1;
        } else if (iwsec_graph[var].size() == 3) {
            local_depth[var] = max({local_depth[iwsec_graph[var][0]], local_depth[iwsec_graph[var][1]], local_depth[iwsec_graph[var][2]]}) + 1;
        }
    }

    return changed;
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        cout << "Pass matrix name\n";
        return 1;
    }
    string mat = argv[1];
    // reading the target matrix in text file
    // setting up the distance
    ReadTargetMatrix();
    // Large value for initialization
    int BestCount = LARGE;
    float BestCost1 = LARGE;
    float BestCost2 = LARGE;
    float BestCost3 = LARGE;
    float BestCost4 = LARGE;
    int BestDepth = LARGE;
    clock_t start = clock();
    int best_time;
    int iterations = TIME_LIMIT;
    ofstream logs;
    string log_file;
    #ifdef LESS_THAN
    	#ifdef XOR4C
	    logs.open(mat + "_bp_xor4_" + to_string(XOR4C) + "_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + "_less_than.log", std::ios_base::app);
	    log_file = mat + "_bp_xor4_" + to_string(XOR4C) + "_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + "_less_than.log";
    	#elif defined(XOR3C)
	    logs.open(mat + "_bp_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + "_less_than.log", std::ios_base::app);
	    log_file = mat + "_bp_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + "_less_than.log";
        #elif defined(IWSEC)
        logs.open(mat + "_bp_xor3_iwsec_" + to_string(XOR2C) + "_less_than.log", std::ios_base::app);
	    log_file = mat + "_bp_xor3_iwsec_" + to_string(XOR2C) + "_less_than.log";
	#else
	    logs.open(mat + "_bp_xor2_" + to_string(XOR2C) + "_less_than.log", std::ios_base::app);
	    log_file = mat + "_bp_xor2_" + to_string(XOR2C) + "_less_than.log";
	#endif
    #else
    	#ifdef XOR4C
	    logs.open(mat + "_bp_xor4_" + to_string(XOR4C) + "_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + ".log", std::ios_base::app);
	    log_file = mat + "_bp_xor4_" + to_string(XOR4C) + "_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + ".log";
    	#elif defined(XOR3C)
	    logs.open(mat + "_bp_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + ".log", std::ios_base::app);
	    log_file = mat + "_bp_xor3_" + to_string(XOR3C) + "_xor2_" + to_string(XOR2C) + ".log";
        #elif defined(IWSEC)
        logs.open(mat + "_bp_xor3_iwsec_" + to_string(XOR2C) + ".log", std::ios_base::app);
	    log_file = mat + "_bp_xor3_iwsec_" + to_string(XOR2C) + ".log";
	#else
	    logs.open(mat + "_bp_xor2_" + to_string(XOR2C) + ".log", std::ios_base::app);
	    log_file = mat + "_bp_xor2_" + to_string(XOR2C) + ".log";
	#endif
    #endif
    // cout << "Log file: " << log_file << "\n\n";
    // Number of rounds
    // while ((clock()-start)/CLOCKS_PER_SEC < TIME_LIMIT) {
    string log_str;
    while (iterations--) {
        BestCount = LARGE;
        BestCost1 = LARGE;
        BestCost2 = LARGE;
        BestCost3 = LARGE;
        BestCost4 = LARGE;
        BestDepth = LARGE;
        XorCount = 0;
        Xor2Count = 0;
        Xor3Count = 0;
        Xor4Count = 0;
        XorCost1 = 0;
        XorCost2 = 0;
        XorCost3 = 0;
        XorCost4 = 0;
        // refreshing the distance and base for subsequent rounds
        refreshDistAndTarget();
        InitBase();
        // refreshDist();

        // main loop
        int _returnVal2 = 0;
        while (TargetsFound < NumTargets) {
            //cout << "Targets Found: " << TargetsFound << " out of " << NumTargets << endl;
            #if defined(XOR3C) || defined(XOR4C)
            int _returnVal = EasyMoveXOR3();
            if (_returnVal == 0){ 
                PickNewBaseElementXOR3();
            }else if (_returnVal == 2){
                _returnVal2 = 2;
                cout << "false: "<< iterations << endl;
                //logs << "Result Invalid"<< endl << endl << endl << "--------------" << endl << endl << endl;
                break;
            }

            #else
            if (!EasyMove()) PickNewBaseElement();
            #endif
        }
        if (((BestCost1-XorCost1) > 1e-3) || ((BestCost2-XorCost2) > 1e-3) || ((BestCost3-XorCost3) > 1e-3) || ((BestCost4-XorCost4) > 1e-3) || ((BestCount-XorCount) > 1e-3))
        {
	    int trial_no = TIME_LIMIT - iterations;
	    cout << "Trial No: "<< trial_no << endl;
        if(_returnVal2 != 2 )
	    logs << "Trial No: "<< trial_no << endl;

        depth = max_element(depth_map.begin(), depth_map.end(),
    [](const pair<string, int>& p1, const pair<string, int>& p2) {
        return p1.second < p2.second; })->second;

            if (TargetsFound == NumTargets)
            {
                
                #ifdef IWSEC
                    while (true) {
                        if (!iwsec_xor3_postprocessing()) {
                            break;
                        }
                    }
                    depth = iwsec_compute_depth();
                #endif
                
                // cout << "SLP Heuristic XorCount: " << XorCount << endl;
                if(_returnVal2 != 2 )
                logs << "SLP Heuristic XorCount: " << XorCount;
                if(_returnVal2 != 2 )
                logs << ", Depth: " << depth << endl;
		#ifdef XOR4C
		    if(_returnVal2 != 2 )
            logs << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << ", XOR4: " << Xor4Count << endl;
		#elif defined(XOR3C)
		    if(_returnVal2 != 2 )
            logs << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << endl;
        #elif defined(IWSEC)
		    if(_returnVal2 != 2 )
            logs << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << endl;
		#else
		    if(_returnVal2 != 2 )
            logs << ", XOR2: " << Xor2Count << endl;
		#endif
                // cout << "SLP Heuristic XorCost ASIC1: " << XorCost1 << endl;
                if(_returnVal2 != 2 )
                logs << "SLP Heuristic XorCost ASIC1: " << XorCost1 << endl;
                // cout << "SLP Heuristic XorCost ASIC2: " << XorCost2 << endl;
                if(_returnVal2 != 2 )
                logs << "SLP Heuristic XorCost ASIC2: " << XorCost2 << endl;
                // cout << "SLP Heuristic XorCost ASIC3: " << XorCost3 << endl;
                if(_returnVal2 != 2 )
                logs << "SLP Heuristic XorCost ASIC3: " << XorCost3 << endl;
                // cout << "SLP Heuristic XorCost ASIC4: " << XorCost4 << endl;
                if(_returnVal2 != 2 )
                logs << "SLP Heuristic XorCost ASIC4: " << XorCost4 << endl;
            }
            
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);
            std::ostringstream oss;
            oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
            string time_stamp = oss.str();
            cout << "Timestamp: " << time_stamp << endl;
            if(_returnVal2 != 2 )
            logs << "Timestamp: " << time_stamp << endl;
            // float acc_time = (clock()-start)/CLOCKS_PER_SEC;
            // cout << "Accumulated Time: " << acc_time << endl;
            // logs << "Accumulated Time: " << acc_time << endl;
            best_time = (clock()-start)/CLOCKS_PER_SEC;
            if (BestDepth > depth) {
                BestDepth = depth;
                cout << "Depth improved: " << BestDepth << endl;
            }
            if ((BestCount-XorCount) > 1e-3) {
                BestCount = XorCount;
            cout << "XORCount improved: " << BestCount;
		#ifdef XOR4C
		    cout << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << ", XOR4: " << Xor4Count << endl;
		#elif defined(XOR3C)
		    cout << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << endl;
        #elif defined(IWSEC)
		    cout << ", XOR2: " << Xor2Count << ", XOR3: " << Xor3Count << endl;
		#else
		    cout << ", XOR2: " << Xor2Count << endl;
		#endif
                cout << "Depth: " << depth << endl;
            }
            if ((BestCost1-XorCost1) > 1e-3) {
                BestCost1 = XorCost1;
                cout << "ASIC1 improved: " << BestCost1 << endl;
            }
            if ((BestCost2-XorCost2) > 1e-3) {
                BestCost2 = XorCost2;
                cout << "ASIC2 improved: " << BestCost2 << endl;
            }
            if ((BestCost3-XorCost3) > 1e-3) {
                BestCost3 = XorCost3;
                cout << "ASIC3 improved: " << BestCost3 << endl;
            }
            if ((BestCost4-XorCost4) > 1e-3) {
                BestCost4 = XorCost4;
                cout << "ASIC4 improved: " << BestCost4 << endl;
            }
            cout << "Logfile: " << log_file << endl;
            
            #ifdef IWSEC
                for(auto &i:insert_order) {
                    if (iwsec_graph[i].size() != 0) {
                    
                        string out  = i;
                        out += " = ";
                        for(auto &j:iwsec_graph[i]) {
                            out += j + " + ";
                        }
                        out = out.substr(0, out.size() - 3);
                        if(_returnVal2 != 2 )
                        logs << out << endl;
                    }
                }
                iwsec_graph.clear();
                insert_order.clear();
            #else
                for (int j = 0; j < XorCount; j++) {
                    if(_returnVal2 != 2 )
                    logs << Program[NumInputs + j] << " (" << depths[j + NumInputs] << ") " << endl;
                }
            #endif

            if (direct.size()) {
                for (auto i:direct) {
                    if(_returnVal2 != 2 )
                    logs << i << endl;
                }
            }

            // cout << "Row:\n";
            if(_returnVal2 != 2 )
            logs << "Row: " << endl;
            log_str = "";
            for(auto i:row) {
                // cout << i << " ";
                log_str += to_string(i) + " ";
            }
            if(_returnVal2 != 2 )
            logs << log_str;
            
            log_str = "";
            // cout << "\nCol:\n";
            if(_returnVal2 != 2 )
            logs << "\nCol: " << endl;
            for(auto i:col) {
                // cout << i << " ";
                log_str += to_string(i) + " ";
            }
            log_str += "\n";
            if(_returnVal2 != 2 )
            logs << log_str << endl;
            cout << "\n\n";
        }
    }
    logs.close();
    // check if it is stablized
    // if (best_time < ((clock()-start)/2)/CLOCKS_PER_SEC) cout << "Stable" << endl;
    // else cout << "Unstable" << endl;

    return 0;
}

void InitBase() {
    TargetsFound = 0;
    Base[0] = 1;
    Program[0] = "x00";
    XorCost1 = 0;
    XorCost2 = 0;
    XorCost3 = 0;
    XorCost4 = 0;
    stringstream ss;
    string s;
    for (int i = 1; i < NumInputs; i++) {
        ss << i;
        s = ss.str();
	if (s.length() < 2)
	    s.insert(s.begin(), 2 - s.length(), '0');
        Base[i] = 2*Base[i-1];
        Program[i] = "x" + s;
        ss.str("");
    }
    BaseSize = NumInputs; //initial base is just the xi's
    for (int i = 0; i < NumTargets; i++) {
        if (Dist[i] == 0) {
            TargetsFound++;
        }
    }
}

void TotalDistance() { //returns the sum of distances to targets
    int D = 0;
    int t;
    for (int i = 0; i < NumTargets; i++) {
        t = NewDistance(i);
        NDist[i] = t;
    }
}

void TotalDistanceXOR3(Gate gate) { //returns the sum of distances to targets
    int D = 0;
    int t;
    for (int i = 0; i < NumTargets; i++) {
        t = NewDistanceXOR3(i, gate);
        NDist[i] = t;
    }
}

bool EasyMove() {
    int t;
    bool foundone = false;

    //see if anything in the distance vector is 1
    for(int i = 0; i < NumTargets; i++) {
        if (Dist[i] == 1) {
            foundone = true;
            t = i;
            break;
        }
    }
    if (!foundone) {
        return false;
    }
    //update Dist array
    NewBase = Target[t];
    for (int u = 0; u < NumTargets; u++) {
        Dist[u] = NewDistance(u);
    }
    //update Base with NewBase
    Base[BaseSize] = NewBase;
    //find which lines in Base caused this
    string a,b;
    for (int i = 0; i < BaseSize; i++) {
        for (int j = i + 1; j < BaseSize; j++) {
            if ((Base[i] ^ Base[j]) == Target[t]) {
                a = Program[i].substr(0, Program[i].find(" "));
                b = Program[j].substr(0, Program[j].find(" "));
                break;
            }
        }
    }
    stringstream ss;
    string s1;
    ss << t;
    s1 = ss.str();
    if (s1.length() < 2)
        s1.insert(s1.begin(), 2 - s1.length(), '0');
    ss.str("");
    Program[BaseSize] = "y" + s1 + " = " + a + " + " + b;
    depth_map["y" + s1] = max(depth_map[a], depth_map[b]) + 1;
    depths[BaseSize] = max(depth_map[a], depth_map[b]) + 1;
    iwsec_graph["y" + s1] = {a, b};
    insert_order.push_back("y" + s1);
    // XorCost += XOR2C;
    XorCost1 += XOR2_ASIC1;
    XorCost2 += XOR2_ASIC2;
    XorCost3 += XOR2_ASIC3;
    XorCost4 += XOR2_ASIC4;
    BaseSize++;
    XorCount++;
    Xor2Count++;
    TargetsFound++;
    return true;
}

int EasyMoveXOR3() {
    int t;
    bool foundone = false;
    Gate gateFound;

    //see if anything in the distance vector is 1
    for(int i = 0; i < NumTargets; i++) {
	#ifdef XOR4C
        if (Dist[i] == 3 || Dist[i] == 2 || Dist[i] == 1) {
            gateFound = Dist[i] == 3 ? Gate::XOR4 : (Dist[i] == 2 ? Gate::XOR3 : Gate::XOR2);
            foundone = true;
            t = i;
            break;
        }
	#elif defined(XOR3C)
        if (Dist[i] == 2 || Dist[i] == 1) {
            gateFound = Dist[i] == 2 ? Gate::XOR3 : Gate::XOR2;
            foundone = true;
            t = i;
            break;
        }
	#endif
    }
    if (!foundone) {
        return 0;
    }
    //update Dist array
    NewBase = Target[t];
    for (int u = 0; u < NumTargets; u++) {
        Dist[u] = NewDistanceXOR3(u, gateFound);
    }
    if (Dist[t] != 0) {
        cout << "SHOULD NOT BE CALLED";
        exit(0);
    }
    //update Base with NewBase
    Base[BaseSize] = NewBase;
    //find which lines in Base caused this
    string a,b,c,d;
    bool _foundone = true;
    //cout << "Gate found: " << endl;
    if (gateFound == Gate::XOR2) {
        for (int i = 0; i < BaseSize; i++) {
            for (int j = i + 1; j < BaseSize; j++) {
                if ((Base[i] ^ Base[j]) == Target[t]) {
                    a = Program[i].substr(0, Program[i].find(" "));
                    b = Program[j].substr(0, Program[j].find(" "));
                    if(depths[i] + 1 > 4 || depths[j] + 1 > 4){
                        _foundone = false;
                    }
                    break;
                }
            }
        }
    } else if (gateFound == Gate::XOR3) {
        bool found = false;
        for(int i = 0; i < BaseSize; i++) {
            for(int j = i + 1; j < BaseSize; j++) {
                for(int k = j + 1; k < BaseSize; k++) {
                    if ((Base[i] ^ Base[j] ^ Base[k]) == Target[t]) {
                        a = Program[i].substr(0, Program[i].find(" "));
                        b = Program[j].substr(0, Program[j].find(" "));
                        c = Program[k].substr(0, Program[k].find(" "));
                        found = true;
                        if(depths[i] + 1 > 4 || depths[j] + 1 > 4 || depths[k] + 1 > 4){
                            _foundone = false;
                        }
                        break; 
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
    } else {
        bool found = false;
        for(int i = 0; i < BaseSize; i++) {
	    for(int j = i + 1; j < BaseSize; j++) {
		for(int k = j + 1; k < BaseSize; k++) {
		    for(int l = k + 1; l < BaseSize; l++) {
			if ((Base[i] ^ Base[j] ^ Base[k] ^ Base[l]) == Target[t]) {
			    a = Program[i].substr(0, Program[i].find(" "));
			    b = Program[j].substr(0, Program[j].find(" "));
			    c = Program[k].substr(0, Program[k].find(" "));
			    d = Program[l].substr(0, Program[l].find(" "));
			    found = true;
                if(depths[i] + 1 > 4 || depths[j] + 1 > 4 || depths[k] + 1 > 4 || depths[l] + 1 > 4){
                    _foundone = false;
                }
			    break;
			}
		    }
		    if (found) break;
		}
		if (found) break;
	    }
	    if (found) break;
	}
    }
    //cout << "SHOULD BE CALLED"<<endl;
    if(_foundone == false){
        /* for(int i = 0; i < BaseSize; i++){
            cout << Program[i] << endl;
        } *//* 
        cout << "SHOULD NOT BE CALLED"<<endl; */
        return 2;
    }else{
        //cout << "SHOULD BE CALLED"<<endl;
    }
    stringstream ss;
    string s1;
    ss << t;
    s1 = ss.str();
    if (s1.length() < 2)
        s1.insert(s1.begin(), 2 - s1.length(), '0');
    ss.str("");
    
    if (gateFound == Gate::XOR2) {
        Program[BaseSize] = "y" + s1 + " = " + a + " + " + b;
        // XorCost += XOR2C;
        XorCost1 += XOR2_ASIC1;
        XorCost2 += XOR2_ASIC2;
        XorCost3 += XOR2_ASIC3;
        XorCost4 += XOR2_ASIC4;
        depth_map["y" + s1] = max(depth_map[a], depth_map[b]) + 1;
        depths[BaseSize] = max(depth_map[a], depth_map[b]) + 1;
    } else if (gateFound == Gate::XOR3) {
        Program[BaseSize] = "y" + s1 + " = " + a + " + " + b + " + " + c;
        // XorCost += XOR3C;
        XorCost1 += XOR3_ASIC1;
        XorCost2 += XOR3_ASIC2;
        XorCost3 += XOR3_ASIC3;
        XorCost4 += XOR3_ASIC4;
        depth_map["y" + s1] = max({depth_map[a], depth_map[b], depth_map[c]}) + 1;
        depths[BaseSize] = max({depth_map[a], depth_map[b], depth_map[c]}) + 1;
    } else {
        Program[BaseSize] = "y" + s1 + " = " + a + " + " + b + " + " + c + " + " + d;
        // XorCost += XOR4C;
        XorCost1 += XOR4_ASIC1;
        XorCost2 += XOR4_ASIC2;
        XorCost3 += XOR4_ASIC3;
        XorCost4 += XOR4_ASIC4;
        depth_map["y" + s1] = max({depth_map[a], depth_map[b], depth_map[c], depth_map[d]}) + 1;
        depths[BaseSize] = max({depth_map[a], depth_map[b], depth_map[c], depth_map[d]}) + 1;
    }
    BaseSize++;
    XorCount++;
    if (gateFound == Gate::XOR2)
	Xor2Count++;
    else if (gateFound == Gate::XOR3)
	Xor3Count++;
    else
	Xor4Count++;
    TargetsFound++;
    return 1;
}

// PickNewBaseElement is only called when there are no 1's in Dist[]
void PickNewBaseElementXOR3() {
    // Allocate memory for all possible bases
    Element* AllElements = new Element[BaseSize*(BaseSize-1)*(BaseSize*BaseSize-4*BaseSize+5)];
    int counter = 0; // counter to track last element

    for (int i = 0; i < BaseSize; i++) {
        for (int j = i+1; j < BaseSize; j++) {
            if(depths[i] + 1 > 4 || depths[j] + 1 > 4){
                //cout<<"Depth limit exceeded"<<endl;
                continue;
            }
            NewBase = Base[i] ^ Base[j];
            TotalDistanceXOR3(Gate::XOR2); //this calculates NDist[]

            // Putting in the data into the array
            for (int k = 0; k < NumTargets; k++)
            {
                AllElements[counter].newDist[k] = NDist[k];
            }
            AllElements[counter].parent_i = i;
            AllElements[counter].parent_j = j;
            AllElements[counter].gate = Gate::XOR2;
            counter++;

            for(int k = j+1; k < BaseSize; k++) {
                if(depths[i] + 1 > 4 || depths[j] + 1 > 4 || depths[k] + 1 > 4){
                    /* cout<<"Depth limit exceeded3"<<endl; */
                    continue;
                }
                NewBase = Base[i] ^ Base[j] ^ Base[k];
                TotalDistanceXOR3(Gate::XOR3); //this calculates NDist[]

                // Putting in the data into the array
                for (int k = 0; k < NumTargets; k++)
                {
                    AllElements[counter].newDist[k] = NDist[k];
                }
                AllElements[counter].parent_i = i;
                AllElements[counter].parent_j = j;
                AllElements[counter].parent_k = k;
                AllElements[counter].gate = Gate::XOR3;
                counter++;

		#ifdef XOR4C
		for(int l = k+1; l < BaseSize; l++) {
                    NewBase = Base[i] ^ Base[j] ^ Base[k] ^ Base[l];
                    TotalDistanceXOR3(Gate::XOR4); //this calculates NDist[]
		        
                    // Putting in the data into the array
                    for (int k = 0; k < NumTargets; k++)
                    {
                        AllElements[counter].newDist[k] = NDist[k];
                    }
                    AllElements[counter].parent_i = i;
                    AllElements[counter].parent_j = j;
                    AllElements[counter].parent_k = k;
		    AllElements[counter].parent_l = l;
                    AllElements[counter].gate = Gate::XOR4;
                    counter++;
		}
		#endif
            }
        }
    }

    // for (int i = 0; i < BaseSize - 2; i++) {
    //     for (int j = i+1; j < BaseSize - 1; j++) {
    //         for (int k = j+1; k < BaseSize; k++) {
    //             NewBase = Base[i] ^ Base[j] ^ Base[k];
    //             TotalDistanceXOR3(Gate::XOR3); //this calculates NDist[]

    //             // Putting in the data into the array
    //             for (int l = 0; l < NumTargets; l++)
    //             {
    //                 AllElements[counter].newDist[l] = NDist[l];
    //             }
    //             AllElements[counter].parent_i = i;
    //             AllElements[counter].parent_j = j;
    //             AllElements[counter].parent_k = k;
    //             AllElements[counter].gate = Gate::XOR3;
    //             counter++;
    //         }
    //     }
    // }

    int chosen = RNBP(AllElements,counter);
    // selecting the best pair to XOR
    // if (OPTION == 1) chosen = RNBP(AllElements,counter);
    // else if (OPTION == 2) chosen = A1(AllElements,counter);
    // else if (OPTION == 3) chosen = A2(AllElements,counter);
    
    Gate chosenGate = AllElements[chosen].gate;

    // Update using the result returned by the criteria
    int bestparent_i, bestparent_j, bestparent_k, bestparent_l;
    if (chosenGate == Gate::XOR2) {
        bestparent_i = AllElements[chosen].parent_i;
        bestparent_j = AllElements[chosen].parent_j;
        Base[BaseSize] = Base[bestparent_i] ^ Base[bestparent_j];
    } else if (chosenGate == Gate::XOR3) {
        bestparent_i = AllElements[chosen].parent_i;
        bestparent_j = AllElements[chosen].parent_j;
        bestparent_k = AllElements[chosen].parent_k;
        Base[BaseSize] = Base[bestparent_i] ^ Base[bestparent_j] ^ Base[bestparent_k];
    } else {
        bestparent_i = AllElements[chosen].parent_i;
        bestparent_j = AllElements[chosen].parent_j;
        bestparent_k = AllElements[chosen].parent_k;
	bestparent_l = AllElements[chosen].parent_l;
        Base[BaseSize] = Base[bestparent_i] ^ Base[bestparent_j] ^ Base[bestparent_k] ^ Base[bestparent_l];
    }

    // Update the Dist Array
    for (int i = 0; i < NumTargets; i++) {
        Dist[i] = AllElements[chosen].newDist[i];
    }


    string a = Program[bestparent_i].substr(0, Program[bestparent_i].find(" "));
    string b = Program[bestparent_j].substr(0, Program[bestparent_j].find(" "));
    string c, d;
    if (chosenGate != Gate::XOR2) {
        c = Program[bestparent_k].substr(0, Program[bestparent_k].find(" "));
	if (chosenGate != Gate::XOR3) {
	    d = Program[bestparent_l].substr(0, Program[bestparent_l].find(" "));
	}
    }
    stringstream ss;
    string s2;
    ss << XorCount;
    s2 = ss.str();
    if (s2.length() < 2)
        s2.insert(s2.begin(), 2 - s2.length(), '0');
    ss.str("");
    if (chosenGate == Gate::XOR2) {
        Program[BaseSize] = "t" + s2 + " = " + a + " + " + b;
        depth_map["t" + s2] = max(depth_map[a], depth_map[b]) + 1;
        depths[BaseSize] = max(depth_map[a], depth_map[b]) + 1;
        // XorCost += XOR2C;
        XorCost1 += XOR2_ASIC1;
        XorCost2 += XOR2_ASIC2;
        XorCost3 += XOR2_ASIC3;
        XorCost4 += XOR2_ASIC4;
    } else if (chosenGate == Gate::XOR3) {
        Program[BaseSize] = "t" + s2 + " = " + a + " + " + b + " + " + c;
        depth_map["t" + s2] = max({depth_map[a], depth_map[b], depth_map[c]}) + 1;
        depths[BaseSize] = max({depth_map[a], depth_map[b], depth_map[c]}) + 1;
        // XorCost += XOR3C;
        XorCost1 += XOR3_ASIC1;
        XorCost2 += XOR3_ASIC2;
        XorCost3 += XOR3_ASIC3;
        XorCost4 += XOR3_ASIC4;
    } else {
        Program[BaseSize] = "t" + s2 + " = " + a + " + " + b + " + " + c + " + " + d;
        depth_map["t" + s2] = max({depth_map[a], depth_map[b], depth_map[c], depth_map[d]}) + 1;
        depths[BaseSize] = max({depth_map[a], depth_map[b], depth_map[c], depth_map[d]}) + 1;
        // XorCost += XOR3C;
        XorCost1 += XOR4_ASIC1;
        XorCost2 += XOR4_ASIC2;
        XorCost3 += XOR4_ASIC3;
        XorCost4 += XOR4_ASIC4;
    }
    BaseSize++;
    XorCount++;
    if (chosenGate == Gate::XOR2)
	Xor2Count++;
    else if (chosenGate == Gate::XOR3)
	Xor3Count++;
    else
	Xor4Count++;

    // free up the memory
    free(AllElements);

    return;
}

// PickNewBaseElement is only called when there are no 1's in Dist[]
void PickNewBaseElement() {
    cout<<"pick" <<endl;
    // Allocate memory for all possible bases
    Element* AllElements = new Element[BaseSize*(BaseSize-1)];
    int counter = 0; // counter to track last element

    for (int i = 0; i < BaseSize - 1; i++) {
        for (int j = i+1; j < BaseSize; j++) {
            NewBase = Base[i] ^ Base[j];
            TotalDistance(); //this calculates NDist[]

            // Putting in the data into the array
            for (int k = 0; k < NumTargets; k++)
            {
                AllElements[counter].newDist[k] = NDist[k];
            }
            AllElements[counter].parent_i = i;
            AllElements[counter].parent_j = j;
	    AllElements[counter].gate = Gate::XOR2;
            counter++;
        }
    }
    int chosen = RNBP(AllElements,counter);
    // selecting the best pair to XOR
    // if (OPTION == 1) chosen = RNBP(AllElements,counter);
    // else if (OPTION == 2) chosen = A1(AllElements,counter);
    // else if (OPTION == 3) chosen = A2(AllElements,counter);
    

    // Update using the result returned by the criteria
    int bestparent_i = AllElements[chosen].parent_i;
    int bestparent_j = AllElements[chosen].parent_j;
    Base[BaseSize] = Base[bestparent_i] ^ Base[bestparent_j];

    // Update the Dist Array
    for (int i = 0; i < NumTargets; i++) {
        Dist[i] = AllElements[chosen].newDist[i];
    }


    string a = Program[bestparent_i].substr(0, Program[bestparent_i].find(" "));
    string b = Program[bestparent_j].substr(0, Program[bestparent_j].find(" "));
    stringstream ss;
    string s2;
    ss << XorCount;
    s2 = ss.str();
    if (s2.length() < 2)
        s2.insert(s2.begin(), 2 - s2.length(), '0');
    ss.str("");
    Program[BaseSize] = "t" + s2 + " = " + a + " + " + b;
    depth_map["t" + s2] = max(depth_map[a], depth_map[b]) + 1;
    depths[BaseSize] = max(depth_map[a], depth_map[b]) + 1;
    iwsec_graph["t" + s2] = {a, b};
    insert_order.push_back("t" + s2);
    // XorCost += XOR2C;
    XorCost1 += XOR2_ASIC1;
    XorCost2 += XOR2_ASIC2;
    XorCost3 += XOR2_ASIC3;
    XorCost4 += XOR2_ASIC4;
    BaseSize++;
    XorCount++;
    Xor2Count++;

    // free up the memory
    free(AllElements);

    return;
}
// Original BP random criteria
int RNBP(Element AllElements[], int counter)
{
    // initialization
    float bestDist = LARGE;
    int bestNorm = -1*LARGE;
    float currentDist, currentNorm;
    vector<int> candidates;
    for (int i = 0; i < counter; i++)
    {   
        // float gateCost = AllElements[i].gate == Gate::XOR2 ? XOR2C : XOR3C;
	#ifdef XOR4C
            float gateCost = AllElements[i].gate == Gate::XOR2 ? 1.0 : (AllElements[i].gate == Gate::XOR3 ? XOR3C / XOR2C : XOR4C / XOR2C);
            currentDist = calculateDist(AllElements[i].newDist,NumTargets) + gateCost;
	#elif defined(XOR3C)
	    float gateCost = AllElements[i].gate == Gate::XOR2 ? 1.0 : XOR3C / XOR2C;
	    currentDist = calculateDist(AllElements[i].newDist,NumTargets) + gateCost;
	#else
	    currentDist = calculateDist(AllElements[i].newDist,NumTargets);
	#endif
        // currentDist = calculateDist(AllElements[i].newDist,NumTargets);
        currentNorm = calculateNorm(AllElements[i].newDist,NumTargets);
        #ifdef LESS_THAN
        if ((currentDist < bestDist) || (currentDist == bestDist && currentNorm < bestNorm))
        #else
        if ((currentDist < bestDist) || (currentDist == bestDist && currentNorm > bestNorm))
        #endif
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            bestNorm = currentNorm;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist && currentNorm == bestNorm)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }

    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0,candidates.size()-1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}

/* int A1(Element AllElements[], int counter)
{
    // Applying the Filter
    int nearest = 1; // change this to relax the filter of nearest targets
    int filter_dist; // keep track of the largest distance that will pass through the filter
    vector<int> filter_indices; // keep track of the indices of Target that satisfy the filter

    // sort the distance array
    int sorted_dist[NumTargets-TargetsFound];
    int next_index = 0;
    for (int i = 0; i < NumTargets; i++)
    {
        if (Dist[i] == 0) continue;
        sorted_dist[next_index++] = Dist[i];
    }
    sort(sorted_dist,sorted_dist+NumTargets-TargetsFound);

    // indices that can pass through the filter
    filter_dist = sorted_dist[min(nearest-1,NumTargets-TargetsFound-1)]; // largest distance that will pass through the filter
    for (int i = 0; i < NumTargets; i++)
    {

        if (Dist[i] <= filter_dist && Dist[i] > 0) 
        {
            filter_indices.push_back(i);
        }
    }

    // initialization
    int bestDist = LARGE;
    int bestNorm = -1*LARGE;
    int currentDist, currentNorm;
    vector<int> candidates;


    for (int i = 0; i < counter; i++)
    {
        // Filtering
        if (!filtering(AllElements[i].newDist, filter_indices)) continue;
        // Normal BP rand
        currentDist = calculateDist(AllElements[i].newDist,NumTargets);
        currentNorm = calculateNorm(AllElements[i].newDist,NumTargets);
        if ((currentDist < bestDist) || (currentDist == bestDist && currentNorm > bestNorm))
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            bestNorm = currentNorm;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist && currentNorm == bestNorm)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }

    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0,candidates.size()-1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}


// BP rand with the filter of choosing the nearest target
int A2(Element AllElements[], int counter)
{
    // Applying the Filter
    int nearest = 1; // change this to relax the filter of nearest targets
    int filter_dist; // keep track of the largest distance that will pass through the filter
    vector<int> filter_indices; // keep track of the indices of Target that satisfy the filter

    // sort the distance array
    int sorted_dist[NumTargets-TargetsFound];
    int next_index = 0;
    for (int i = 0; i < NumTargets; i++)
    {
        if (Dist[i] == 0) continue;
        sorted_dist[next_index++] = Dist[i];
    }
    sort(sorted_dist,sorted_dist+NumTargets-TargetsFound);

    // indices that can pass through the filter
    filter_dist = sorted_dist[min(nearest-1,NumTargets-TargetsFound-1)]; // largest distance that will pass through the filter
    for (int i = 0; i < NumTargets; i++)
    {

        if (Dist[i] <= filter_dist && Dist[i] > 0) 
        {
            filter_indices.push_back(i);
        }
    }

    // initialization
    int bestDist = LARGE;
    int currentDist;
    vector<int> candidates;


    for (int i = 0; i < counter; i++)
    {
        // Filtering
        if (!filtering(AllElements[i].newDist, filter_indices)) continue;
        // Normal BP rand
        currentDist = calculateDist(AllElements[i].newDist,NumTargets);
        if (currentDist < bestDist)
        {
            // Updating the best distance and norm
            bestDist = currentDist;
            // clear previous candidates
            candidates.clear();
            // inputting the new candidates
            candidates.push_back(i);
        }
        else if (currentDist == bestDist)
        {
            // equal candidates
            candidates.push_back(i);
        }
    }
    // randomly choose one of the candidates
    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0,candidates.size()-1);
    int rand_num = rand_distribution(rand_generator);
    return candidates[rand_num];
}
 */

bool filtering(int tempDist[], vector<int> filter_indices)
{
    for (int i = 0; i < filter_indices.size(); i++)
    {
        // if any of the acceptable distance is reduced, return true immediately
        if (tempDist[filter_indices[i]] < Dist[filter_indices[i]]) return true;
    }
    return false;
}

 
int calculateDist(int A[],int length)
{
    int s = 0;
    for (int i = 0; i < length; i++) s += A[i];
    return s;
}
int calculateNorm(int A[],int length)
{
    int s = 0;
    for (int i = 0; i < length; i++) s += A[i]*A[i];
    return s;
}

void refreshDist(){
   for(int k = 0; k < 32; k++) Dist[k] = InitDist[k];
}

void refreshDistAndTarget(){
    
    rand_generator.seed(time(0));

    shuffle(row.begin(), row.end(), rand_generator);
    shuffle(col.begin(), col.end(), rand_generator);

    vector<vector<int> > shuffled_mat = TargetMatrix;

    ShuffleFromIndices(shuffled_mat, row);

    for(int i = 0; i < NumTargets; i++) {
        ShuffleFromIndices(shuffled_mat[i], col);
    }

    direct.clear();
    
    int bit;
    int dj;
    int bitcnt = 0;
    for (int i = 0; i < NumTargets; i++) { //read row i
        long long int PowerOfTwo  = 1;
        Target[i] = 0;
        Dist[i] = -1;
        int cnt = 0;
        for (int j = 0; j < NumInputs; j++) {
            bit = shuffled_mat[i][j];
            if (bit) {
                Dist[i]++;
                Target[i] = Target[i] + PowerOfTwo;
                bitcnt += 1;
                dj = j;
            }
            PowerOfTwo = PowerOfTwo * 2;
        }
        // cout << cnt << "\n";
        if (bitcnt == 1) {
            // cout << dj << "\n";
            string ys, xs;
            ys = to_string(i);
            xs = to_string(dj);

            if (ys.length() < 2)
	            ys.insert(ys.begin(), 2 - ys.length(), '0');
            
            if (xs.length() < 2)
	            xs.insert(xs.begin(), 2 - xs.length(), '0');
            
            direct.push_back("y" + ys + " = " + "x" + xs);
        }
        bitcnt = 0;
    }

    // for(int k = 0; k < 32; k++) Dist[k] = InitDist[k];
}

void ReadTargetMatrix() {
    cin >> NumTargets;
    cin >> NumInputs;
    //check that NumInputs is < wordsize
    if (NumInputs >= 8*sizeof(long long int)) {
        cout << "too many inputs" << endl;
        exit(0);
    }

    row.resize(NumTargets);
    col.resize(NumInputs);
    iota(row.begin(), row.end(), 0);
    iota(col.begin(), col.end(), 0);
    TargetMatrix.resize(NumTargets);

    int bit;
    for (int i = 0; i < NumTargets; i++) { //read row i
        long long int PowerOfTwo  = 1;
        Target[i] = 0;
        Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
        TargetMatrix[i].resize(NumInputs);
        for (int j = 0; j < NumInputs; j++) {
            cin >> bit;
            TargetMatrix[i][j] = bit;
            if (bit) {
                Dist[i]++;
                Target[i] = Target[i] + PowerOfTwo;
            }
            PowerOfTwo = PowerOfTwo * 2;
        }
    }
    // Update the InitDist for subsequent rounds
    for (int k = 0; k < 32; k++) InitDist[k] = Dist[k];
}

bool is_base(long long int x) {
    //sanity check, shouldn't ask if 0 is base
    if (x==0) {
        cout << "asking if 0 is in Base " << endl;
        exit(0);
    }

    for (int i = 0; i < BaseSize; i++) {
        if (x == Base[i]) {
            return true;
        }
    }
    return false;
}

// Distance is 1 less than the number of elements
// in the base that I need to add in order to get Target[u].
// The next function calculates the distance from the base,
// augmented by NewBase, to Target[u]. Uses the following observations:
// Adding to the base can only decrease distance.
// Also, since NewBase is the sum of two old base
// elements, the distance from the augmented base
// to Target[u] can decrease at most by 1. If the
// the distance decreases, then NewBase must be one
// of the summands.

int NewDistance(int u) {
    //if Target[u] is in augmented base return 0;
    if (is_base(Target[u]) || (NewBase == Target[u])) {
        return 0;
    }

    // Try all combinations of Dist[u]-1 base elements until one sums
    // to Target[u] + NewBase. If this is true, then Target[u] is the
    // sum of Dist[u] elements in the augmented base, and therefore
    // the distance decreases by 1.

    if (reachable(Target[u] ^ NewBase, Dist[u] - 1, NumInputs)) {
        return (Dist[u]-1);
    } else {
        return Dist[u]; //keep old distance
    }
}

int NewDistanceXOR3(int u, Gate gate) {
    //if Target[u] is in augmented base return 0;
    if (is_base(Target[u]) || (NewBase == Target[u])) {
        return 0;
    }

    // Try all combinations of Dist[u]-1 base elements until one sums
    // to Target[u] + NewBase. If this is true, then Target[u] is the
    // sum of Dist[u] elements in the augmented base, and therefore
    // the distance decreases by 1.

    if (gate == Gate::XOR4) {
        if (reachable(Target[u] ^ NewBase, Dist[u] - 3, NumInputs)) {
           return (Dist[u]-3);
	} else if (reachable(Target[u] ^ NewBase, Dist[u] - 2, NumInputs)) {
           return (Dist[u]-2);
        } else if (reachable(Target[u] ^ NewBase, Dist[u] - 1, NumInputs)) {
            return (Dist[u]-1);
        }
    } else if (gate == Gate::XOR3) {
        if (reachable(Target[u] ^ NewBase, Dist[u] - 2, NumInputs)) {
           return (Dist[u]-2);
        } else if (reachable(Target[u] ^ NewBase, Dist[u] - 1, NumInputs)) {
            return (Dist[u]-1);
	}
    } else {
        if (reachable(Target[u] ^ NewBase, Dist[u] - 1, NumInputs)) {
           return (Dist[u]-1);
        }
    }
    return Dist[u]; //keep old distance
}


//return true if T is the sum of K elements among Base[S..BaseSize-1]
bool reachable(long long int T, int K, int S) {
    if (__builtin_popcount(T) <= K)
    {
        return true;
    }
    if (S > BaseSize-1) {
        return false; //exceeded count
    }

    if (K==0) {
        return false; //this is probably not reached
    }

    if (K==1) {
        for (int i=S; i < BaseSize; i++) if (T == Base[i]) {
            return true;
        }
        return false;
    }
    
    //consider those sums containing Base[S]
    if (reachable(T^Base[S], K-1, S+1)) {
        return true;
    }

    //consider those sums not containing Base[S]
    if (reachable(T, K, S+1)) {
        return true;
    }


    //not found
    return false;
}
