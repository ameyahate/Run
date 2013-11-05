/*
	Copyright (c) 2010 Cindy Hui 

	This file is part of DIFNet.

    DIFNet is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DIFNet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DIFNet.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <string>
#include "network.h"

using namespace std;

const int MAX = 51;				///< Default: Maximum of 50 node types
double sourceTrustTable[MAX][MAX];
double nodeTrustTable[MAX][MAX];

/**	Insert trust values into the map
 *	\param senderType
 *	\param receiverType
 *	\param value
 */
void insertTrustValue(int senderType, int receiverType, double value)
{
	nodeTrustTable[senderType][receiverType] = value;
}

/**	Returns the trust value between a sender and receiver
 *	\param senderType
 *	\param receiverType
 */
double getTrustValue(int senderType, int receiverType)
{
    return nodeTrustTable[senderType][receiverType];
}
	

/**	Main function
 *	\param parameter file
 * 	\param edgelist file
 *	\param average trust
 * 	\param network scenario
 *	\param trust differential e
 * 	\param output suffix
 * 	\param output folder
 * 	\param nodelist file
 * 	\param seed file
 */
map<string, int> setup(char* paramfile, char* edgefile, char* nodefile, double networkAvgTrust, int networkScenario, double networkTrustDiff, char * seedfile, bool show)
{
	Network graph;
	int u, v, w, x, y;
	double lp1, lp2, lq1, lq2;
	double high = 0.0;
	double low = 0.0;
	double LB, UB, sigma, val, pval, weight;
	double r;
	
	int totalRuns = 30;				// Number of time steps to run
	int contain_weights = 0;		// Edgefile contains weights or not
	int use_edge_weight = 0;		// Use the weights found in the edgefile or not

	int totalNodes = 0;				// Number of nodes + source nodes
	int numNodes =0;				// Number of nodes
	int numSources =0;				// Number of source nodes
	int numNodeTypes =0;			// Number of node types (characteristics)
	
	/********************************************//**
	* Read parameters from file
	* 
	***********************************************/



	//Source types
	map<int, int> sourceMap;		// Maps the source id to an internal source index
	int numSeedEdges[MAX];			// Number of seeds from each source
	int totalSeedEdges = 0;
	
	//Node types
	map<int,int> nodeList;			// Maps the internal node index to its node type characteristics
	int nodeGroups[MAX];			// Node group determines trust
	double nodeLB[MAX];				// Node lower bound
	double nodeUB[MAX];				// Node upper bound
	double psendList[MAX];			// Stores probability of successful communication psend
	double nodeSigma[MAX];			// Node spreading parameter
	int removedTimeList[MAX];		// Stores time between Believed and Removed states
	int abortTimeList[MAX];			// Stores times for spreading in Disbelieved state
	
	/********************************************//**
	* Read parameters from file
	* 
	***********************************************/
	


	ifstream input(paramfile);
	if(!input) {
		cerr << "Error opening " << paramfile << endl;
		exit(1);
	}	
	if (show)
		cout << "Read parameters from " << paramfile << endl;
	char c;
	input >> c;	
	while (!input.eof()){
		switch(c){
 			case 's':{		//source [id] [group] [info value] [psend] [evac] [activate] [edges]
				input >> u >> v >> val >> pval >> w >> x >> y;		
				
				// 	/********************************************//**
				// 	* Insert source nodes into Network graph
				// 	***********************************************/
				graph.nodeMap[u] = totalNodes;
 				sourceMap[u] = numSources;
 				graph.addSourceNode(u, u, v,val,pval,w, x);
				numSeedEdges[sourceMap[u]] = y;
				totalSeedEdges = totalSeedEdges + y;
				
				totalNodes++;
				numSources = numSources + 1;
				break;
 			}
			case 'v':		//node [type] [group] [LB] [UB] [sigma] [psend] [timer_removed] [timer]
				input >> u >> v >> LB >> UB >> sigma>> pval >> w >> x ;

				// 	/********************************************//**
				// 	* Store node parameters
				// 	***********************************************/
				nodeGroups[u] = v;
				nodeLB[u]=LB;
				nodeUB[u]=UB;
				nodeSigma[u]=sigma;
				psendList[u]=pval;
				removedTimeList[u]=w;
				abortTimeList[u]=x;

				numNodeTypes++;

				break;
			case 't':{				//trust values between node and source
				input >> u >> v >> val;
				sourceTrustTable[sourceMap[u]][v] = val; //node type v trusts source u with value val


				break;}
			case 'w':		//flag: edge file contains/does not contain weights
				input >> u;
				contain_weights = u;
				break;				
			case 'e':		//flag: use/ignore weights in edge file
				input >> u;
				use_edge_weight = u;
				break;
			case 'l':		//lambda parameter
				//set the value of lambda
				input >> lp1 >> lp2 >> lq1 >> lq2;
				if(show)
					cout << "lambda1: " << lp1 << "\tlambda2: " << lp2<< "\tlambdaAbort1: " << lq1 << "\tlambdaAbort2: " << lq2 << endl;
				break;
			case 'r':		//time steps
				//set the total number of time steps
				input >> totalRuns;
				break;
			case '%':		//comment
				input.ignore(1000, '\n');
				break;			
			default:
				break;
		}
		input >> c;
	}
	input.close();

	/********************************************//**
	* Read node list from file and insert into Network graph
	* 
	***********************************************/
	if (show)
		cout << "Read node list from " << nodefile << endl;

	ifstream input1(nodefile);
	if(!input1) {
		cerr << "Error opening " << nodefile << endl;
		exit(1);
	}		
	
 	input1 >> u >> v;
 	while (!input1.eof()){
		graph.nodeMap[u] = totalNodes;
		graph.addNode(u, v, nodeGroups[v], nodeLB[v],nodeUB[v],nodeSigma[v]);

		//set time to evac 
		graph.setNodeRemovedTime(u,removedTimeList[v]); 
		//set time for spread
		graph.setNodeAbortTime(u,abortTimeList[v]);			
		//set psend
		graph.setNodepSend(u,psendList[v]);			
		//set lambda parameters
		graph.setNodeLambdaParam1(u, lp1);
		graph.setNodeLambdaParam2(u, lp2);	
		numNodes++;	
 		totalNodes++;
 		nodeList[graph.nodeMap[u]] = v;
 		input1 >> u >> v;		
 	}
 	input1.close();	
 	if(show)
 	{
		cout << "nodeMap.size:\t" << (int)(graph.nodeMap).size() << endl;
		cout << "numNodes:\t" << numNodes << endl;
		cout << "numSources:\t" << numSources << endl;
		cout << "totalNodes:\t" << totalNodes << endl;
		cout << "totalSeedEdges:\t" << totalSeedEdges << endl;
 	}
	/********************************************//**
	* Count number of edges
	* 
	***********************************************/
	int edgesSame = 0;
	int edgesDiff = 0;
	val = 0.0;
	
	ifstream input4(edgefile);
	if(!input4) {
		cerr << "Error opening edges.txt" << endl;
		exit(1);
	}
	if(show)
		cout << "Read edges from " << edgefile << "..." ;
	if(contain_weights > 0){
		input4 >> u >> v >> weight;
		while (!input4.eof()){
			//read edges
			if(nodeGroups[nodeList[graph.nodeMap[u]]] == nodeGroups[nodeList[graph.nodeMap[v]]]){
				edgesSame++;
			} else {
				edgesDiff++;
			}	
			input4 >> u >> v >> weight;
		}		
	} else {
		input4 >> u >> v;	
		while (!input4.eof()){
			//read edges
			if(nodeGroups[nodeList[graph.nodeMap[u]]] == nodeGroups[nodeList[graph.nodeMap[v]]]){
				edgesSame++;
			} else {
				edgesDiff++;
			}			
			input4 >> u >> v;
		}	
	}
	input4.close();
	if(show)
		cout << "done"<< endl;
	
	int totalEdges = edgesSame + edgesDiff;
	//List of random edge weights
	vector<double> randomEdgeWeights;
	//Random weight
	double r1;
	//Counters
	int front_counter = 0;
	int end_counter = totalEdges - 1;
	
	if(show)
	{
		cout << "totalEdges:\t" << totalEdges << endl;
		cout << "edgesSame:\t" << edgesSame << endl;
		cout << "edgesDiff:\t" << edgesDiff << endl;
	}

	/********************************************//**
	* Assign edge weights
	* 
	***********************************************/
	double networkLowTrust = networkAvgTrust;
	double networkHighTrust = networkAvgTrust;
	double networkHighProp = 0.5;
	
	if(networkLowTrust < 0 ){
		networkLowTrust = 0;
	}
	if( networkHighTrust > 1){
		networkHighTrust = 1;
	}

	switch(networkScenario){
		case 0:	//Scenario - Equal trust
			if(show)
				cout << "Trust networkScenario 0"<< ": Equal trust between all nodes " << endl;

			networkLowTrust = networkAvgTrust;
			networkHighTrust = networkAvgTrust;
			if (show)
				cout << "Avg/Low/High:\t" << networkAvgTrust << "\t" << networkLowTrust << "\t" << networkHighTrust << endl;

			for(int i=1; i<=numNodeTypes; i++){
				for(int j=1; j<=numNodeTypes; j++){
					insertTrustValue(i,j,networkAvgTrust);
				}
			}
			
			break;
		case 1: //Scenario - High trust in same type
			cout << "Trust networkScenario 1"
			     << ": High trust in same group (High Trust = Avg + e; Low Trust = Avg - e)" << endl;
			networkLowTrust = networkAvgTrust - networkTrustDiff;
			networkHighTrust = networkAvgTrust + networkTrustDiff;
			if(networkLowTrust < 0 ){
				networkLowTrust = 0;
			}
			if( networkHighTrust > 1){
				networkHighTrust = 1;
			}
			cout << "Avg trust = " << networkAvgTrust << "\t" 
				 << (edgesSame*networkHighTrust+edgesDiff*networkLowTrust)/(edgesSame+edgesDiff) << endl;
			cout << "Avg/Low/High:\t" << networkAvgTrust << "\t" << networkLowTrust << "\t" << networkHighTrust << endl;	

			for(int i=1; i<=numNodeTypes; i++){
				for(int j=1; j<=numNodeTypes; j++){
					if( nodeGroups[i]==nodeGroups[j]){					
						insertTrustValue(i,j,networkHighTrust);
					} else {
						insertTrustValue(i,j,networkLowTrust);
					}
				}
			}		
			break;
		case 3: //Scenario - High trust in same type
			//Compute Low and High trust values based on number of edges
			if (show)
				cout << "Trust networkScenario 3"<< ": High trust in same group (High Trust = Avg + e; Low Trust computed) " << endl;
			if(edgesDiff > 0){
				networkLowTrust = networkAvgTrust - (edgesSame)*(networkTrustDiff)/(edgesDiff);
			} else {
				networkLowTrust = networkAvgTrust;
			}
			networkHighTrust = networkAvgTrust + networkTrustDiff;
			if(networkLowTrust < 0 ){
				networkLowTrust = 0;
			}
			if( networkHighTrust > 1){
				networkHighTrust = 1;
			}
			if (show)
				cout << "Avg/Low/High:\t" << networkAvgTrust << "\t" << networkLowTrust << "\t" << networkHighTrust << endl;
//			cout << e1 << "\t" << e2 << "\t" << low << "\t" << high << endl;			 	 
			for(int i=1; i<=numNodeTypes; i++){
				for(int j=1; j<=numNodeTypes; j++){
					if( nodeGroups[i]==nodeGroups[j]){					
						insertTrustValue(i,j,networkHighTrust);
					} else {
						insertTrustValue(i,j,networkLowTrust);
					}
				}
			}		
			break;			
		case 4: //Scenario - Randomly assign LOW and HIGH trust edges with prob 0.5
			cout << "Trust networkScenario 4" 
				 << ": Randomly distribute trust with prob 0.5" << endl;
			networkLowTrust = networkAvgTrust - networkTrustDiff;
			networkHighTrust = networkAvgTrust + networkTrustDiff;
			if(networkLowTrust < 0 ){
				networkLowTrust = 0;
			}
			if( networkHighTrust > 1){
				networkHighTrust = 1;
			}
			cout << "Avg trust = " << networkAvgTrust << "\t" 
				 << (edgesSame*high+edgesDiff*low)/(edgesSame+edgesDiff) << endl;
			cout << "Avg/Low/High:\t" << networkAvgTrust << "\t" << networkLowTrust << "\t" << high << endl;
			break;
		case 5: //Scenario - Randomly assign LOW and HIGH trust edges proportionally
			cout << "Trust networkScenario 5"
			  	 << ": Randomly distribute trust proportionally " << endl;	
			networkLowTrust = networkAvgTrust - networkTrustDiff;
			networkHighTrust = networkAvgTrust + networkTrustDiff;
			if(networkLowTrust < 0 ){
				networkLowTrust = 0;
			}
			if( networkHighTrust > 1){
				networkHighTrust = 1;
			}
			cout << "Avg trust = " << networkAvgTrust << "\t" 
				 << (edgesSame*high+edgesDiff*low)/(edgesSame+edgesDiff) << endl;
//			cout << "Compute proportion of high trust edges using edge counts\n";			
			networkHighProp = (double)edgesSame/(double)(edgesSame+edgesDiff);
			cout << "Avg/Low/High:\t" << networkAvgTrust << "\t" << networkLowTrust << "\t" << networkHighTrust << endl;
			cout << "Proportion of high trust edges: " << networkHighProp << endl;			
			break;	
		case 6: //Scenario - 
			cout << "Trust networkScenario 6"	
				 << ": Randomly assign trust values between range "
				 << "[low = " << networkLowTrust << ", high = " << networkHighTrust << "] uniform"<< endl;
//			cout << edgesSame << "\t" << edgesDiff << endl;
			//For average trust avg, we have the range [ low = avg - e1, high = avg + e1]
			//Randomly generate tn values in [low, high], where tn is the total # of edges
			for(int i=0;i<totalEdges;i++){

				//random double between 0 and 0.99999
				r1 = rand()/(double(RAND_MAX)+1);
				//random double between low and high
        		r1 = r1*(networkHighTrust-networkLowTrust)+networkLowTrust; 

				//Store values in randomEdgeWeights[totalEdges]				
				randomEdgeWeights.push_back(r1);
			}
			//sort the tn edge weights in ascending order
			sort(randomEdgeWeights.begin(), randomEdgeWeights.end());
			
			//random shuffle - weak tie values
			random_shuffle ( randomEdgeWeights.begin(), randomEdgeWeights.begin() + edgesDiff);
			//random shuffle - strong tie values
			random_shuffle ( randomEdgeWeights.end() - edgesSame, randomEdgeWeights.end() );
			
			//assign the first w edge weights to the weak ties, where w = edgesDiff
			//assign remaining edge weights, n - w, to the strong ties
			break;			
		default:
			break;
	}

// 	for(int i=1; i<=numNodeTypes; i++){
// 		for(int j=1; j<=numNodeTypes; j++){
// 			cout << i << "\t" << j << "\t" << getTrustValue(i,j) << endl;
// 		}
// 	}
	
	srand((unsigned)time(0));	
	
	/********************************************//**
	* Insert edges and weights into Network graph
	* 
	***********************************************/
	if(contain_weights > 0){
		if (show)
			cout << "Edge file contains weights; ";
		if(use_edge_weight > 0){
			if(show)
				cout << "Use edge weights" << endl;
		} else {
			if (show)
				cout << "Ignore edge weights" << endl;
		}
	}
		
	ifstream input2(edgefile);
	if(!input2) {
		cerr << "Error opening edges.txt" << endl;
		exit(1);
	}		
	if (show)
		cout << "Read edges from " << edgefile << "...";
	if(contain_weights > 0 && use_edge_weight > 0){			//edge list contains weights
		input2 >> u >> v >> weight;
		while (!input2.eof()){			
			graph.addEdge(u,v,  weight);
			graph.addEdge(v,u,  weight);		//comment out this line if edge list represented a directed graph

			input2 >> u >> v >> weight;
		}				
	} else {							//edge list does not contain weights
		if(contain_weights > 0 ){
			input2 >> u >> v >> weight;
		} else {
			input2 >> u >> v;	
		}
		while (!input2.eof()){

			switch(networkScenario){
				case 0:	//Scenario - Equal trust
					weight = networkAvgTrust;
					break;
				case 1: //High trust in same type
					weight = nodeTrustTable[nodeList[graph.nodeMap[u]]][nodeList[graph.nodeMap[v]]];
					break;
				case 3: //High trust in same type
					//Compute Low and High trust values based on number of edges		
					weight = nodeTrustTable[nodeList[graph.nodeMap[u]]][nodeList[graph.nodeMap[v]]];
					break;			
				case 4: //Randomly assign LOW and HIGH trust edges with prob 0.5
					r = rand()/(float)RAND_MAX;
					if(r<=0.5){
						weight = networkHighTrust;
					} else {
						weight = networkLowTrust;
					}
					break;
				case 5: //Randomly assign LOW and HIGH trust edges proportionally
					r = rand()/(float)RAND_MAX;
					if(r<=networkHighProp){
						weight = networkHighTrust;
					} else {
						weight = networkLowTrust;
					}
					break;	
				case 6: //Assign the first w edge weights to the weak ties, where w = edgesDiff
					//Assign remaining edge weights, n - w, to the strong ties
					if(nodeGroups[nodeList[graph.nodeMap[u]]] == nodeGroups[nodeList[graph.nodeMap[v]]]){
						//is a strong tie				
						weight = randomEdgeWeights.at(end_counter);
						end_counter--;
					} else {
						//is a weak tie
						weight = randomEdgeWeights.at(front_counter);
						front_counter++;
					}							
					break;			
				default:
					cout << "Wrong Network scenario : "<< networkScenario << endl;
					break;
			}	
	
			graph.addEdge(u,v, weight);
			graph.addEdge(v,u, weight);


			if(contain_weights > 0 ){
				input2 >> u >> v >> weight;
			} else {
				input2 >> u >> v;	
			}
		}		
	}
	input2.close();
	if (show)
		cout << "done"<< endl;
		
	/********************************************//**
	* Insert edges from source nodes into Network graph
	* 
	***********************************************/	

	// Read the edges from source nodes from seedfile
	if (show)
		cout << "Read edges from " << seedfile << endl;
	ifstream input3(seedfile);
	if(!input3) {
		cerr << "Error opening " << seedfile << endl;
		exit(1);
	}		
	int temp = 0;
	while (input3 >> u >> v){
		
//		if(numSeedEdges[sourceMap[u]]>0){
			weight = 0.9;
			graph.addEdge(u,v, weight);	
			temp++;
//			cout << u << "\t" << v << "\t>\t"<< graph.nodeMap[u] << "\t" << graph.nodeMap[v] << "\t" << weight << endl;
//			cout << u << "\t" << v << "\t" << weight << endl;
//			numSeedEdges[sourceMap[u]]--;
//		}
	}
	input3.close();
	if (show)
		cout << "Seeds:\t" << temp << endl;
	
// 	cout << 100001 << "\t" << graph.nodeMap[100001] << endl;
// 	cout << 100002 << "\t" << graph.nodeMap[100002] << endl;
// 	cout << 1 << "\t" << graph.nodeMap[1] << endl;
// 	cout << 2 << "\t" << graph.nodeMap[2] << endl;
// 	cout << 3 << "\t" << graph.nodeMap[3] << endl;
// 	cout << 4 << "\t" << graph.nodeMap[4] << endl;
	


	
//	/********************************************//**
//	* Run
//	*
//	***********************************************/
	return (graph.run(totalRuns));
}



int main(int argc, char *argv[])
{
	int iterations = 2;
	//command line arguments
	char * paramfile = "param.txt";			//parameter file
	char * edgefile = "edges.txt";			//edgelist file
	char * nodefile = "nodes.txt";			//node file
	double networkAvgTrust = 0.0; 			//network average trust
	int networkScenario = 0;				//trust network scenario
	double networkTrustDiff = 0.0;			//trust differential
	char * seedfile = "seed.txt";			//seed file
	char * outputfolder = "./";

	if(argc>=9){
		paramfile = argv[1];
		edgefile = argv[2];
		nodefile = argv[3];
		networkAvgTrust = atof(argv[4]);
		networkScenario = atoi(argv[5]);
		networkTrustDiff = atof(argv[6]);
		seedfile = argv[7];
		outputfolder = argv[8];
	} else {
		cout << "Usage:\n";
		cout << "  ./main argv[1]... argv[9]" << endl;
		cout << "argv[1] = parameter file" << endl;
		cout << "argv[2] = edge file" << endl;
		cout << "argv[3] = node file " << endl;
		cout << "argv[4] = average trust (e.g. 0.50)" << endl;
		cout << "argv[5] = trust network scenario" << endl;
		cout << "argv[6] = trust differential e (e.g. 0.05)" << endl;
		cout << "argv[7] = seed file " << endl;
		cout << "argv[8] = output folder " << endl;
		cout << " Scenarios: " << endl
			 << "   0: Homogeneous trust (Avg trust between all nodes) " << endl
		     << "   1: High trust in same group (High Trust = Avg + e; Low Trust = Avg - e)" << endl
		     << "   3: High trust in same group (High Trust = Avg + e; Low Trust computed) " << endl
		     << "   4: Randomly distribute trust (High Trust = Avg + e; Low Trust = Avg - e) with prob 1/2 " << endl
		     << "   5: Randomly distribute trust (High Trust = Avg + e; Low Trust = Avg - e) proportionally " << endl
			 << "   6: Randomly assign trust values between range [Avg - e, Avg + e] uniform"<< endl;

		exit(1);
	}


	/********************************************//**
	* Construct identifier for output files
	***********************************************/
	char * ident = NULL;
	char pident[50];
//	cout << "edgefile: " << edgefile << endl;
//	cout << "paramfile: " << paramfile << endl;

	string s, s1, s2;
	size_t s_start, s_end;

	s  = paramfile;
	s_start = s.find_last_of("/") + 7;
	s_end = s.find_last_of(".");
	s1 = s.substr(s_start, s_end - s_start);

	s = seedfile;
	s_start = s.find_last_of("/") + 1;
	s_end = s.find_last_of(".");
	s2 = s.substr(s_start, s_end - s_start);

	// param-seedfile-networkAverageTrust-networkScenario-networkTrustDiff*10-suffix
	sprintf(pident, "%s-%s.txt",s1.c_str(),s2.c_str());
	ident=pident;
//	cout << "Identifier: " << ident << endl;

	//Run 100 times, store values and print average
	map<string, int> tempmap;


	double avg = 0.0;
	bool show = false;
	for(int count = 1; count <= iterations; count++)
	{
		tempmap = setup(paramfile, edgefile, nodefile, networkAvgTrust, networkScenario, networkTrustDiff, seedfile, show);
		avg += (double)tempmap["evacuated"] + (double)tempmap["believer"];
	}

	avg = avg/iterations;

	ofstream output(ident,ios::app);
	if(!ident)
	{
		cerr<< "Error opening file: "<< ident<<endl;
		exit(1);
	}
	output << avg << endl;
	output.close();
}



