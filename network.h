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
#include <vector>
#include <map>
#include "node.h"

class Network
{
	public:
		Network(){};
		~Network(){};

		//Add regular node
		void addNode(int nid, int ntype, int ngroup, double nLB, double nUB, double nsigma);
		//Add source node
		void addSourceNode(int sid, int stype, int sgroup, double svalue,  double pval, int sremove, int sactive);
		
		void addEdge(int uid, int vid, double wgt = 0.0);
		void removeInEdgesNode(int nid);			//remove all incoming edges

		Node * getNode(int id);						//return pointer to node with id
		int getNumNodes();							//return total number of nodes

		map<string,int>  run(int numtime);	//run

		//Print a list of node id, type, group, status, and fused values at the given timestep
		void printNodes(int timestep);
		//Print a list of all the active edges between informed nodes at the given timestep
		void printEdges(int timestep);
 		// Print a list of all the edges and its trust weight at the given timestep
		void printGraph(int timestep);
		//make .dot file for given timestep
		void printDotGraph(int timestep);	

 		//Print number of nodes in each state by node group
 		void printNodesSummaryByGroup(char * outputfolder, char * ident);
 		
 		//Print a list of node id, type, group, status, fused value
		void printNodesDetailed(char * outputfolder, char * ident);		
		
		//Print number of Uninformed, Disbelieved, Undecided, Believed, and Evacuated nodes at given time step
		void printCounts(FILE * out,int time);

		//set time between believed state and removed state
		void setNodeRemovedTime(int ident, int etime);
		//set time for spreading in disbelieved state
		void setNodeAbortTime(int ident, int stime);
		
		int removedTime;				/**< time between Believed state and Removed */
		int abortTime;					/**< time for spreading in Disbelieved state */

		//set probability of successfully communicating with others
		void setNodepSend(int id, double p);	
		
		//Lambda parameters used in computeInfoUnionMax and computeInfoFusion
		void setNodeLambdaParam1(int id, double lp);
		void setNodeLambdaParam2(int id, double lp);
		//Lambda parameters used in computeInfoUnionMax and computeAbortFusion
		void setNodeLambdaAbortParam1(int id, double lq);
		void setNodeLambdaAbortParam2(int id, double lq);
		
		map<int, int> nodeMap; 		/**< Maps the node ids to their index in nodeSet */
		vector<Node *> nodeSet;		/**< store node list in a vector */
		multimap<int,int> edgeSet;	/**< store edge list (used for graph) in a multimap */
};

