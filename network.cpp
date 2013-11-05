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
#include "network.h"
using namespace std;

/**	Add node to nodeSet list
 *	\param nid
 * 	\param ntype
 *	\param ngroup
 *	\param nLB
 *	\param nUB
 *	\param nsigma
 */
void Network::addNode(int nid, int ntype, int ngroup, double nLB, double nUB, double nsigma)
{
	Node * newNode = new Node(nid,ntype,ngroup,nLB,nUB,nsigma);
//	cout << nid << "\t" << ntype << "\t" << nodeMap[nid] << endl;
	nodeSet.push_back(newNode);
}

/**	Return a pointer to the Node with the specified id
 *	\param id
 */
Node * Network::getNode(int id)
{
	return nodeSet[nodeMap[id]];
}

/*	Set time steps between Believed state and Removed state for Node with the specified id
 *	\param id
 *	\param etime
 */
void Network::setNodeRemovedTime(int id, int etime)
{
	Node * temp = getNode(id);
	temp->setRemovedTime(etime);
}

/*	Set time steps for Disbelieved node to spread Abort message for Node with the specified id
 *	\param id
 *	\param stime 
 */
void Network::setNodeAbortTime(int id, int stime)
{
	Node * temp = getNode(id);
	temp->setAbortTime(stime);
}

/**	Set probability of successfully communicating with others for Node with the specified id
 *	\param id
 *	\param p 
 */
void Network::setNodepSend(int id, double p)
{
	Node * temp = getNode(id);
	temp->setpSend(p);
}

/**	Set lambda parameter used in computeInfoUnionMax for Node with the specified id
 *	\param id
 *	\param lp 
 */
void Network::setNodeLambdaParam1(int id, double lp)
{
	Node * temp = getNode(id);
	temp->setLambdaParam1(lp);
}

/**	Set lambda parameter used in computeInfoFusion for Node with the specified id
 *	\param id
 *	\param lp 
 */
void Network::setNodeLambdaParam2(int id, double lp)
{
	Node * temp = getNode(id);
	temp->setLambdaParam2(lp);
}

/**	Set lambda parameter used in computeInfoUnionMax for Node with the specified id
 *	\param id
 *	\param lq
 */
void Network::setNodeLambdaAbortParam1(int id, double lq)
{
	Node * temp = getNode(id);
	temp->setLambdaAbortParam1(lq);
}

/**	Set lambda parameter used in computeAbortFusion for Node with the specified id
 *	\param id
 *	\param lq
 */
void Network::setNodeLambdaAbortParam2(int id, double lq)
{
	Node * temp = getNode(id);
	temp->setLambdaAbortParam2(lq);
}

/**	Return the number of nodes in nodeSet	
 */
int Network::getNumNodes()
{
	return (int)nodeSet.size();
}

/**	Add a source node
 *	\param nid
 * 	\param sid
 *	\param sgroup
 *	\param svalue
 *	\param sactive
 */
void Network::addSourceNode(int sid, int stype, int sgroup, double svalue,  double pval, int sremove, int sactive)
{
	Node * newNode = new Node(sid,stype,sgroup,svalue, pval, sremove, sactive);
//	cout << sid << "\t" << stype << "\t" << newNode->getValue()<< endl;
	nodeSet.push_back(newNode);
}

/**	Remove all incoming edges to node with the specified id
 *	\param nid
 */
void Network::removeInEdgesNode(int nid)
{
	Node * sNode = getNode(nid);
	sNode->removeInEdges();
}

/**	Add an edge from node uid to vid with weight wgt
 *	\param uid
 * 	\param vid
 *	\param wgt
 */
void Network::addEdge(int uid, int vid, double wgt)
{
	Node * u = getNode(uid);
	Node * v = getNode(vid);
	
	u->neighbors.push_back(v);
	u->neighbors_weights.push_back(wgt);
	u->received_from.push_back(0);
	edgeSet.insert(pair<int,int>(uid,vid));
}

/**	Run 
 *	\param outputfolder
 *	\param ident
 *	\param numtime number of time steps to run
 */
map<string,int>  Network::run(int numtime)
{
//	char cFileName[100];
	
//	sprintf(cFileName, "%soutput-%s-counts.txt", outputfolder, ident);
//	FILE * cOut = fopen(cFileName,"w");
//	fprintf(cOut, "#Time\tNonInf\tNonbel\tUndec\tBeliev\tEvac\n");

	int timestep=0;
	srand((unsigned)time(0));
	
	//printNodes(0);
	//printEdges(0);
	//printGraph(0);
//	printCounts(cOut,0);
	bool show = false;
		
    for(timestep=1; timestep<=numtime; timestep++){
		if(show)
			cout << "Time step: "<< timestep << endl;
		for(int i=0; i<(int)nodeSet.size(); i++){
			if(show)
				cout << "nodeID: " << nodeSet[i]->getID() << "\t Type: " <<  nodeSet[i]->getType() << "\t Status: " << nodeSet[i]->getStatus() << endl;			
			if((nodeSet[i])->getStatus() == 3){	//is a source node
				if((nodeSet[i])->getActivate() == timestep){
					if((nodeSet[i])->getGroup() == 1){
						(nodeSet[i])->sendMessage(1);
						if(show)
							cout << (nodeSet[i])->getType() << " send action msg" << endl;						
					} else if ((nodeSet[i])->getGroup() == 2){
						(nodeSet[i])->sendMessage(2);
						if(show)
							cout << (nodeSet[i])->getType() << " send abort msg" << endl;						
					} else {}
					(nodeSet[i])->setStatus(-3);
				}
			}
			//check status, if status is 2 (believer), send message to neighbors
			else if((nodeSet[i])->getStatus() == 2 ){
				if(nodeSet[i]->timeRemoved <= 0){		//time to evacuate
					nodeSet[i]->setStatus(-2);
					nodeSet[i]->removeInEdges();
					nodeSet[i]->removeOutEdges();
					if(show)	
						cout << "Node " << nodeSet[i]->getID() << " has evacuated w/ status " << nodeSet[i]->getStatus() << endl;

				} else {				
					(nodeSet[i])->sendMessage(1);
					nodeSet[i]->timeRemoved -= 1;
				}
			}
			//check status, if status is 1 (undecided), query neighbors
			else if((nodeSet[i])->getStatus() == 1 ){
				(nodeSet[i])->queryNeighbors();
			}
			//check status, if status is 0 (disbelieved) and spreadAbortflag is TRUE, 
			//send abort message to neighbors
			else if((nodeSet[i])->getStatus() == 0 && (nodeSet[i])->spreadAbortflag == 1){
				if(nodeSet[i]->timeAbort <= 0){		//time to spread abort
				} else {
					(nodeSet[i])->sendMessage(2);
					nodeSet[i]->timeAbort -= 1;
				}
			}
		}

//		cout << "\nCompute node Information Fused Value and update node Status" << endl;
		for(int i=0; i<(int)nodeSet.size(); i++){
			if((abs((nodeSet[i])->getStatus()) != 3) && ((nodeSet[i])->getStatus() >-2)){
				//is not a source node and has not evacuated
				(nodeSet[i])->computeInfoFusion();
				(nodeSet[i])->computeAbortFusion();
				(nodeSet[i])->updateStatus();
			}
		}

//		printNodes(timestep);		
//		printEdges(timestep);
		//printGraph(timestep);
//		printCounts(cOut, timestep);
	}
//	cout << "Wrote to " << cFileName << endl;
//	printNodesSummaryByGroup(outputfolder, ident);
//	printNodesDetailed(outputfolder,ident);
//	fclose(cOut);

	int status;
	map<string, int> states;

	states["noninformed"]=0;
	states["nonbeliever"]=0;
	states["undecided"]=0;
	states["believer"]=0;
	states["evacuated"]=0;
	for(int i=0; i<(int)nodeSet.size(); i++){
		status = (nodeSet[i])->getStatus();
		if(status==-1){
			states["noninformed"]++;
		} else if(status==0){
			states["nonbeliever"]++;
		} else if (status==1){
			states["undecided"]++;
		} else if(status==2){
			states["believer"]++;
		} else if(status==-2){
			states["evacuated"]++;
		}
	}
	return states;
}

/**	Print a list of node id, type, group, status, and fused values at the given timestep
 * \param timestep
 */
void Network::printNodes(int timestep)
{
	char nFileOut[20];
	sprintf(nFileOut,"output-states-%02d.txt",timestep);
	FILE * nOut = fopen(nFileOut,"w");

	int id, status, type, group;
	double value, abortvalue;
	fprintf(nOut,"NodeID\tType\tGroup\tStatus\tFusedValue\tAbortFusedValue\n");
	for(int i=0; i<(int)nodeSet.size(); i++){
		id = (nodeSet[i])->getID();
		type = (nodeSet[i])->getType();
		group = (nodeSet[i])->getGroup();
		status = (nodeSet[i])->getStatus();
		value = (nodeSet[i])->getValue();
		abortvalue = (nodeSet[i])->getAbortValue();		
		fprintf(nOut, "%2d\t%2d\t%2d\t%2d\t%ft%f\n", id, type, group, status, value, abortvalue);
	}
	fprintf(nOut,"\n");
}

/**	Print number of nodes in each state by node group
 *	\param outputfolder
 *	\param ident
 */
void Network::printNodesSummaryByGroup(char * outputfolder, char * ident)
{
	char nFileName[100];
	sprintf(nFileName, "%soutput-%s-nodes-summary.txt",outputfolder, ident);
	FILE * out = fopen(nFileName,"w");

	int MAX = 5;							//max number of groups of nodes
	int nodegroupstates[MAX][5];			//number of nodes of each group in each state
	for(int i=0;i<=MAX;i++){
		for(int j=0;j<5;j++){
		nodegroupstates[i][j]=0;	
		}
	}
	int status, type, group;
	fprintf(out,"NodeGroup\tStatus\tCount\n");
	for(int i=0; i<(int)nodeSet.size(); i++){
		status = (nodeSet[i])->getStatus();
		if((nodeSet[i])->getStatus() > -3 && (nodeSet[i])->getStatus() < 3 ){	
			type = (nodeSet[i])->getType();
			group = (nodeSet[i])->getGroup();
			nodegroupstates[group][status+2]++;
		}
	}
	
	for(int i=1;i<MAX;i++){
		int Unin = nodegroupstates[i][1];
		int Disb = nodegroupstates[i][2];
		int Unde = nodegroupstates[i][3];
		int Beli = nodegroupstates[i][4];
		int Evac = nodegroupstates[i][0];
		fprintf(out, "%2d\t%2d\t%2d\t%2d\t%2d\t%2d\n", i, Unin, Disb, Unde, Beli, Evac );			
	}	
		

	cout << "Wrote to " << nFileName << endl;
}
/**	Print a list of node id, type, group, status, fused value
 * 	\param outputfolder
 *	\param ident
 */
void Network::printNodesDetailed(char * outputfolder, char * ident)
{
	char nFileName[100];
	sprintf(nFileName, "%soutput-%s-nodes-detailed.txt",outputfolder, ident);
	
	FILE * out = fopen(nFileName,"w");
	
	int id, status, type, group, active, degree;
	double value, abortvalue, psend, lambda1, lambda2, lb, ub, timeRemoved;
cout << "nodeSet.size: " << (int)nodeSet.size() << endl;
//	fprintf(out,"NodeID\tType\tDegree\tLB\tUB\tStatus\tFusedValue\tActive\tTime\tpSend\tlambda1\tlambda2\n");
//	fprintf(out,"NodeID\tType\tStatus\tFusedValue\tTime\n");
	fprintf(out,"NodeID\tType\tGroup\tStatus\tFusedValue\tAbortFusedValue\n");
	for(int i=0; i<(int)nodeSet.size(); i++){
		status = (nodeSet[i])->getStatus();
		//print only the Uninformed, Disbelieved, Undecided, and Believed
//		if(status >= -1 && status <= 2){
			id = (nodeSet[i])->getID();
			type = (nodeSet[i])->getType();
			group = (nodeSet[i])->getGroup();
			degree = (int)(nodeSet[i]->neighbors).size();
			lb = (nodeSet[i])->getLB();
			ub = (nodeSet[i])->getUB();
			value = (nodeSet[i])->getValue();
			abortvalue = (nodeSet[i])->getAbortValue();			
			timeRemoved = (nodeSet[i])->timeRemoved;
			psend = (nodeSet[i])->psend;
			lambda1 = (nodeSet[i])->lambda1;
			lambda2 = (nodeSet[i])->lambda2;
			active = (nodeSet[i])->timeActivate;
//			fprintf(out, "%2d\t%2d\t%2d\t%0.2f\t%0.2f\t%2d\t%f\t%2d\t%3.0f\t%.2f\t%.2f\t%.2f\n", id, type, degree, lb, ub, status, value, active, timeRemoved, psend, lambda1, lambda2);
//			fprintf(out, "%2d\t%2d\t%2d\t%f\t%3.0f\n", id, type, status, value, timeRemoved);
//			cout << i << "\t" << id << "\t" << type << "\t" << status << endl;
			fprintf(out, "%2d\t%2d\t%2d\t%2d\t%f\t%f\n", id, type, group, status, value, abortvalue);
//		}
	}
	cout << "Wrote to " << nFileName << endl;
}

/**	Print a list of all the edges and its trust weight at the given timestep
 * \param timestep
 */
void Network::printGraph(int timestep)
{
	char eFileOut[20];
	sprintf(eFileOut,"output-graph-%02d.txt",timestep);
	FILE * eOut = fopen(eFileOut,"w");
	for(int i=0; i<(int)nodeSet.size(); i++){
	 	// Print nodes that are not sources
		if(nodeSet[i]->getStatus()<3 && nodeSet[i]->getStatus()>-3){

		int uid, vid;
		double w;
		uid = nodeSet[i]->getID();
		//cout << "uid " << "\t" << uid << "\t" << (nodeSet[i]->neighbors).size() <<endl;
		for(int j=0; j<(int)(nodeSet[i]->neighbors).size();j++){
			if((nodeSet[i]->neighbors).at(j)!=NULL){
				vid = ((nodeSet[i]->neighbors).at(j))->getID();
				//fprintf(eOut, "%d\t%d\n", uid, vid);
				//w is the trust from u to v
				w = (nodeSet[i]->neighbors_weights).at(j);
				fprintf(eOut, "%d\t%d\t%0.2f\n", uid, vid, w);

			}
		}
		}
	}
//	fprintf(eOut,"\n");
}


/**	Print a list of all the active edges between informed nodes at the given timestep
 * 	\param timestep
 */
void Network::printEdges(int timestep)
{
	char eFileOut[20];
	sprintf(eFileOut,"output-edges-%02d.txt",timestep);
	FILE * eOut = fopen(eFileOut,"w");
	for(int i=0; i<(int)nodeSet.size(); i++){
	 	// Print nodes that were informed
		if(nodeSet[i]->getStatus()<3 && nodeSet[i]->getStatus()>-1){
		
		int uid, vid;
		//double w;
		uid = nodeSet[i]->getID();
		//cout << "uid " << "\t" << uid << "\t" << (nodeSet[i]->neighbors).size() <<endl;
		for(int j=0; j<(int)(nodeSet[i]->neighbors).size();j++){

			if((nodeSet[i]->neighbors).at(j)!=NULL){

				int neighbor_status = ((nodeSet[i]->neighbors).at(j))->getStatus();

				if((neighbor_status<3 && neighbor_status>-1)){
					vid = ((nodeSet[i]->neighbors).at(j))->getID();
					fprintf(eOut, "%d\t%d\n", uid, vid);
				}
			}
		}
		}
	}
//	fprintf(eOut,"\n");
}

/**	Creates a .dot file for the graph at the specified time step
 * \param timestep
 */
void Network::printDotGraph(int timestep)
{
	int id, status;
	const char * color;
	char gFileOut[20];
	sprintf(gFileOut,"graph%d.dot",timestep);

	//cout << "Create .dot file for timestep "<< timestep << endl;

	FILE * gOut = fopen(gFileOut,"w");
	fprintf(gOut,"digraph graph%d {\n",timestep);
//	fprintf(gOut,"graph graph%d {\n",timestep);	
	fprintf(gOut,"graph [normalize=true, size=\"6.5,10\"]\n");
	fprintf(gOut,"node [style=\"filled\"]\n");
	for(int i=0; i<(int)nodeSet.size(); i++){
		id = (nodeSet[i])->getID();
		status = (nodeSet[i])->getStatus();
		switch(status){
			case -1:
				color = "white";			//Noninformed
				break;
			case 0:
				color = "grey";				//Nonbeliever
				break;
			case 1:
				color = "green";			//Undecided
				break;
			case 2:
				color = "blue";				//Believer
				break;
			case -2:
				color = "black";			//Evacuated
				break;
			case 3:
				color = "red";				//Source
				break;
			case -3:
				color = "red";				//Removed source
				break;
			default:
				color = "white";

		}
	//	fprintf(gOut, "N%d [label=\"%d\" fillcolor=\"%s\"]\n", id, id, color);
	//	fprintf(gOut, "N%d [label=\"%d\" shape=\"circle\" fillcolor=\"%s\"]\n", id, id, color);
	//	fprintf(gOut, "N%d [label=\"%d\" shape=\"point\" fillcolor=\"%s\"]\n", id, id, color);	
		fprintf(gOut, "N%d [label=\"%d\" shape=\"point\" height=\"0.20\" fillcolor=\"%s\"]\n", id, id, color);		
	}
/*	for(int i=0; i<(int)nodeSet.size(); i++){
		uid = nodeSet[i]->getID();
		for(int j=0; j<(int)(nodeSet[i]->neighbors).size();j++){
			if((nodeSet[i]->neighbors).at(j)!=NULL){
				vid = ((nodeSet[i]->neighbors).at(j))->getID();
				fprintf(gOut, "N%d->N%d\n", uid, vid);
			}
		}
	}
*/

	multimap<int,int>::iterator p=edgeSet.begin();
	while(p!=edgeSet.end()){
		if((p->first)!=101 && (p->second)!=101){
		fprintf(gOut, "N%d->N%d [arrowhead = none] \n",p->first,p->second);
//		fprintf(gOut, "N%d--N%d  \n",p->first,p->second);		
		}
		p++;
	}
	fprintf(gOut,"}\n");
}

/**	Print number of Uninformed, Disbelieved, Undecided, Believed, and Evacuated nodes at given time step
 * 	\param out
 * 	\param time
 */
void Network::printCounts(FILE * out, int time)
{
	int status;
	int noninformed=0;
	int nonbeliever=0;
	int undecided=0;
	int believer=0;
	int evacuated=0;
	for(int i=0; i<(int)nodeSet.size(); i++){
		status = (nodeSet[i])->getStatus();
		if(status==-1){
			noninformed++;
		} else if(status==0){
			nonbeliever++;
		} else if (status==1){
			undecided++;
		} else if(status==2){
			believer++;
		} else if(status==-2){
			evacuated++;
		}
	}
	fprintf(out, "%d\t%d\t%d\t%d\t%d\t%d\n", time, noninformed, nonbeliever, undecided, believer, evacuated);
}

