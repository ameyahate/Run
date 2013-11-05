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
#include <iomanip>
#include <list>
#include <vector>
#include <map>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "node.h"
using namespace std;

/**	Returns the trust value for edge between node and its neighbor
 *	\param receiver
 */
double Node::getTrustValue(int receiver)
{
	return neighbors_weights.at(receiver);
}

/**	Set time steps between Believed state and Removed state
 *	\param etime time steps
 */
void Node::setRemovedTime(int etime)
{
    //timeRemoved = getPoissonRN(1.0/etime);
	timeRemoved = etime;
	timeReset = etime;
}

/**	Set time steps for Disbelieved node to spread Abort message
 *	\param stime time steps to spread Abort message
 */
void Node::setAbortTime(int stime)
{
	timeAbort = stime;
}

/**	Set node's probability of successfully communicating with others 
 *	\param p value between 0 and 1
 */
void Node::setpSend(double p)
{
    psend = p;
}

/**	Remove all incoming edges to this node
 * 	
 */
void Node::removeInEdges()
{
    //cout << "remove incoming edges to " << getID() << endl;
    for(int i=0; i<(int) neighbors.size(); i++){
        Node * temp = neighbors[i];
        if(temp!=NULL){
            ////cout << "Neighbors: " << temp->getID() << endl;
            for(int j=0; j<(int) (temp->neighbors).size(); j++){
                if(temp->neighbors[j]!=NULL){
                    int t = (temp->neighbors[j])->getID();
                    if(t==getID()){
                        temp->neighbors[j] = NULL;
                        ////cout << "removed " << t <<  endl;
                    }
                }
            }
        }
    }

}

/**	Remove all outgoing edges from this node
 * 	
 */
void Node::removeOutEdges()
{
    //cout << "remove outgoing edges from " << getID() << endl;
    for(int i=0; i<(int) neighbors.size(); i++){
        neighbors[i] = NULL;
    }
}

/**	Send message to neighbors with probability psend
 * 	\param msgtype
 */
void Node::sendMessage(int msgtype)
{
	double r=0.0;			//random number
    for(int i=0; i<(int) neighbors.size(); i++){
        Node * toNode = neighbors[i];
        if(toNode!= NULL){
        	if(toNode->getStatus() > -2){        	
				r = rand()/(float)RAND_MAX;			
				if(r<psend){
					//was able to reach neigbor node
					//sender computes info update and sends sourceSet data structure to receiver
					sourceMap sendSet;
					computeInfoUpdate(sendSet, msgtype, i);
					toNode->receiveMessage(sendSet, msgtype);
				} else {
					//was unable to reach neighbor node
				}
			}
        }
    }
}

/**	Send message to the specified node
 * 	\param msgtype
 *	\param nodeID nodeID of the recipient node
 */
void Node::sendMessageToNode(int msgtype, int nodeID)
{
    for(int i=0; i<(int) neighbors.size(); i++){
		Node * toNode = neighbors[i];
		if(toNode!=NULL){
			if(toNode->getID() == nodeID && toNode->getStatus() > -2){
				//sender computes info update and sends sourceSet data structure to receiver
				sourceMap sendSet;
				computeInfoUpdate(sendSet, msgtype, i);
				toNode->receiveMessage(sendSet, msgtype);
				break;
			}
		}
	}
}

/**	Receive message. 
 *	If node was Uninformed, then insert message set into sourceSet. 
 *	If node was previously informed, then call computeInfoUnionMax function.
 *	\param senderSet
 *	\param msgtype
 */
void Node::receiveMessage(sourceMap senderSet, int msgtype)
{

//	cout << "Node " << getID() << " received message \n";
	if(msgtype==1){
		receivedflag = 1;
	} else if (msgtype ==2){
		receivedAbortflag = 1;
	} else {}
	
    if(getStatus()==-1){        //is Uninformed
        sourceMap::iterator iter = senderSet.begin();
        while(iter!=senderSet.end()){
			if(msgtype==1){
            	sourceSet.insert(sPair(iter->first,iter->second));
            } else if (msgtype==2){
            	abortSet.insert(sPair(iter->first, iter->second));
            } else {}
            iter++;
        }
    } else {                    //was informed before
		computeInfoUnionMax(senderSet,msgtype);
    }
	hasUpdated = 1;
}

/**	Query neighbors for their set of source value pairs 
 *	with probability of successful communication psend
 */
void Node::queryNeighbors()
{
	double r=0.0;			//random number

//	cout << "QUERY: Node " << getID() << " query neighbors" << endl;
    for(int i=0; i<(int) neighbors.size(); i++){
        Node * askNode = neighbors[i];
		if(askNode!= NULL){

			r = rand()/(float)RAND_MAX;
			if(r<psend){
				if(askNode->getStatus() == -2){
					//Node has already evacuated
//					cout << "QUERY: Neighbor Node " << askNode->getID() << " has already evacuated" << endl;					
				} else if(askNode->getStatus() == -1){
					//Node has not been informed
//					cout << "QUERY: Neighbor Node " << askNode->getID() << " is uninformed" << endl;
				} else if(askNode->getStatus() == 0 ){
					//Node is a disbelieved
//					cout << "QUERY: Neighbor Node " << askNode->getID() << " is disbelieved" << endl;
					if((askNode->hasUpdated == 1) || received_from[i] == 0){
						if(!askNode->sourceSet.empty()){
							askNode->sendMessageToNode(1, getID());
						}

						if(!askNode->abortSet.empty()){
							askNode->sendMessageToNode(2, getID());
						}
						received_from[i] = 1;
					}
				
				} else if(askNode->getStatus() == 1){
					//Node is a undecided
					//Ask the neighbor node to send a message to this node
//					cout << "QUERY: Neighbor Node " << askNode->getID() << " is undecided" << endl;
					if((askNode->hasUpdated == 1) || received_from[i] == 0){
						if(!askNode->sourceSet.empty()){
							askNode->sendMessageToNode(1, getID());
						}
						if(!askNode->abortSet.empty()){
							askNode->sendMessageToNode(2, getID());
						}
						received_from[i] = 1;
					}
				} else if(askNode->getStatus() == 2){
					//Node is believed
					//Ask the neighbor node to send a message to this node
//					cout << "QUERY: Neighbor Node " << askNode->getID() << " is believed" << endl;
					askNode->sendMessageToNode(1, getID());
				}
			} else {
				//was unable to reach neighbor
//				cout << "QUERY: Node " << getID() << " was unable to reach Node " << askNode->getID() << endl;
			}
		}
	}
}

/**	Recipient node: Fused information values for a singleSource using parameter lambda
 *	\param infoSet
 *	\param lambda
 */
double Node::fusedSingleSource(valueMap infoSet, double lambda)
{
    valueMap::iterator iterS = infoSet.begin();
	double sum = 0.0;
	double max = 0.0;
	double combinedValue = 0.0;
	while(iterS != infoSet.end()){

		sum = sum + iterS->second;	//sum of info value for source s
		if(iterS->second > max){
			max = iterS->second;
		}
		iterS++;
	}

	combinedValue=lambda*sum+(1.0-lambda)*max;		//compute val
	if(combinedValue>1.0){							//take min of 1 and val
		combinedValue = 1.0;
	} 
	return combinedValue;
}

/**	Sender node: Prepare the message set to send to recipient node by updating information value using trust value.
 * 	\param sendSet stores the node's message set to send
 * 	\param msgtype
 *	\param receiver
 */
void Node::computeInfoUpdate(sourceMap &sendSet, int msgtype, int receiver)
{
    double t, alpha, tUpdated;
    alpha = getTrustValue(receiver);
//    cout << "InfoUpdate sender: " << getID() << "\treceiver: " << (neighbors.at(receiver))->getID() << "\ttrust: " << alpha << endl;
    
    if(msgtype == 1){	//is an action message
		valueMap::iterator iter = sourceSetValue.begin();
		while(iter!=sourceSetValue.end()){	
			t = iter->second;
			tUpdated = t*alpha;
			valueMap nSet;
			nSet.insert(nPair(getID(), tUpdated));
			sendSet.insert(sPair(iter->first,nSet));   //sender tUpdated(s)
			iter++;
		}

	} else {			//is an abort message
		valueMap::iterator iter = abortSetValue.begin();
		while(iter!=abortSetValue.end()){	
			t = iter->second;
			tUpdated = t*alpha;
			valueMap nSet;
			nSet.insert(nPair(getID(), tUpdated));
			sendSet.insert(sPair(iter->first,nSet));   //sender tUpdated(s)
			iter++;
		}
	}

	return;
}

/**	Recipient node: Source set is the union of the sender set with it's current set
 *	\param senderSet
 *	\param msgtype
 */
void Node::computeInfoUnionMax(sourceMap senderSet, int msgtype)
{
	double max=0.0;
    sourceMap::iterator iterS = senderSet.begin();
	
	if(msgtype==1){
		while(iterS != senderSet.end()){
			sourceMap::iterator iterR = sourceSet.find(iterS->first);
			if(iterR == sourceSet.end()){
				sourceSet.insert(sPair(iterS->first,iterS->second));
			} else {
				valueMap::iterator iterR2 = (iterR->second).find((iterS->second).begin()->first);
				if(iterR2 == iterR->second.end()){
					(iterR->second).insert(nPair((iterS->second).begin()->first,(iterS->second).begin()->second));
				} else {
 					max = (iterR->second).begin()->second;	//max of current info value for s
 					if((iterS->second).begin()->second > (iterR->second).begin()->second){
 						max = (iterS->second).begin()->second;
 					}
					iterR2->second = max;
 					max=0.0;							
				}
			}
			iterS++;
		}
    } else {
		while(iterS != senderSet.end()){
			sourceMap::iterator iterR = abortSet.find(iterS->first);
			if(iterR == abortSet.end()){
				abortSet.insert(sPair(iterS->first,iterS->second));
			} else {
				valueMap::iterator iterR2 = (iterR->second).find((iterS->second).begin()->first);
				if(iterR2 == iterR->second.end()){
					(iterR->second).insert(nPair((iterS->second).begin()->first,(iterS->second).begin()->second));
				} else {
 					max = (iterR->second).begin()->second;	//max of current info value for s
 					if((iterS->second).begin()->second > (iterR->second).begin()->second){
 						max = (iterS->second).begin()->second;
 					}
					iterR2->second = max;
 					max=0.0;							
				}
			}
			iterS++;
		}    
    }
}

/**	Recipient node: Compute new fusedInfoValue
 */
double Node::computeInfoFusion()
{
    sourceMap::iterator iter = sourceSet.begin();
    double t, prodT=1.0;
	double sum = 0.0;
	double max = 0.0;
	double newfusedInfoValue = 0.0;
    //go through the values
    while(iter!=sourceSet.end()){
        t = fusedSingleSource(iter->second,lambda1);
        
        valueMap::iterator iter2 = sourceSetValue.find(iter->first);
        if(iter2 == sourceSetValue.end()){
        	sourceSetValue.insert(nPair(iter->first, t));        
        } else {
        	iter2->second = t;
        }
        
        //Reliability
        prodT = prodT * (1.0-t);
	    //Sum of the values
    	sum = sum + t;
    	//Max of the values
		if(t>max){
			max = t;
		}
        iter++;
    }

    if(lambda2 > 1){
	    //Use Reliability
	    newfusedInfoValue = 1.0-prodT;
    } else {
    	//Use SUM or MAX
    	newfusedInfoValue = lambda2*sum+(1.0-lambda2)*max;
		if(newfusedInfoValue>1.0){							
			newfusedInfoValue = 1.0;
		}   
	}
	
	if(newfusedInfoValue > fusedInfoValue){
		hasUpdated = 1;
	} else {
		hasUpdated = 0;
	}
	fusedInfoValue = newfusedInfoValue;
//	cout << "InfoFusion NodeID: " << getID() << "\tfusedInfoValue: " << fusedInfoValue << endl;
    return fusedInfoValue;
}

/**	Recipient node: Compute new fusedAbortValue
 */
double Node::computeAbortFusion()
{
    sourceMap::iterator iter = abortSet.begin();
    double t, prodT=1.0;
	double sum = 0.0;
	double max = 0.0;
	double newfusedAbortValue = 0.0;
    //go through the values
    while(iter!=abortSet.end()){
        t = fusedSingleSource(iter->second,lambdaAbort1);

        valueMap::iterator iter2 = abortSetValue.find(iter->first);
        if(iter2 == abortSetValue.end()){
        	abortSetValue.insert(nPair(iter->first, t));        
        } else {
        	iter2->second = t;
        }
        
        //Reliability
        prodT = prodT * (1.0-t);
	    //Sum of the values
    	sum = sum + t;
    	//Max of the values
		if(t>max){
			max = t;
		}
        iter++;
    }

    if(lambdaAbort2 > 1){
	    //Use Reliability
	    newfusedAbortValue = 1.0-prodT;
    } else {
    	//Use SUM or MAX
    	newfusedAbortValue = lambdaAbort2*sum+(1.0-lambdaAbort2)*max;
		if(newfusedAbortValue>1.0){							
			newfusedAbortValue = 1.0;
		}   
	}
	
	if(newfusedAbortValue > fusedAbortValue){
		hasUpdated = 1;
	} else {
		hasUpdated = 0;
	}
	fusedAbortValue = newfusedAbortValue;
//	cout << "InfoFusion NodeID: " << getID() << "\tfusedInfoValue: " << fusedInfoValue << endl;
    return fusedAbortValue;
}

/**	Recipient node: Update status of the node
 *
 */
void Node::updateStatus()
{
    if(status == 3 || status == -3){    //do nothing if this is a source node
        return;
    }
    if(fusedInfoValue <= 0 && fusedAbortValue <= 0){
        status = -1;            //still NotInformed
    } else if (fusedInfoValue-fusedAbortValue <= LB){
    	if(status==2){
//    		cout << "Reset:" << timeRemoved << "\t" << timeReset << endl;
    		timeRemoved = timeReset;
    	}
        status = 0;             //Disbelieved
        if(fusedInfoValue-fusedAbortValue <= sigma){
        	spreadAbortflag = 1;
        }
    } else if (fusedInfoValue-fusedAbortValue < UB){
    	if(status==2){
//    		cout << "Reset:" << timeRemoved << "\t" << timeReset << endl;
    		timeRemoved = timeReset;
    	}
        status = 1;             //Undecided
    } else if (fusedInfoValue-fusedAbortValue >= UB){
        status = 2;             //Believed
    }	
}
