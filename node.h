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
using namespace std;

typedef map<int,double> valueMap;
typedef pair<int,double> nPair;
typedef map<int,valueMap> sourceMap;
typedef pair<int,valueMap> sPair;

class Node
{
    public:
		Node (int nid){
            id = nid;
            type = -1;
            group = -1;
            LB = -1;
            UB = -1;
            sigma = -1;
            fusedInfoValue = 0.0;
            fusedAbortValue = 0.0;
            status = -1;
            timeAbort = 999;	
            timeRemoved = 999;    
            timeReset = 999;
            timeActivate = 0;
			psend = 1.0;
			lambda1 = 0.0;
			lambda2 = 0.0;	
			lambdaAbort1 = 0.0;
			lambdaAbort2 = 0.0;				
			hasUpdated = 0;
			receivedflag = 0;
			receivedAbortflag = 0;
			spreadAbortflag = 0;			

		}
        
        //regular node
        Node (int nid, int ntype, int ngroup, double nLB, double nUB, double nsigma){
            id = nid;
            type = ntype;
            group = ngroup;
            LB = nLB;
            UB = nUB;
            sigma = nsigma;
            fusedInfoValue = 0.0;
            fusedAbortValue = 0.0;
            status = -1;
            timeAbort = 999;
            timeRemoved = 999;  
            timeReset = 999;  
            timeActivate = 0;
			psend = 1.0;
			lambda1 = 0.0;
			lambda2 = 0.0;	
			lambdaAbort1 = 0.0;
			lambdaAbort2 = 0.0;				
			hasUpdated = 0;
			receivedflag = 0;
			receivedAbortflag = 0;
			spreadAbortflag = 0;
        }

		//source node
        Node (int sid, int stype, int sgroup, double svalue, double pval, int sremove, int sactive){
            id = sid;
            type = stype;
            group = sgroup;
            LB = -1;
            UB = -1;
            sigma = -1;
            status = 3;
            timeAbort = 999;
            timeRemoved = sremove;  
            timeReset = 999;  
            timeActivate = sactive;
			psend = pval;
			lambda1 = 0.0;
			lambda2 = 0.0;	
			lambdaAbort1 = 0.0;
			lambdaAbort2 = 0.0;				
			hasUpdated = 0;
			receivedflag = 0;
			receivedAbortflag = 0;
			spreadAbortflag = 0;
			
			if(sgroup == 1 ){
				sourceSetValue.insert(nPair(sid,svalue));
				fusedInfoValue = svalue; fusedAbortValue = 0;
			} else if(sgroup == 2){
				abortSetValue.insert(nPair(sid,svalue));
				fusedInfoValue = 0; fusedAbortValue = svalue;
			} else {}
        }

        ~Node(){};
        int getID(){ return id;}
        int getStatus(){ return status;}
        int getType(){ return type;}
        int getGroup(){ return group;}
        double getValue() {return fusedInfoValue;}
        double getAbortValue() {return fusedAbortValue;}
		double getLB() {return LB;}
		double getUB() {return UB;}

		void setNode(int ntype, double nLB, double nUB, double nsigma){
			type = ntype; LB = nLB; UB = nUB; sigma = nsigma;};	//set node member values
        void setStatus(int nstatus){ status = nstatus;};		//set node state

		void setRemovedTime(int etime);		//set time between believed state and removed state
		void setAbortTime(int stime);		//set time for spreading in disbelieved state

		void setpSend(double p);		//set probability of successfully communicating with others

        int getActivate(){ return timeActivate;}	//time that source node sends initial message
        void setActivate(int sactive){ timeActivate = sactive;};
		
		double getTrustValue(int receiver);			//returns the trust weight for neighbor
		
		//Lambda parameters used in computeInfoUnionMax and computeInfoFusion
		void setLambdaParam1(double lp){lambda1 = lp;};
		void setLambdaParam2(double lp){lambda2 = lp;};
		
		//Lambda parameters used in computeInfoUnionMax and computeAbortFusion
		void setLambdaAbortParam1(double lq){lambdaAbort1 = lq;};
		void setLambdaAbortParam2(double lq){lambdaAbort2 = lq;};

        void removeInEdges();                       //remove all incoming edges
        void removeOutEdges();                      //remove all outgoing edges
        
        void sendMessage(int msgtype);                 		//send message to neighbors
		void sendMessageToNode(int msgtype, int nodeID);	//send message to specifed node
        
        //receive the message and process it
        void receiveMessage(sourceMap senderSet, int msgtype);
		//query neighbors
		void queryNeighbors();						

		//fuse single source value
		double fusedSingleSource(valueMap infoSet, double lambda);
		
		//update information value
        void computeInfoUpdate(sourceMap &sendSet, int msgtype, int receiver);
        
        //compute information union
        void computeInfoUnionMax(sourceMap senderSet, int msgtype);
        
        double computeInfoFusion();		//compute fusedInfoValue for this node
        double computeAbortFusion();    //compute fusedAbortValue for this node
        void updateStatus();

        int timeRemoved;				//time until node is removed
        int timeReset;					//time until node is removed (saved value)
        int timeActivate;				//time that source node sends initial message
		int timeAbort;					//time until node stops spreading Abort message
        									
		double psend;					//node's probability of successfully communicating with others

		double lambda1;					//Lambda parameter used in InfoUnionMax function
		double lambda2;					//Lambda parameter used in InfoFusion function
		
		double lambdaAbort1;			//parameter for AbortUnionMax function
		double lambdaAbort2;			//parameter for AbortFusion function
		
        sourceMap sourceSet;            //stores the action message (sourceID,value) pair
        valueMap sourceSetValue;		
        sourceMap abortSet;				//stores the abort message (sourceID,value) pair
        valueMap abortSetValue;			

        vector<Node *> neighbors;			//neighbor nodes
        vector<double> neighbors_weights; 	//neighbor nodes outgoing trust weights
        vector<int> received_from;			//received information from neighbor
        
        int hasUpdated;					/**< fusedinfovalue has changed */
        int receivedflag;				/**< received an Action message */
        int receivedAbortflag;			/**< received an Abort message */
        int spreadAbortflag;			/**< in Disbelieved state and spreads abort information */

    private:
        int id;							/**< node internal id */
        int status;						/**< node state */
        int type;						/**< node type (used for assigning node characteristics) */
        int group;						/**< node group (used for assigning trust values) */
        double LB;						/**< lower bound threshold */
        double UB;						/**< upper bound threshold */
        double sigma;					/**< threshold for spreading abort */
        double fusedInfoValue;			/**< fused information value of action information */
        double fusedAbortValue;			/**< fused information value of abort information */

};

