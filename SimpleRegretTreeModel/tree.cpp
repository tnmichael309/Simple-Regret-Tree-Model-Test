#include "tree.h"

bool simpleNodeData::sort(simpleNodeData& node1, simpleNodeData& node2){
	return node1.mean > node2.mean;
};


treeNode::treeNode(){
	parent = NULL;
	init();
};

treeNode::~treeNode(){};

void treeNode::init(){
	moveValue = 0;
	cumulativeMoveValue = 0;
	heuristic = 0;

	if(m_vAlgorithmData.size() == 0){
		m_vAlgorithmData.push_back(NULL);
		m_vAlgorithmData.push_back(NULL);
		m_vAlgorithmData.push_back(NULL);
	}else{
		for(int i = 0; i < (int)MAX_TREE_NODE_DATA_NUM; i++){
			if(m_vAlgorithmData[i]){
				delete m_vAlgorithmData[i];
				m_vAlgorithmData[i] = NULL;
			}
		}
	}
};

void treeNode::genMoveValue(bool isPositive){
	moveValue = rand()%128;
	moveValue = isPositive ? moveValue : -1*moveValue;
}

void treeNode::genSpecialMoveValue(bool isPositive, bool isIncreasingRewards, int depth){
	double startingReward = 128.0;
	int range = 0;
	if(isIncreasingRewards) range = (increaseRate == 1.0) ? startingReward : (int)(startingReward*pow(increaseRate,ceil(depth*1.0/2.0)));
	else range = startingReward;
	
	if(range < 2) range = 2;
	moveValue = rand()%range;
	moveValue = isPositive ? moveValue : -1*moveValue;
}


void treeNode::showInfo(){
	cerr << " move value: " << moveValue << " cumulative rewards: " << cumulativeMoveValue << " heuristics: " << heuristic << endl;
}

baseTreeNodeData* treeNode::getAlgorithmNodeData( TREENODEDATA treeNodeData )
{
	return (m_vAlgorithmData[(int)treeNodeData]);
}

void treeNode::updateAlgorithmNodeData( baseTreeNodeData* nodeData,TREENODEDATA treeNodeData )
{
	m_vAlgorithmData[(int)treeNodeData] = nodeData;
}


tree::tree(int bf, int depth, int expType){
	m_depth = depth;
	m_bf = bf;
	m_root = new treeNode();

	/*cout << sizeof(*m_root) << endl;
	cout << 16*sizeof(treeNode*)<< endl;
	cout << m_root->getAlgorithmNodeData(MCTS_DATA) << "\t" << sizeof(*(m_root->getAlgorithmNodeData(MCTS_DATA)))<< endl;
	cout << m_root->getAlgorithmNodeData(SIMUST_DATA) << "\t" << sizeof(*(m_root->getAlgorithmNodeData(SIMUST_DATA)))<< endl;
	cout << m_root->getAlgorithmNodeData(SHOT_DATA) << "\t" << sizeof(*(m_root->getAlgorithmNodeData(SHOT_DATA)))<< endl;
*/


	m_expType = expType;

	double max,min;
	max = 0.0;
	min = 1.0;
	buildTree(m_root,0,max,min);
}

tree::~tree(){
	deleteTree(m_root,0);
};

// reset tree information
void tree::resetTree(treeNode* startNode, int depth, double& max, double& min){

	if(startNode == NULL){
		startNode = m_root;
		startNode->init();
	}

	double tempMax,tempMin;
	for(int i = 0; i < m_bf; i++){

		treeNode* newNode = startNode->children[i];
		
		// reset tree node information
		newNode->init();

#if isUseNormal
		// normal random values
		if(depth % 2 == 0) newNode->genMoveValue(true);
		else newNode->genMoveValue(false);
#else
		//special random values
		// left -> right: 0.5 ... 1
		if(startNode == m_root) {
			if(i == 0) newNode->increaseRate = 0.5;
			else if(i > 0 && i < m_bf-1) newNode->increaseRate = 0.5 + i*0.5/m_bf;
			else newNode->increaseRate = 1.0;
		}else {
			newNode->increaseRate = startNode->increaseRate;
		}
		if(depth % 2 == 0) newNode->genSpecialMoveValue(true,true,depth);
		else newNode->genSpecialMoveValue(false,true,depth);
#endif

		newNode->cumulativeMoveValue = newNode->moveValue + startNode->cumulativeMoveValue;
		//newNode->heuristic = newNode->cumulativeMoveValue;

		tempMax = 0.0;
		tempMin = 1.0;
		if(depth < m_depth-1) resetTree(newNode, depth+1,tempMax, tempMin);
		else;

		if(depth == m_depth - 1) {

#if isUseNormal
			double maxReward = m_depth*128.0/2.0;
#else
			double maxReward = 0.0;
			double round = ceil(m_depth*1.0/2.0);

			maxReward = 128.0*round;
			/*if(newNode->increaseRate == 1){
				maxReward = 128.0*round;
			}else{
				maxReward = 128.0*(1.0-pow(newNode->increaseRate,round))/(1-newNode->increaseRate);
			}*/
#endif

#if isZeroOrOneHeuristic
			if(newNode->cumulativeMoveValue > 0) newNode->heuristic = 1;
			if(newNode->cumulativeMoveValue == 0) newNode->heuristic = 0.5;
			if(newNode->cumulativeMoveValue < 0) newNode->heuristic = 0;
#else
			newNode->heuristic = (newNode->cumulativeMoveValue+maxReward)/(2.0*maxReward);
			newNode->varianceHeuristic = minVariance;

			tempMax = pow(newNode->heuristic,2.0);
#endif
		}
		startNode->heuristic +=  newNode->heuristic;
		max += tempMax;
	}
	startNode->heuristic /= (double)m_bf;
	max /= (double)m_bf;
	double vheuristic = max-pow(startNode->heuristic,2.0);
	if(vheuristic <= minVariance) vheuristic = minVariance;
	startNode->varianceHeuristic = vheuristic;
	//startNode->varianceHeuristic = max-pow(startNode->heuristic,2.0);
}


// 
// root -> depth : 0
// min -> depth : 1
// depth: depth of the startNode
//
void tree::buildTree(treeNode* startNode, int depth, double& max, double& min){
	
	double tempMax,tempMin;

	for(int i = 0; i < m_bf; i++){
		treeNode* newNode = new treeNode();
		
#if isUseNormal
		// normal random values
		if(depth % 2 == 0) newNode->genMoveValue(true);
		else newNode->genMoveValue(false);
#else
		//special random values
		// left -> right: 0.5 ... 1
		if(startNode == m_root) {
			if(i == 0) newNode->increaseRate = 0.5;
			else if(i > 0 && i < m_bf-1) newNode->increaseRate = 0.5 + i*0.5/m_bf;
			else newNode->increaseRate = 1.0;
		}else {
			newNode->increaseRate = startNode->increaseRate;
		}
		if(depth % 2 == 0) newNode->genSpecialMoveValue(true,true,depth);
		else newNode->genSpecialMoveValue(false,true,depth);
#endif

		newNode->parent = startNode;
		startNode->children.push_back(newNode);
		newNode->cumulativeMoveValue = newNode->moveValue + startNode->cumulativeMoveValue;
		//newNode->heuristic = newNode->cumulativeMoveValue;

		tempMax = 0.0;
		tempMin = 1.0;
		if(depth < m_depth-1) buildTree(newNode, depth+1,tempMax,tempMin);
		else;

		if(depth == m_depth - 1) {

#if isUseNormal
			double maxReward = m_depth*128.0/2.0;
#else
			double maxReward = 0.0;
			double round = ceil(m_depth*1.0/2.0);

			maxReward = 128.0*round;
			/*if(newNode->increaseRate == 1){
				maxReward = 128.0*round;
			}else{
				maxReward = 128.0*(1.0-pow(newNode->increaseRate,round))/(1-newNode->increaseRate);
			}*/
#endif

#if isZeroOrOneHeuristic
			if(newNode->cumulativeMoveValue > 0) newNode->heuristic = 1;
			if(newNode->cumulativeMoveValue == 0) newNode->heuristic = 0.5;
			if(newNode->cumulativeMoveValue < 0) newNode->heuristic = 0;
#else
			newNode->heuristic = (newNode->cumulativeMoveValue+maxReward)/(2.0*maxReward);
			newNode->varianceHeuristic = minVariance;

			tempMax = pow(newNode->heuristic,2.0);

#endif
		}
		startNode->heuristic +=  newNode->heuristic;
		max += tempMax;
	}
	startNode->heuristic /= (double)m_bf;
	max /= (double)m_bf;
	double vheuristic = max-pow(startNode->heuristic,2.0);
	if(vheuristic <= minVariance) vheuristic = minVariance;
	startNode->varianceHeuristic = vheuristic;

};


void tree::dumpTree(){
	showTree(m_root,0);
}

void tree::deleteTree(treeNode* startNode, int depth){
	
	static int counter = 0;
	for(int i = 0; i < m_bf; i++){
		treeNode* nextNode = startNode->children[i];
		
		counter++;
		if(depth < m_depth-1) deleteTree(nextNode, depth+1);
		else;
		counter--;

		delete nextNode;
	}

	if(counter == 0) delete startNode;
}

int tree::getBestMove(bool isShowInfo){
	int maxMoveScore = 0;
	int bestBranch = 0;
	for(int i = 0; i < m_bf; i++){
		treeNode* nextNode = m_root->children[i];

		// run a job for 200 simulations
		/*runMCTSJob(nextNode,50);
		nextNode->heuristic = nextNode->getAlgorithmNodeData(MCTS_DATA)->getWinRate();*/


		int score = 0;
		if(0 < m_depth-1) score = traverseForBestMove(nextNode, 1);
		else score = nextNode->cumulativeMoveValue;

		if(i == 0 || score > maxMoveScore){
			maxMoveScore = score;
			bestBranch = i;
		}
	}

	double realScore = maxMoveScore;
	double maxReward = 128.0*m_depth/2.0;
	m_root->realScore = (realScore+maxReward)/(2.0*maxReward);

	//clearMCTSInfo(m_root,0);

	if(isShowInfo) cerr << "Real Best Move is: " << bestBranch << "\tScore: " <<  (double)(maxMoveScore+128.0*m_depth/2.0)/(128.0*m_depth) << endl;
	return bestBranch;
};

int tree::traverseForBestMove(treeNode* startNode, int depth){
	
	int maxMoveScore = 0;
	for(int i = 0; i < m_bf; i++){
		treeNode* nextNode = startNode->children[i];

		// run a job for 200 simulations
		/*runMCTSJob(nextNode,50);
		nextNode->heuristic = nextNode->getAlgorithmNodeData(MCTS_DATA)->getWinRate();*/

		int score = 0;
		if(depth < m_depth-1) score = traverseForBestMove(nextNode, depth+1);
		else score = nextNode->cumulativeMoveValue;
	
		if(depth % 2 == 0);
		else score *= -1;

		if(i == 0 || score > maxMoveScore){
			maxMoveScore = score;
		}

	}
	

	double realScore = depth%2==0 ? maxMoveScore : maxMoveScore*-1;
	double maxReward = 128.0*m_depth/2.0;
	startNode->realScore = (realScore+maxReward)/(2.0*maxReward);

	if(depth % 2 == 0) return maxMoveScore;
	else return maxMoveScore*-1;
};
