#include "mcts.h"

void mcts::runAlgorithm( int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove )
{
	mctsData* rootStoData = (mctsData*)m_root->getAlgorithmNodeData(MCTS_DATA);
	rootStoData = new mctsData;
	rootStoData->init();
	m_root->updateAlgorithmNodeData(rootStoData, rootStoData->getDataType());

	int currentRound = 0;
	while(currentRound < budget){
		treeNode* selectNode = mctsSelect();
		if(isUseHeuristic) mctsUpdate(selectNode, ((mctsData*)selectNode->getAlgorithmNodeData(MCTS_DATA))->getWinRate());
		else mctsUpdate(selectNode, (double)mctsSimulate(selectNode));

		currentRound++;
	}

	bestMove = mctsRecommed();
}

treeNode* mcts::mctsSelect( treeNode* node /*= NULL*/ )
{
	if(node == NULL) node = m_root;
	else;
	int depth = 0;
	if(node != NULL){
		treeNode* temp = node;
		while(temp->parent) {
			temp = temp->parent;
			depth++;
		}
	}

	while(((mctsData*)node->getAlgorithmNodeData(MCTS_DATA))->isExpanded == true){

		// leaf node: break
		if(depth == m_depth) break;

		double maxActionScore = 0;
		treeNode* nextNode = NULL;

		for(int i = 0; i < m_bf; i++){
			treeNode* child = node->children[i];
			mctsData* childNodeData = (mctsData*) child->getAlgorithmNodeData(MCTS_DATA);
			double score = 0;

			if(childNodeData->isExpanded) {

				// win rate
				score = childNodeData->getWinRate();
				if(depth % 2 == 0); // max node
				else score *= -1.0; // min node
				
				// add exploration term
				score += 1.0/sqrt(2.0)*sqrt(2.0*log((double)((mctsData*)node->getAlgorithmNodeData(MCTS_DATA))->visitCount)/childNodeData->visitCount);

			}else score = 100000000.0;

			if(score > maxActionScore || nextNode == NULL) {
				maxActionScore = score;
				nextNode = child;
			}
		}


		depth ++;
		node = nextNode;
	};

	treeNode* bestChild = NULL;
	if(depth == m_depth){
		//node->getAlgorithmNodeData(MCTS_DATA)->totalRewards += node->heuristic;
		bestChild = node;
	}
	else{
		double maxMeanOfUnexpandedNode = -100.0;
		
		for(int i = 0; i < m_bf; i++){

			treeNode* child = node->children[i];
			mctsData* childNodeData = (mctsData*) child->getAlgorithmNodeData(MCTS_DATA);

			double mean = child->heuristic;

			if(depth % 2 == 0);
			else mean = 1 - mean; //*= -1.0;

			// set up mcts data
			childNodeData = new mctsData;
			childNodeData->init();
			child->updateAlgorithmNodeData(childNodeData, childNodeData->getDataType());

			// expand the i'th child
			childNodeData->totalRewards = (double)child->heuristic;
			childNodeData->visitCount = 1;

			if(mean > maxMeanOfUnexpandedNode) {
				maxMeanOfUnexpandedNode = mean;
				bestChild = child;
			}
		}
		//node->getAlgorithmNodeData(MCTS_DATA)->totalRewards = bestChild->heuristic;
		//node->getAlgorithmNodeData(MCTS_DATA)->visitCount = 0;
	}

	((mctsData*)node->getAlgorithmNodeData(MCTS_DATA))->isExpanded = true;
	//node->getAlgorithmNodeData(MCTS_DATA)->visitCount ++;

	return bestChild;
}

double mcts::mctsSimulate( treeNode* node )
{
	while(node->children.size() != 0){
		node = node->children[rand()%m_bf];
	}
	return node->heuristic;
	//return (double)(node->cumulativeMoveValue > 0);
}

void mcts::mctsUpdate( treeNode* node, double reward, treeNode* endNode /*= NULL*/ )
{
	if(endNode == NULL) endNode = m_root;
	else;

	// find depth for 'node'
	treeNode* tempNode = node;
	int depth = 0;
	while(tempNode->parent){
		depth++;
		tempNode = tempNode->parent;
	};


	while(node){

		// uct update
		((mctsData*)node->getAlgorithmNodeData(MCTS_DATA))->visitCount += 1;
		((mctsData*)node->getAlgorithmNodeData(MCTS_DATA))->totalRewards += reward;

		node = node->parent;
		depth--;
	};

}

int mcts::mctsRecommed()
{
	int maxVisitCount = 0;
	int selectArm = 0;


	for(int i = 0; i < m_bf; i++){
		int visitCount = ((mctsData*)(m_root->children[i])->getAlgorithmNodeData(MCTS_DATA))->visitCount;
		if(visitCount > maxVisitCount){
			maxVisitCount = visitCount;
			selectArm = i;
		}
	}

	return selectArm;
}

void mcts::dumpAllocationInfoToFile( string fileName )
{
	
	fstream f;
	f.open(fileName.c_str(), ios::out | ios::app);

	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		s.mean = ((mctsData*)m_root->children[i]->getAlgorithmNodeData(MCTS_DATA))->visitCount; //();
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	for(int i = 0; i < vNodeData.size(); i++){
		int originalIndex = vNodeData[i].originalIndex;
		f << vNodeData[i].mean << " " << m_root->children[originalIndex]->realScore << " ";
		f <<  ((mctsData*)m_root->children[originalIndex]->getAlgorithmNodeData(MCTS_DATA))->visitCount;
		f << " ";
	}

	f << endl;
}

int mcts::getNthBestMoveOfAlgorithms( int n )
{
	
	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		s.mean = ((mctsData*)m_root->children[i]->getAlgorithmNodeData(MCTS_DATA))->visitCount; //getWinRate();
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	if(n < 0 || n >= m_bf) return -1;
	else return vNodeData[n].originalIndex;
}

void mcts::showTree( treeNode* startNode, int depth )
{
	if(startNode == m_root || startNode->parent == m_root);
	else return;


	for(int i = 0; i < m_bf; i++){
		treeNode* nextNode = startNode->children[i];
		
		string emptyString = "";
		for(int k = 0; k < depth+1; k++) emptyString += string("      ");
		if(depth % 2 == 0) cerr << emptyString << "min node:\n";
		else cerr << emptyString << "max node:\n";

		if(nextNode->getAlgorithmNodeData(MCTS_DATA)){
			//cerr << "=========== MCTS Data===============\n";
			cerr << emptyString << "UCT win rate: " << ((mctsData*)nextNode->getAlgorithmNodeData(MCTS_DATA))->getWinRate()
				<< " visit: " << ((mctsData*)nextNode->getAlgorithmNodeData(MCTS_DATA))->visitCount << "\n";
			//cerr << "=====================================\n";
		}
		
		
		/*cerr << "mean: " << nextNode->getAlgorithmNodeData(SHOT_DATA)->totalRewards/nextNode->getAlgorithmNodeData(SHOT_DATA)->visitCount
			<< " visit: " << nextNode->getAlgorithmNodeData(SHOT_DATA)->visitCount << " ";
*/


		//nextNode->showInfo();
		 
		if(depth < m_depth-1) showTree(nextNode, depth+1);
		else;
	}
}
