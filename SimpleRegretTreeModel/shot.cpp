#include "shot.h"

void shot::runAlgorithm( int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove )
{
	shotData* rootStoData = (shotData*)m_root->getAlgorithmNodeData(SHOT_DATA);

	rootStoData = new shotData;
	rootStoData->init();
	m_root->updateAlgorithmNodeData(rootStoData, rootStoData->getDataType());

	rootStoData->totalRewards = m_root->heuristic;
	rootStoData->visitCount = 1;

	int b = 0,p = 0;
	double r = 0;
	bestMove = doSHOT(m_root,budget, b,p,r, isUseHeuristic);
}

double shot::allocateABudget( treeNode* startNode, bool isUseHeuristic )
{
	// find depth for 'node'
	treeNode* tempNode = startNode;
	int depth = 0;
	while(tempNode->parent){
		depth++;
		tempNode = tempNode->parent;
	};

	double addedReward = 0.0;
	shotData* startNodeData = (shotData*)startNode->getAlgorithmNodeData(SHOT_DATA);

	if(startNodeData == NULL){
		startNodeData = new shotData;
		startNodeData->init();
		startNode->updateAlgorithmNodeData(startNodeData, startNodeData->getDataType());

		startNodeData->totalRewards = startNode->heuristic;
		startNodeData->visitCount = 1;
	}


	if(depth == m_depth){
		addedReward = startNode->heuristic;
	}
	else{
		double maxMeanOfUnexpandedNode = -100.0;
		treeNode* bestChild = NULL;

		for(int i = 0; i < m_bf; i++){

			treeNode* child = startNode->children[i];
			shotData* childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA);


			if(childNodeData == NULL){
				// init shot data
				childNodeData = new shotData;
				childNodeData->init();
				child->updateAlgorithmNodeData(childNodeData, childNodeData->getDataType());

				// expand the i'th child
				childNodeData->totalRewards = child->heuristic;
				childNodeData->visitCount = 1;
			}
				

			double mean = childNodeData->totalRewards/childNodeData->visitCount;
			if(depth % 2 == 0);
			else mean = 1 - mean; //*= -1.0;

			if(mean > maxMeanOfUnexpandedNode) {
				maxMeanOfUnexpandedNode = mean;
				bestChild = child;
			}
		}

		if(isUseHeuristic){
			addedReward = maxMeanOfUnexpandedNode;
			if(depth % 2 == 0);
			else addedReward = (1.0 - addedReward);
		}
		else addedReward = shotSimulate(bestChild);

	}
	
	startNodeData->visitCount += 1;
	startNodeData->totalRewards += addedReward;

	return addedReward;
}

int shot::doSHOT( treeNode* startNode, int budget, int& budgetUsed, int& playouts, double& rewards, bool isUseHeuristic )
{
	// find depth for 'node'
	treeNode* tempNode = startNode;
	int depth = 0;
	while(tempNode->parent){
		depth++;
		tempNode = tempNode->parent;
	};

	if(startNode->children.size() == 0 || budget == 1){ // leaf node or only one budget left
		
		while(budgetUsed < budget){
			budgetUsed+=1;
			playouts+=1;

			rewards += allocateABudget(startNode,isUseHeuristic);
			
		}

		return -1;
	}
	shotData* startNodeData = (shotData*)startNode->getAlgorithmNodeData(SHOT_DATA);


	if( startNodeData->visitCount <= m_bf){ // for those never played branches, play once
		for(int i = 0; i < m_bf; i++){
			treeNode* child = startNode->children[i];
			shotData* childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA);

			if(childNodeData == NULL){
				budgetUsed++;
				playouts++;

				double addedReward = allocateABudget(child,isUseHeuristic);
				childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA); // childNodeData will be reallocated after allocating a budget

				rewards+= addedReward;

				startNodeData->visitCount++;
				startNodeData->totalRewards += addedReward;

				if(playouts == budget) return -1;
			}
		}
	}

	vector<simpleNodeData> vSimpleNodes;
	for(int i = 0; i < m_bf; i++){
		treeNode* child = startNode->children[i];
		shotData* childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA);

		simpleNodeData sm;
		sm.mean = depth % 2 == 0 ? childNodeData->totalRewards/childNodeData->visitCount : 1-childNodeData->totalRewards/childNodeData->visitCount;
		sm.originalIndex = i;
		vSimpleNodes.push_back(sm);
	}
	sort(vSimpleNodes.begin(),vSimpleNodes.end(),simpleNodeData::sort);


	int b = 0;
	while(vSimpleNodes.size() > 1){
		int possibleBudget = floor((startNodeData->visitCount+budget-budgetUsed)/(vSimpleNodes.size() * ceil(log((double)m_bf)/log(2.0))));
		b += (1 > possibleBudget ? 1 : possibleBudget);
		
		/*if(startNode == m_root || startNode->parent == m_root){

			if(startNode == m_root) cerr << "In root" << endl;
			else cerr << "In second layer" << endl;

			cerr << "parent visit count: " << startNode->getAlgorithmNodeData(SHOT_DATA)->visitCount << endl
				<< "total budget allowed: " << budget << endl
				<< "current rest children number: " << vSimpleNodes.size() << endl
				<< "total number of rounds: " << ceil(log((double)m_bf)/log(2.0))<< endl;

			cerr << "rest budget: " << budget - budgetUsed << endl;
			cerr << "budget to be reached for each child: " << b << endl;
		}*/
		for(int i = 0; i < vSimpleNodes.size(); i++){

			treeNode* child = startNode->children[vSimpleNodes[i].originalIndex];
			shotData* childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA);

			if(childNodeData->visitCount < b){

				int b1 = b - childNodeData->visitCount; // extra budget to play
				
				if(startNode == m_root && vSimpleNodes.size() == 2 && i == 0){ // basically not run in our test cases
					treeNode* secondChild = startNode->children[vSimpleNodes[1].originalIndex];
					shotData* secondChildNodeData = (shotData*)secondChild->getAlgorithmNodeData(SHOT_DATA);
					b1 = budget - budgetUsed - (b - secondChildNodeData->visitCount);
				}

				b1 = b1 < (budget - budgetUsed) ? b1 : (budget - budgetUsed);

				int budgetUsedByChild = 0;
				int childPlayouts = 0;
				double childRewards = 0;
				doSHOT(child, b1, budgetUsedByChild, childPlayouts, childRewards,isUseHeuristic);
				// updat child node data
				childNodeData = (shotData*)child->getAlgorithmNodeData(SHOT_DATA);
				startNodeData = (shotData*)startNode->getAlgorithmNodeData(SHOT_DATA);

				// update playouts, budgetUsed and rewards
				playouts += childPlayouts;
				budgetUsed += budgetUsedByChild;
				rewards += childRewards;

				
				startNodeData->visitCount += childPlayouts;
				startNodeData->totalRewards += childRewards;

				// update mean data
				vSimpleNodes[i].mean = depth % 2== 0 ? childNodeData->getWinRate(): 1-childNodeData->getWinRate();
			}

			//if(startNode == m_root) cout << "budget used: " << budgetUsed << " budgets allocated: " << budget << " will break: " << (budgetUsed >= budget) << endl;


			if(budgetUsed >= budget) break;
		}

		sort(vSimpleNodes.begin(),vSimpleNodes.end(),simpleNodeData::sort);
		/*if(startNode == m_root || startNode->parent == m_root){

			if(startNode == m_root) cerr << "In root" << endl;
			else cerr << "In second layer" << endl;

			for(int i = 0; i < vSimpleNodes.size(); i++){
				treeNode* child = startNode->children[vSimpleNodes[i].originalIndex];
				cout << "child: " << vSimpleNodes[i].originalIndex 
					<< "\tvisit: " << child->getAlgorithmNodeData(SHOT_DATA)->visitCount 
					<< "\ttotal rewards: " << child->getAlgorithmNodeData(SHOT_DATA)->totalRewards
					<< "\tmean: " << child->getAlgorithmNodeData(SHOT_DATA)->totalRewards/child->getAlgorithmNodeData(SHOT_DATA)->visitCount << endl;
			}
			cerr << "child size: " << vSimpleNodes.size() << endl;
			getBestMove(true);
			dumpTree();

			cin.get();
		}*/
		/*if(startNode == m_root) {
			cerr << "child size: " << vSimpleNodes.size() << endl;
			cin.get();
			dumpTree();
			cin.get();
		}*/

		vector<simpleNodeData> newSimpleNodes;
		int newSize = ceil(vSimpleNodes.size()*1.0/2.0);
		newSimpleNodes.assign(vSimpleNodes.begin(), vSimpleNodes.begin()+newSize);
		vSimpleNodes = newSimpleNodes;

		if(budgetUsed >= budget) break;
	}
	
	//cerr << "SHOT recommended: " << vSimpleNodes[0].originalIndex << endl;
	return vSimpleNodes[0].originalIndex;
}

void shot::dumpAllocationInfoToFile( string fileName )
{
	
	fstream f;
	f.open(fileName.c_str(), ios::out | ios::app);

	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		s.mean = ((shotData*)m_root->children[i]->getAlgorithmNodeData(SHOT_DATA))->visitCount; //getWinRate();
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	for(int i = 0; i < vNodeData.size(); i++){
		int originalIndex = vNodeData[i].originalIndex;
		f << vNodeData[i].mean << " " << m_root->children[originalIndex]->realScore << " ";
		f <<  ((shotData*)m_root->children[originalIndex]->getAlgorithmNodeData(SHOT_DATA))->visitCount;
		f << " ";
	}

	f << endl;
}

int shot::getNthBestMoveOfAlgorithms( int n )
{
	
	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		s.mean = ((shotData*)m_root->children[i]->getAlgorithmNodeData(SHOT_DATA))->visitCount; //getWinRate();
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	if(n < 0 || n >= m_bf) return -1;
	else return vNodeData[n].originalIndex;
}

void shot::showTree( treeNode* startNode, int depth )
{

}

double shot::shotSimulate( treeNode* node )
{
	while(node->children.size() != 0){
		node = node->children[rand()%m_bf];
	}
	return node->heuristic;
	//return (double)(node->cumulativeMoveValue > 0);
}
