#include "SIMuST.h"

void SIMuST::runAlgorithm( int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove )
{
	// init root node
	stoData* rootStoData = (stoData*)m_root->getAlgorithmNodeData(SIMUST_DATA);
	rootStoData = new stoData;
	m_root->updateAlgorithmNodeData(rootStoData, rootStoData->getDataType());

	rootStoData->init();
	rootStoData->mean = m_root->heuristic;
	rootStoData->meanOfMean = m_root->heuristic;
	rootStoData->meanOfSquareSum = m_root->heuristic*m_root->heuristic;
	rootStoData->visitCount = 1;
	rootStoData->bestProbability = 1.0/(double)m_bf;
	rootStoData->bestProbabilityMean = rootStoData->mean;
	rootStoData->entropy = 1.0;
	rootStoData->sensitivity = sqrt(1.0/(double)m_bf);
	rootStoData->variance = 10.0;


	int currentRound = 0;
	while(currentRound < budget){
		treeNode* selectNode = stoSelect();
		stoUpdate(selectNode);
		currentRound += 1;
	}

	bestMove = stoRecommed();
}


treeNode* SIMuST::stoSelect()
{
	treeNode* node = m_root;
	int depth = 0;

	static int counter = 0;
	counter++;

	bool isParentSelectBest = true;

	while(((stoData*)node->getAlgorithmNodeData(SIMUST_DATA))->isExpanded == true){

		// leaf node: break
		if(depth == m_depth) {
			//cerr << "select leaf" << endl;
			break;
		}


		int maxMeanIndex = 0;
		double maxMean = 0;
		vector<double> vMeans; 
		vector<int> vMoveCount;
		vector<double> vVariance;
		vector<double> vMeanOfMeans; 
		double totalVar = 0.0;
		double maxVar= 0.0;
		double minVar= 100.0;
		for(int i = 0; i < m_bf; i++){

			treeNode* child = node->children[i];
			stoData* childStoData = (stoData*)child->getAlgorithmNodeData(SIMUST_DATA);

			// find child with max mean 
			double mean = childStoData->mean;
			if(depth % 2 == 0);
			else mean = 1 - mean; //*= -1.0;
			if(i == 0 || mean > maxMean){
				maxMean = mean;
				maxMeanIndex = i;
			}

			vMeans.push_back(mean);
			vMoveCount.push_back(childStoData->visitCount);


			double inputVariance = childStoData->variance / (double)childStoData->visitCount;
			if(inputVariance <= minVariance){
				vVariance.push_back(minVariance);
				minVar = minVariance;
				totalVar += minVariance;
			}else {
				vVariance.push_back(inputVariance);
				totalVar += inputVariance;
			}
			if(inputVariance > maxVar) maxVar = inputVariance;
			if(inputVariance < minVar) minVar = inputVariance;

			vMeanOfMeans.push_back(childStoData->meanOfMean);
		}


		double dMaxZValue = 0.0;
		double dMinZValue = 10000.0;
		for(int i = 0; i < m_bf; i++){

			if(i == maxMeanIndex) continue;

			double dSquareSumOfTwoDeviations = vVariance[maxMeanIndex] + vVariance[i];
			double dDifferenceOfWinRate = vMeans[maxMeanIndex]-vMeans[i];
			double dZValue = dDifferenceOfWinRate/sqrt(dSquareSumOfTwoDeviations);

			if(dZValue > dMaxZValue) dMaxZValue = dZValue;
			if(dZValue < dMinZValue) dMinZValue = dZValue;
		}
		if(dMaxZValue == dMinZValue) {
			dMaxZValue = 1.0;
			dMinZValue = 0.0;
		}

		/*double maxDev = sqrt(maxVar);
		double minDev = sqrt(minVar);
		for(int i = 0; i < vVariance.size(); i++){

			double dev = sqrt(vVariance[i]);

			dev = (dev - minVar + 1)/(maxVar - minVar + sqrt(2.0));
			vVariance[i] = pow(dev,2.0);
		}
		cerr << maxVar << "\t" << minVar << endl;
		cin.get();*/


#if isNormalizeVariance == 1
		double factor = 0.25/maxVar;
		for(int i = 0; i < vVariance.size(); i++){
			vVariance[i] = vVariance[i]*factor;
			//vVariance[i] = vVariance[i]/(totalVar*2.0);
		}
#endif


		treeNode* nextNode = NULL;

		double dLargestGradient = -100.0;
		double dGradientOfBestChild = 0.0;
		double sumOfSensitivity = 0;
		vector<double> vSensitivities;
		for(int i = 0; i < vMeans.size(); i++){


			if(i == maxMeanIndex) {
				vSensitivities.push_back(0.0);
				continue;
			}

			double dSquareSumOfTwoDeviations = vVariance[maxMeanIndex] + vVariance[i];
			double dDifferenceOfWinRate = vMeans[maxMeanIndex]-vMeans[i];
			double dZValue = dDifferenceOfWinRate/sqrt(dSquareSumOfTwoDeviations);
			
			// scale the zValue to 0 ~ 1, where z value has better linearity
			//dZValue = (dZValue - dMinZValue)/(dMaxZValue-dMinZValue);

			int n1 = vMoveCount[maxMeanIndex];
			int n2 = vMoveCount[i];
			double degreeOfFreedom = ((n1<=1) || (n2<=1)) ? 1 : pow(vVariance[maxMeanIndex]+vVariance[i],2.0)/(pow(vVariance[maxMeanIndex],2.0)/(n1-1)+pow(vVariance[i],2.0)/(n2-1));
			degreeOfFreedom = degreeOfFreedom <= 1 ? 1 : degreeOfFreedom;


#if isUsePropagate == 0
			double dGradient = calculateStudentTPDF(dZValue,degreeOfFreedom)/calculateStudentTCDF(dZValue,degreeOfFreedom)/dSquareSumOfTwoDeviations
				*(-1*sqrt(dSquareSumOfTwoDeviations) - 1/sqrt(dSquareSumOfTwoDeviations)*(vMeans[i] - vMeanOfMeans[i])*dDifferenceOfWinRate/(vMoveCount[i]*(vMoveCount[i]+1)));
#else
			
			// propagate
			double dGradient = calculateStudentTPDF(dZValue,degreeOfFreedom)/calculateStudentTCDF(dZValue,degreeOfFreedom)/dSquareSumOfTwoDeviations
				*(-1*sqrt(dSquareSumOfTwoDeviations));
#endif


			//node->children[i]->getAlgorithmNodeData(SIMUST_DATA)->sensitivity = fabs(dGradient);
			double inputGradient = fabs(dGradient);
			//if(inputGradient<=minVariance) inputGradient = minVariance;
			sumOfSensitivity += inputGradient;
			vSensitivities.push_back(inputGradient);
			dGradient = sqrt(vVariance[i])*inputGradient;



			if(dGradient > dLargestGradient){
				nextNode = node->children[i];
				dLargestGradient = dGradient;
			}


#if isUsePropagate == 0
			dGradientOfBestChild += calculateStudentTPDF(dZValue,degreeOfFreedom)/calculateStudentTCDF(dZValue,degreeOfFreedom)/dSquareSumOfTwoDeviations
				*(sqrt(dSquareSumOfTwoDeviations) 
				- 1/sqrt(dSquareSumOfTwoDeviations)*(vMeans[maxMeanIndex] - vMeanOfMeans[maxMeanIndex])*dDifferenceOfWinRate/(vMoveCount[maxMeanIndex]*(vMoveCount[maxMeanIndex]+1)));
#else
			
			// propagate
			dGradientOfBestChild += calculateStudentTPDF(dZValue,degreeOfFreedom)/calculateStudentTCDF(dZValue,degreeOfFreedom)/dSquareSumOfTwoDeviations
				*(sqrt(dSquareSumOfTwoDeviations));
			
#endif


		}

		double inputGradient = fabs(dGradientOfBestChild);
		sumOfSensitivity += inputGradient;
		vSensitivities[maxMeanIndex] = inputGradient;
		dGradientOfBestChild = sqrt(vVariance[maxMeanIndex])*inputGradient;
/*

		// normalize sensitivity to [0,1]
		for(int i = 0; i < vSensitivities.size(); i++){
			//if(vSensitivities[i]<0) cout << (i==maxMeanIndex) << endl;
			node->children[i]->getAlgorithmNodeData(SIMUST_DATA)->sensitivity = vSensitivities[i] / sumOfSensitivity;

			if(sumOfSensitivity == 0){
				if(i != maxMeanIndex) node->children[i]->getAlgorithmNodeData(SIMUST_DATA)->sensitivity = 0;
				else node->children[i]->getAlgorithmNodeData(SIMUST_DATA)->sensitivity = 1.0;
			}

		}
*/

		
		if(dGradientOfBestChild > dLargestGradient){
			nextNode = node->children[maxMeanIndex];
			dLargestGradient = dGradientOfBestChild;

			isParentSelectBest = true;
			/*if(depth == 0) {
				treeNode* temp = node->children[maxMeanIndex]->parent;
				cerr << "select best mean child" << endl;
			}*/
		}else isParentSelectBest = false;

		if(dLargestGradient == 0){
			int id = rand()%m_bf;
			nextNode = node->children[id];

			if(id == maxMeanIndex) isParentSelectBest = true;
			else isParentSelectBest = false;
		}
		depth ++;
		node = nextNode;
	};

	
	treeNode* bestChild = NULL;


	if(depth == m_depth) bestChild = node;
	else{
		double maxMeanOfUnexpandedNode = -100.0;

		double maxVariance = 0.0;
		for(int i = 0; i < m_bf; i++){

			treeNode* child = node->children[i];
			stoData* childStoData = (stoData*)child->getAlgorithmNodeData(SIMUST_DATA);

			double mean = child->heuristic;
			if(depth % 2 == 0);
			else mean = 1 - mean; //*= -1.0;

			// init sto data
			childStoData = new stoData;
			childStoData->init();
			child->updateAlgorithmNodeData(childStoData, childStoData->getDataType());

			// expand the i'th child
			childStoData->mean =(double)child->heuristic;
#if isUseMaxOfflineVariance == 1
			childStoData->variance = 10;
#else
			childStoData->variance = fabs(child->realScore-child->heuristic);//child->varianceHeuristic; //10.0;
#endif
			childStoData->meanOfSquareSum =childStoData->mean*childStoData->mean;
			childStoData->meanOfMean = childStoData->mean;
			childStoData->visitCount = 1;
			childStoData->bestProbability = 1.0/(double)m_bf;
			childStoData->bestProbabilityMean = childStoData->mean;
			childStoData->entropy = 1.0;
			childStoData->sensitivity = sqrt(1.0/(double)m_bf);
			childStoData->isBest = false;

			if(mean > maxMeanOfUnexpandedNode) {
				maxMeanOfUnexpandedNode = mean;
				bestChild = child;
			}

			if(childStoData->variance > maxVariance){
				maxVariance = childStoData->variance;
			}
		}

#if isUseMaxOfflineVariance == 0
		for(int i = 0; i < m_bf; i++){
			treeNode* child = node->children[i];
			stoData* childStoData = (stoData*)child->getAlgorithmNodeData(SIMUST_DATA);
			childStoData->variance = childStoData->variance*10.0/maxVariance;
		}
#endif

		((stoData*)bestChild->getAlgorithmNodeData(SIMUST_DATA))->isBest = true;

		//前一版沒助解掉
		//node->getAlgorithmNodeData(SIMUST_DATA)->mean = bestChild->heuristic;
	}

	((stoData*)node->getAlgorithmNodeData(SIMUST_DATA))->isExpanded = true;

	if(node == m_root){
		node = node;
	}

	//node->getAlgorithmNodeData(SIMUST_DATA)->visitCount ++;
	return bestChild;
}

void SIMuST::stoUpdate( treeNode* node )
{
	// find depth for 'node'
	treeNode* tempNode = node;
	int depth = 0;
	while(tempNode->parent){
		depth++;
		tempNode = tempNode->parent;
	};


	while(node->parent){

		treeNode* parent = node->parent;
		
		// find max mean value
		double maxMean = 0;
		vector<double> vBestProbabilityMeans; 
		vector<int> vMoveCount;
		vector<double> vVariance;

		treeNode* bestChild = NULL;
		treeNode* originBestChild = NULL;
		double childPropagateVariance = 0;
		double totalVar = 0.0;
		double maxVar = 0.0;
		for(int i = 0; i < m_bf; i++){
			treeNode* child = parent->children[i];
			stoData* childStoData = (stoData*)child->getAlgorithmNodeData(SIMUST_DATA);
			
			// from only expanded children
				
			double mean = childStoData->mean;
			if(depth%2 == 0) mean = 1 - mean; //*= -1.0; // child is 'max' node, parent is 'min' node
			else;

			childPropagateVariance += pow(childStoData->sensitivity,2.0)*childStoData->variance;

			if(i == 0 || mean > maxMean){
				maxMean = mean;
				bestChild = child;
			}

			if(childStoData->isBest){
				originBestChild = child;
			}

			if(depth%2 == 0) vBestProbabilityMeans.push_back(1-childStoData->bestProbabilityMean); // child is 'max' node, parent is 'min' node
			else vBestProbabilityMeans.push_back(childStoData->bestProbabilityMean);

			vMoveCount.push_back(childStoData->visitCount);

			//if(child->getAlgorithmNodeData(SIMUST_DATA)->variance != 0) vVariance.push_back(child->getAlgorithmNodeData(SIMUST_DATA)->variance / (double)child->getAlgorithmNodeData(SIMUST_DATA)->visitCount);
			//else vVariance.push_back(0.00000000001 / (double)child->getAlgorithmNodeData(SIMUST_DATA)->visitCount);

			double inputVariance =childStoData->variance / childStoData->visitCount;
			if(inputVariance <= minVariance) vVariance.push_back(minVariance);
			else vVariance.push_back(inputVariance);

			totalVar += inputVariance;

			if(inputVariance > maxVar) maxVar = inputVariance;

			
		}

#if isNormalizeVariance == 1
		/*double factor = 0.25/maxVar;
		for(int i = 0; i < vVariance.size(); i++){
			vVariance[i] = vVariance[i]*factor;
			//vVariance[i] = vVariance[i]/(totalVar*2.0);
		}*/
#endif

#if isResetVariance == 1
		// reset variance if from best to not the best
		if(bestChild != originBestChild){
			/*double difference = fabs(bestChild->getAlgorithmNodeData(SIMUST_DATA)->mean - originBestChild->getAlgorithmNodeData(SIMUST_DATA)->mean);
			if(difference == 0){
				
			}else{
				bestChild->getAlgorithmNodeData(SIMUST_DATA)->variance = 10;//pow(difference*100.0,2.0);
				originBestChild->getAlgorithmNodeData(SIMUST_DATA)->variance = 10;//pow(difference*100.0,2.0);
				/ *bestChild->getAlgorithmNodeData(SIMUST_DATA)->visitCount = 1;
				originBestChild->getAlgorithmNodeData(SIMUST_DATA)->visitCount = 1;* /
			}*/
			bestChild->getAlgorithmNodeData(SIMUST_DATA)->variance = 10.0;
			originBestChild->getAlgorithmNodeData(SIMUST_DATA)->variance = 10.0;
			/*bestChild->getAlgorithmNodeData(SIMUST_DATA)->visitCount = 1;
			originBestChild->getAlgorithmNodeData(SIMUST_DATA)->visitCount = 1;*/

			bestChild->getAlgorithmNodeData(SIMUST_DATA)->isBest = true;
			originBestChild->getAlgorithmNodeData(SIMUST_DATA)->isBest = false;
		}
#endif

		double maxBestProbability= 0;
		double maxBestProbabilityMean = 0;
		double newEntropy = 0;

		for(int i = 0; i < m_bf; i++){
			treeNode* child = parent->children[i];
			stoData* childStoData = (stoData*)child->getAlgorithmNodeData(SIMUST_DATA);

			// from only expanded children

			double bestProbability = calculateStudentTBestProbability(i, vBestProbabilityMeans, vVariance, vMoveCount);
			childStoData->bestProbability = bestProbability;

			newEntropy += (bestProbability == 1.0 || bestProbability == 0.0) ? 0 : bestProbability*log(bestProbability)/log((double)m_bf);

			/*if(bestProbability == 1.0 || bestProbability == 0.0) {
				parent = parent;
			}*/

			//cerr << i << " " << vBestProbabilityMeans[i] << endl;
			if(i == 0 || bestProbability > maxBestProbability){
				maxBestProbability = bestProbability;
				maxBestProbabilityMean = vBestProbabilityMeans[i];
			}
			
		}
		newEntropy *= -1;


		// update parent data
		stoData* parentStoData = (stoData*)parent->getAlgorithmNodeData(SIMUST_DATA);
		double lastMean = parentStoData->mean;
		parentStoData->mean = depth%2 == 0 ? 1.0-maxMean : maxMean;
		//parent->getAlgorithmNodeData(SIMUST_DATA)->meanOfMean = (1-parent->getAlgorithmNodeData(SIMUST_DATA)->entropy)*parent->getAlgorithmNodeData(SIMUST_DATA)->meanOfMean + parent->getAlgorithmNodeData(SIMUST_DATA)->entropy*parent->getAlgorithmNodeData(SIMUST_DATA)->mean;

		parentStoData->meanOfMean = (parentStoData->meanOfMean * parentStoData->visitCount + parentStoData->mean) 
			/ (parentStoData->visitCount + 1);


		// new variance: average squared difference bt the newest mean


#if isUsePropagate == 0
		// original variance
		if(parentStoData->visitCount >= 2){
			parentStoData->variance = parentStoData->variance * (parentStoData->visitCount - 1) / (parentStoData->visitCount)
				+ pow(parentStoData->mean - parentStoData->meanOfMean, 2.0)/(parentStoData->visitCount + 1);
		}
#else

		// propagated variance
		if(parentStoData->visitCount >= 2){
			parentStoData->variance = childPropagateVariance;
		}
#endif
		if(parentStoData->variance <= minVariance) parentStoData->variance = minVariance;


		parentStoData->visitCount += 1;
		parentStoData->bestProbabilityMean = depth%2 == 0 ? 1.0-maxBestProbabilityMean : maxBestProbabilityMean;
		
		//cerr << parent->getAlgorithmNodeData(SIMUST_DATA)->mean << endl;

		node = parent;
		depth--;
	}
}

int SIMuST::stoRecommed()
{
	double dMaxBestProbability = 0.0;
	int selectedIndex = 0;
	for(int i = 0; i < m_bf; i++){
		double bestProbability = ((stoData*)m_root->children[i]->getAlgorithmNodeData(SIMUST_DATA))->bestProbability;

		if(bestProbability > dMaxBestProbability){
			dMaxBestProbability = bestProbability;
			selectedIndex = i;
		}
	}


	return selectedIndex;
}

double SIMuST::calculateStudentTPDF( double tValue, double degreeOfFreedom )
{
	if(isnan(degreeOfFreedom)) degreeOfFreedom = numeric_limits<double>::max();

	students_t dist(degreeOfFreedom);
	return pdf(dist, fabs(tValue));
}

double SIMuST::calculateStudentTCDF( double tValue, double degreeOfFreedom )
{
	if(isnan(degreeOfFreedom)) degreeOfFreedom = numeric_limits<double>::max();

	students_t dist(degreeOfFreedom);
	return cdf(dist, fabs(tValue));
}

double SIMuST::calculateStudentTWinRate( double tValue, double degreeOfFreedom )
{
	if(isnan(degreeOfFreedom)) degreeOfFreedom = numeric_limits<double>::max();

	students_t dist(degreeOfFreedom);
	return cdf(dist, tValue);
}

double SIMuST::calculateStudentTBestProbability( int iIndex, vector<double> &vMeans, vector<double> &vDeviations, vector<int> &vMoveCounts )
{
	double dResult = 1.0;

	for(int i = 0; i < vMeans.size(); i++){
		if(i == iIndex) continue;
		if(vMeans[i] < 0) continue;

		double zValue = (vMeans[iIndex]-vMeans[i])/sqrt(vDeviations[iIndex]+vDeviations[i]);

		int n1 = vMoveCounts[iIndex];
		int n2 = vMoveCounts[i];
		double degreeOfFreedom = ((n1<=1) || (n2<=1)) ? 1 : pow(vDeviations[iIndex]+vDeviations[i],2.0)/(pow(vDeviations[iIndex],2.0)/(n1-1)+pow(vDeviations[i],2.0)/(n2-1));
		degreeOfFreedom = degreeOfFreedom <= 1 ? 1 : degreeOfFreedom;

		dResult *= calculateStudentTWinRate(zValue,degreeOfFreedom);
	}
	return dResult;
}

double SIMuST::showDetailStudentTBestProbability( int iIndex, vector<double> &vMeans, vector<double> &vDeviations, vector<int> &vMoveCounts )
{
	double dResult = 1.0;

	for(int i = 0; i < vMeans.size(); i++){
		if(i == iIndex) continue;
		if(vMeans[i] < 0) continue;

		double zValue = (vMeans[iIndex]-vMeans[i])/sqrt(vDeviations[iIndex]+vDeviations[i]);

		int n1 = vMoveCounts[iIndex];
		int n2 = vMoveCounts[i];
		double degreeOfFreedom = ((n1<=1) || (n2<=1)) ? 1 : pow(vDeviations[iIndex]+vDeviations[i],2.0)/(pow(vDeviations[iIndex],2.0)/(n1-1)+pow(vDeviations[i],2.0)/(n2-1));
		degreeOfFreedom = degreeOfFreedom <= 1 ? 1 : degreeOfFreedom;

		double temp = calculateStudentTWinRate(zValue,degreeOfFreedom);

		cerr << "\t\t\tzValue: " << zValue << "degreeOfFreedom: " << degreeOfFreedom << "\tWin rate of " << iIndex << " v.s. " << i << ": " << temp << " bestProb: " << dResult << endl; 
		dResult *= temp;
	}
	return dResult;
}

void SIMuST::dumpAllocationInfoToFile( string fileName )
{
	
	fstream f;
	f.open(fileName.c_str(), ios::out | ios::app);

	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		stoData* childStoData = (stoData*)m_root->children[i]->getAlgorithmNodeData(SIMUST_DATA);
		s.mean = childStoData->mean;
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	for(int i = 0; i < vNodeData.size(); i++){
		int originalIndex = vNodeData[i].originalIndex;
		f << vNodeData[i].mean << " " << m_root->children[originalIndex]->realScore << " ";
		stoData* childStoData = (stoData*)m_root->children[originalIndex]->getAlgorithmNodeData(SIMUST_DATA);
		f <<  childStoData->visitCount;
		f << " ";
	}

	f << endl;
}

int SIMuST::getNthBestMoveOfAlgorithms( int n )
{
	
	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		s.mean = ((stoData*)m_root->children[i]->getAlgorithmNodeData(SIMUST_DATA))->mean;
		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	if(n < 0 || n >= m_bf) return -1;
	else return vNodeData[n].originalIndex;
}

void SIMuST::showTree( treeNode* startNode, int depth )
{
	if(startNode == m_root || startNode->parent == m_root);
	else return;


	for(int i = 0; i < m_bf; i++){
		treeNode* nextNode = startNode->children[i];
		
		string emptyString = "";
		for(int k = 0; k < depth+1; k++) emptyString += string("      ");
		if(depth % 2 == 0) cerr << emptyString << "min node:\n";
		else cerr << emptyString << "max node:\n";

		stoData* nextNodeStoData = (stoData*)nextNode->getAlgorithmNodeData(SIMUST_DATA);
		if(nextNodeStoData){

			//cerr << "=========== Model Data===============\n";
			cerr << emptyString << "Model mean: " << nextNodeStoData->mean 
				<< " bestProbMean: " << nextNodeStoData->bestProbabilityMean 
				<< " visit: " << nextNodeStoData->visitCount 
				<< " bestProbability: " << nextNodeStoData->bestProbability 
				<< " sensitivity: " << nextNodeStoData->sensitivity
				<< " variance: " << nextNodeStoData->variance 
				<< " heuristic: " << nextNode->heuristic 
				<< " v-heuristic: " << nextNode->varianceHeuristic << "\n";
			//cerr << "======================================\n";
		}

		//nextNode->showInfo();
		 
		if(depth < m_depth-1) showTree(nextNode, depth+1);
		else;
	}
}
