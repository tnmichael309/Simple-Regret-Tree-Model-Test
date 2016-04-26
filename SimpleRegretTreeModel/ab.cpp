#include "ab.h"


void ab::runAlgorithm( int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove )
{
	bestMove = 0;
	optionalMove = 0;

	int consumedBudget = 0;
	int searchDepth = floor(log((double)budget)/log((double)(m_bf)));
	if(searchDepth > maxDepth) searchDepth = maxDepth - 1;
	searchDepth--;
	if(searchDepth < 0) searchDepth = 0;

#if isPrune == 0
	searchDepth = 1;
#endif

	/*cout << searchDepth << endl;
	cin.get();
	cin.get();*/

	ab_pre_mean.clear();
	ab_pre_count.clear();
	ab_post_mean.clear();
	ab_post_count.clear();

	do{
		//cout << searchDepth << endl;
		
		searchDepth++;
		consumedBudget = 0;
		bestMove = optionalMove;

		ab_pre_mean.assign(ab_post_mean.begin(), ab_post_mean.end());
		ab_pre_count.assign(ab_post_count.begin(), ab_post_count.end());

		ab_post_mean.clear();
		ab_post_count.clear();

		optionalMove = runSingleAlphaBeta(m_root, consumedBudget, -1000, 1000, searchDepth);


		//cout << consumedBudget << endl;
		//cout << move1 << " " << move2 << endl;
	}while(consumedBudget < budget && searchDepth <= maxDepth);
	//cin.get();

}

double ab::runSingleAlphaBeta( treeNode* startNode, int& consumedBudget, double alpha, double beta, int remaningDepth )
{
	// find depth for 'startNode'
	treeNode* tempNode = startNode;
	int depth = 0;
	while(tempNode->parent){
		depth++;
		tempNode = tempNode->parent;
	};

	if(remaningDepth == 0 || depth == m_depth){
		consumedBudget++;

	
		double maxScore = -100;
		double bestMove = 0;
		if(depth == m_depth){
			maxScore = startNode->heuristic;
		}
		else{
			for(int i = 0; i < m_bf; i++){

				double score = startNode->children[i]->heuristic;

				//cout << score << "\t";

				if(depth%2==0);
				else score *= -1;

				if(score > maxScore){
					maxScore = score;
					bestMove = i;
				}
			}
		}


/*

		cout << endl << "depth: " << depth << "\t max score: ";
		cout << maxScore << "\n";*/
		if(startNode == m_root) return bestMove;
		else return maxScore;

	}else{
		double m = alpha;
		int bestMove = 0;
		for(int i = 0; i < m_bf; i++){

			int lastBudget = consumedBudget;
			double score = -1* runSingleAlphaBeta(startNode->children[i], consumedBudget, -1.0*beta, -1.0*m, remaningDepth-1);
			
			if(depth == 0) {
				ab_post_mean.push_back(score);
				ab_post_count.push_back(consumedBudget-lastBudget);
				lastBudget = consumedBudget;
			}

			//if(startNode == m_root) cout << score << "\t";

			if(score > m){
				m = score;
				bestMove = i;
			}
#if isPrune == 1
			if(m >= beta) return m;
#endif
		}

		if(startNode == m_root){
			return bestMove;
		}else{
			return m;
		}
	}
}

void ab::dumpAllocationInfoToFile( string fileName )
{
	static bool firstCompute = true;
	fstream f;
	f.open(fileName.c_str(), ios::out | ios::app);

	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		if(firstCompute){
			s.mean = ab_pre_mean[i];
		}else if(!firstCompute){
			s.mean = ab_post_mean[i];
		}

		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	for(int i = 0; i < vNodeData.size(); i++){
		int originalIndex = vNodeData[i].originalIndex;
		f << vNodeData[i].mean << " " << m_root->children[originalIndex]->realScore << " ";

		if(firstCompute){
			f << ab_pre_count[originalIndex];
		}else if(!firstCompute){
			f << ab_post_count[originalIndex];
		}

		f << " ";
	}

	if(firstCompute){
		firstCompute = false;
		dumpAllocationInfoToFile(fileName+string("_psot"));
	}else if(!firstCompute){
		firstCompute = true;
	}

	f << endl;
}

int ab::getNthBestMoveOfAlgorithms( int n )
{
	static bool firstCompute = true;
	vector<simpleNodeData> vNodeData;
	for(int i = 0; i < m_bf; i++){

		simpleNodeData s;
		if(firstCompute){
			s.mean = ab_pre_mean[i];
		}else if(!firstCompute){
			s.mean = ab_post_mean[i];
		}

		s.originalIndex = i;

		vNodeData.push_back(s);
	}

	sort(vNodeData.begin(),vNodeData.end(),simpleNodeData::sort);

	if(firstCompute){
		firstCompute = false;
		getNthBestMoveOfAlgorithms(n);
	}else if(!firstCompute){
		firstCompute = true;
	}

	if(n < 0 || n >= m_bf) return -1;
	else return vNodeData[n].originalIndex;
}

void ab::showTree( treeNode* startNode, int depth )
{

}
