#include "tree.h"


#define isPrune 1

class ab:public tree{
public:

	ab(int bf, int depth, int expType) : tree(bf,depth,expType){};
	~ab(){};

	void runAlgorithm(int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove);
	void dumpAllocationInfoToFile(string fileName);
	int getNthBestMoveOfAlgorithms(int n);

protected:
	void showTree(treeNode* startNode, int depth);

private:
	// Alpha beta
	double runSingleAlphaBeta(treeNode* startNode, int& consumedBudget, double alpha, double beta, int remaningDepth);

	vector<double> ab_pre_mean, ab_post_mean;
	vector<int> ab_pre_count, ab_post_count;
};