#include "tree.h"

class shot:public tree{
public:

	shot(int bf, int depth, int expType) : tree(bf,depth,expType){};
	~shot(){};

	void runAlgorithm(int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove);
	void dumpAllocationInfoToFile(string fileName);
	int getNthBestMoveOfAlgorithms(int n);

protected:
	void showTree(treeNode* startNode, int depth);

private:
	double allocateABudget(treeNode* startNode, bool isUseHeuristic);
	int doSHOT(treeNode* startNode, int budget, int& budgetUsed, int& playouts, double& rewards, bool isUseHeuristic);
	double shotSimulate(treeNode* node);
};