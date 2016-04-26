#include "tree.h"

class mcts:public tree{
public:

	mcts(int bf, int depth, int expType) : tree(bf,depth,expType){};
	~mcts(){};

	void runAlgorithm(int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove);
	void dumpAllocationInfoToFile(string fileName);
	int getNthBestMoveOfAlgorithms(int n);

protected:
	void showTree(treeNode* startNode, int depth);

private:
	treeNode* mctsSelect(treeNode* node = NULL);
	double mctsSimulate(treeNode* node);
	void mctsUpdate(treeNode* node, double reward, treeNode* endNode = NULL);
	int mctsRecommed();

};