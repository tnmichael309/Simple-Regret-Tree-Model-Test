#include "tree.h"

#define isResetVariance 0
#define isUsePropagate 0
#define isUseMaxOfflineVariance 1
#define isNormalizeVariance 0

class SIMuST : public tree{
public:

	SIMuST(int bf, int depth, int expType) : tree(bf,depth,expType){};
	~SIMuST(){};

	void runAlgorithm(int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove);
	void dumpAllocationInfoToFile(string fileName);
	int getNthBestMoveOfAlgorithms(int n);

protected:
	void showTree(treeNode* startNode, int depth);

private:


	//Stochastic Model
	treeNode* stoSelect();
	void stoUpdate(treeNode* node);
	int stoRecommed();


	////////
	//student-t's helper functions
	////////
	double calculateStudentTPDF(double tValue, double degreeOfFreedom);
	// separability
	double calculateStudentTCDF(double tValue, double degreeOfFreedom);
	// win rate
	double calculateStudentTWinRate(double tValue, double degreeOfFreedom);
	double calculateStudentTBestProbability(int iIndex, vector<double> &vMeans, vector<double> &vDeviations, vector<int> &vMoveCounts);
	double showDetailStudentTBestProbability(int iIndex, vector<double> &vMeans, vector<double> &vDeviations, vector<int> &vMoveCounts);
};




