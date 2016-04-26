#ifndef TREE_H
#define TREE_H


#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <boost/math/distributions/students_t.hpp>

#include "derivedTreeNodeData.h"

using namespace std;
using namespace boost::math;

#define isUseNormal 1
#define isZeroOrOneHeuristic 0

#define minVariance 0.0000000001

class simpleNodeData{
public:
	double mean;
	int originalIndex;
	static bool sort(simpleNodeData& node1, simpleNodeData& node2);
};


#pragma pack(push, 1)
class treeNode{
public:

	treeNode();
	~treeNode();

	void init();
	baseTreeNodeData* getAlgorithmNodeData(TREENODEDATA treeNodeData);
	void updateAlgorithmNodeData(baseTreeNodeData* nodeData,TREENODEDATA treeNodeData);
	void genMoveValue(bool isPositive);
	void genSpecialMoveValue(bool isPositive, bool isIncreasingRewards, int depth);
	void showInfo();

	double increaseRate;
	int moveValue;
	int cumulativeMoveValue;
	double heuristic;
	double varianceHeuristic;
	double realScore;
	treeNode* parent;
	vector<treeNode*> children;

	vector<baseTreeNodeData*> m_vAlgorithmData;
	//mctsData* getAlgorithmNodeData(MCTS_DATA);
	//shotData* getAlgorithmNodeData(SHOT_DATA);
	//stoData* getAlgorithmNodeData(SIMUST_DATA);
};

class tree{
public:
	tree(int bf, int depth, int expType);
	~tree();

	//// implement by derived class (algorithm)
	virtual void runAlgorithm(int budget, bool isUseHeuristic, int maxDepth, int& bestMove, int& optionalMove){};
	virtual void dumpAllocationInfoToFile(string fileName){};
	virtual int getNthBestMoveOfAlgorithms(int n){return 0;};
	//////////////////////////////////////////

	int getBestMove(bool isShowInfo);
	void resetTree(treeNode* startNode, int depth, double& max, double& min);
	void dumpTree();

protected:
	treeNode* m_root;
	int m_depth;
	int m_bf;
	int m_expType;

	//// implement by derived class (algorithm)
	virtual void showTree(treeNode* startNode, int depth){};
	//////////////////////////////////////////

private:
	void buildTree(treeNode* startNode, int depth, double& max, double& min);
	void deleteTree(treeNode* startNode, int depth);
	int traverseForBestMove(treeNode* startNode, int depth);


};

#endif