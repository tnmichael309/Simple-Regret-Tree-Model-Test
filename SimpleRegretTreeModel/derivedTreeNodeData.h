#include "baseTreeNodeData.h"

#pragma pack(push, 1)
class mctsData : public baseTreeNodeData{
public:

	mctsData();
	~mctsData();
	double getWinRate();

	// mcts
	bool isExpanded;
	int visitCount;
	double totalRewards;

	void init();
	TREENODEDATA getDataType();
};

#pragma pack(push, 1)
class shotData : public baseTreeNodeData{
public:
	shotData();
	~shotData();

	double getWinRate();

	// shot
	bool isExpanded;
	int visitCount;
	double totalRewards;

	void init();
	TREENODEDATA getDataType();
};

#pragma pack(push, 1)
class stoData : public baseTreeNodeData{
public:
	stoData();
	~stoData();

	// sto
	bool isExpanded;
	int visitCount;
	double mean;
	double variance;
	double meanOfMean;
	double meanOfSquareSum;
	double bestProbability;
	double bestProbabilityMean;
	double sensitivity;
	double entropy;
	bool isBest;

	void init();
	TREENODEDATA getDataType();
};