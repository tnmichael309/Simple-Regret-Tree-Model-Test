#include "derivedTreeNodeData.h"

mctsData::mctsData(){
	init();
}

mctsData::~mctsData(){};

void mctsData::init(){
	isExpanded = false;
	visitCount = 0;
	totalRewards = 0;
}

double mctsData::getWinRate(){
	if(visitCount == 0) return 0.0;
	return (double)totalRewards/(double)visitCount;
};

TREENODEDATA mctsData::getDataType()
{
	return MCTS_DATA;
}


shotData::shotData(){
	init();
};

shotData::~shotData(){};

void shotData::init(){
	isExpanded = false;
	visitCount = 0;
	totalRewards = 0;
};

double shotData::getWinRate(){
	if(visitCount == 0) return 0.0;
	return (double)totalRewards/(double)visitCount;
};

TREENODEDATA shotData::getDataType()
{
	return SHOT_DATA;
}


stoData::stoData(){
	init();
};

stoData::~stoData(){};

void stoData::init(){
	isExpanded = false;
	visitCount = 0;
	mean = 0;
	variance = 0;
	meanOfMean = 0;
	meanOfSquareSum = 0;
	bestProbability = 0;
	bestProbabilityMean = 0;
	entropy = 0;
};

TREENODEDATA stoData::getDataType()
{
	return SIMUST_DATA;
}
