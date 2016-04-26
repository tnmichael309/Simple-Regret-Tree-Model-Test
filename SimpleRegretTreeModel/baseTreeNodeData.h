#include <iostream>

enum TREENODEDATA{
	SIMUST_DATA,
	MCTS_DATA,
	SHOT_DATA,

	MAX_TREE_NODE_DATA_NUM
};

class baseTreeNodeData{
public:
	virtual void init() = 0;
	virtual TREENODEDATA getDataType() = 0;


private:

};