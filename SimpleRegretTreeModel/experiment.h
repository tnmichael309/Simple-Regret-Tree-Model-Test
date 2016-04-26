#include "SIMuST.h"
#include "ab.h"
#include "mcts.h"
#include "shot.h"

#include <iostream>
using namespace std;

class experiment{
public:

	// disallow copying
	experiment(experiment const&); //= delete; /*c++11*/
	void operator=(experiment const&);// = delete; /*c++11*/ 

	static experiment& getInstance(){
		static experiment m_experiment;

		return m_experiment;
	};

	void run();

private:
	experiment(){};


};