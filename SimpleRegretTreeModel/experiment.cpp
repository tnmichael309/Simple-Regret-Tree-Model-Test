#include "experiment.h"


void experiment::run()
{
	int typeOfExp = 0;
	cout << "Enter type of experiments:\n <1> Stochastic Model <2> MCTS+Heuristic <3> SHOT+Heuristic" << endl
		<< " <4> Alpha-beta <5> MCTS+SIM+Heuristic <6> SHOT+SIM+Heuristic" << endl;
	cin >> typeOfExp;

	int nExps = 0;
	cout << "Enter # of experiments: " << endl;
	cin >> nExps;

	int nBudgets = 0;
	cout << "Enter # of budgets: " << endl;
	cin >> nBudgets;

	int bf = 0;
	cout << "Enter # of branches: " << endl;
	cin >> bf;

	int depth = 0;
	cout << "Enter depth of tree: " << endl;
	cin >> depth;

	double correctNum  = 0;
	double optionCorrectNum  = 0;
	double *topNaccuracies;
	topNaccuracies = new double[bf];
	for(int i = 0; i < bf; i++) topNaccuracies[i] = 0;


	//tree t(bf,depth,typeOfExp);
	tree* t;
	if(typeOfExp == 1) t = new SIMuST(bf,depth,typeOfExp);
	else if(typeOfExp == 2 || typeOfExp == 5) t = new mcts(bf,depth,typeOfExp);
	else if(typeOfExp == 3 || typeOfExp == 6) t = new shot(bf,depth,typeOfExp);
	else if(typeOfExp == 4) t = new ab(bf,depth,typeOfExp);

	
	for(int i = 0; i < nExps; i++){
		if(i%10 == 0) cout << "Running game: " << i << "                   \r\n" << endl;
		if(nExps>=10 && i%(nExps/10)==0) cerr << i/(nExps/10)*10 << "% ...";
		//cerr << "Processing experiment ID: " << i << endl;

		int realBestMove = t->getBestMove(false);

		int bestMove = 0;
		int optionalMove = 0;

		if(typeOfExp == 2 || typeOfExp == 3) t->runAlgorithm(nBudgets,true,depth,bestMove,optionalMove);
		else t->runAlgorithm(nBudgets,false,depth,bestMove,optionalMove);

		if(realBestMove == bestMove) correctNum++;
		else {
			/*t.getBestMove(true);
			cerr << "Suggested move: " << bestMove << endl;
			t.dumpTree();*/
		}

		if(realBestMove == optionalMove) optionCorrectNum++;
		//cerr << i << endl;

		for(int j = 0; j < bf; j++){
			if(realBestMove == t->getNthBestMoveOfAlgorithms(j)) topNaccuracies[j]++;
		}

		//t->dumpAllocationInfoToFile(to_string((long long)typeOfExp)+string("_")+to_string((long long)bf)+string("_")+to_string((long long)depth)+string("_")+to_string((long long)nBudgets));
		//cin.get();
		double max,min;
		max = 0.0;
		min = 1.0;
		t->resetTree(NULL,0,max,min);
	}
	cerr << "100% finished" << endl;

	

	double accumulatedAccuracy = 0.0;
	for(int j = 0; j < bf; j++){
		topNaccuracies[j] /= nExps;
		accumulatedAccuracy += topNaccuracies[j];
		cerr << accumulatedAccuracy << "\t";
	}
	cerr << endl;


	double errorProbability = (nExps-correctNum)/nExps;
	cout << "Experiment type: " << typeOfExp 
		<< " budget: " << nBudgets
		<< " repeats: " << nExps
		<< "\nError Prob.: " << errorProbability << "\n95% C.I. bound: "<< 1.96*sqrt(errorProbability*(1.0-errorProbability)/nExps) << endl;

	if(typeOfExp == 4){
		errorProbability = (nExps-optionCorrectNum)/nExps;
		cout << "Experiment type (post): " << typeOfExp 
			<< " budget: " << nBudgets
			<< " repeats: " << nExps
			<< "\nError Prob.: " << errorProbability << "\n95% C.I. bound: "<< 1.96*sqrt(errorProbability*(1.0-errorProbability)/nExps) << endl;
	}

	
}
