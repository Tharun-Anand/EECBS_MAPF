/* Copyright (C) Jiaoyang Li
* Unauthorized copying of this file, via any medium is strictly prohibited
* Confidential
* Written by Jiaoyang Li <jiaoyanl@usc.edu>, May 2020
*/

/*driver.cpp
* Solve a MAPF instance on 2D grids.
*/
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include "ECBS.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <string>

using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

using namespace std;

struct Agent {
  int id;
  vector<pair<int, int>> path;
};


/* Main function */
int main(int argc, char** argv)
{   

	//auto start_time1 = Clock::now();
	namespace po = boost::program_options;
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")

		// params for the input instance and experiment settings
		("map,m", po::value<string>()->required(), "input file for map")
		("agents,a", po::value<string>()->required(), "input file for agents")
		("output,o", po::value<string>(), "output file for statistics")
		("outputPaths", po::value<string>(), "output file for paths")
		("agentNum,k", po::value<int>()->default_value(0), "number of agents")
		("cutoffTime,t", po::value<double>()->default_value(7200), "cutoff time (seconds)")
		("screen,s", po::value<int>()->default_value(1), "screen option (0: none; 1: results; 2:all)")
		("stats", po::value<bool>()->default_value(false), "write to files some detailed statistics")

		// params for CBS node selection strategies
		("highLevelSolver", po::value<string>()->default_value("EES"), "the high-level solver (A*, A*eps, EES, NEW)")
		("lowLevelSolver", po::value<bool>()->default_value(true), "using suboptimal solver in the low level")
		("inadmissibleH", po::value<string>()->default_value("Global"), "inadmissible heuristics (Zero, Global, Path, Local, Conflict)")
		("suboptimality", po::value<double>()->default_value(1.2), "suboptimality bound")

		// params for CBS improvement
		("heuristics", po::value<string>()->default_value("WDG"), "admissible heuristics for the high-level search (Zero, CG,DG, WDG)")
		("prioritizingConflicts", po::value<bool>()->default_value(true), "conflict prioirtization. If true, conflictSelection is used as a tie-breaking rule.")
		("bypass", po::value<bool>()->default_value(true), "Bypass1")
		("disjointSplitting", po::value<bool>()->default_value(false), "disjoint splitting")
		("rectangleReasoning", po::value<bool>()->default_value(true), "rectangle reasoning")
		("corridorReasoning", po::value<bool>()->default_value(true), "corridor reasoning")
		("targetReasoning", po::value<bool>()->default_value(true), "target reasoning")
		("sipp", po::value<bool>()->default_value(0), "using SIPPS as the low-level solver")
		("restart", po::value<int>()->default_value(0), "rapid random restart times")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	po::notify(vm);
	if (vm["suboptimality"].as<double>() < 1)
	{
		cerr << "Suboptimal bound should be at least 1!" << endl;
		return -1;
	}

	high_level_solver_type s;
	if (vm["highLevelSolver"].as<string>() == "A*")
		s = high_level_solver_type::ASTAR;
	else if (vm["highLevelSolver"].as<string>() == "A*eps")
		s = high_level_solver_type::ASTAREPS;
	else if (vm["highLevelSolver"].as<string>() == "EES")
		s = high_level_solver_type::EES;
	else if (vm["highLevelSolver"].as<string>() == "NEW")
		s = high_level_solver_type::NEW;
	else
	{
		cout << "WRONG high level solver!" << endl;
		return -1;
	}

	if (s == high_level_solver_type::ASTAR && vm["suboptimality"].as<double>() > 1)
	{
		cerr << "A* cannot perform suboptimal search!" << endl;
		return -1;
	}

    heuristics_type h;
	if (vm["heuristics"].as<string>() == "Zero")
		h = heuristics_type::ZERO;
	else if (vm["heuristics"].as<string>() == "CG")
		h = heuristics_type::CG;
	else if (vm["heuristics"].as<string>() == "DG")
		h = heuristics_type::DG;
	else if (vm["heuristics"].as<string>() == "WDG")
		h = heuristics_type::WDG;
	else
	{
		cout << "WRONG heuristics strategy!" << endl;
		return -1;
	}

    if ((h == heuristics_type::CG || h == heuristics_type::DG) && vm["lowLevelSolver"].as<bool>())
    {
        cerr << "CG or DG heuristics do not work with low level of suboptimal search!" << endl;
        return -1;
    }

	heuristics_type h_hat; // inadmissible heuristics
	if (s == high_level_solver_type::ASTAR ||
	    s == high_level_solver_type::ASTAREPS ||
	    vm["inadmissibleH"].as<string>() == "Zero")
		h_hat = heuristics_type::ZERO;
	else if (vm["inadmissibleH"].as<string>() == "Global")
		h_hat = heuristics_type::GLOBAL;
	else if (vm["inadmissibleH"].as<string>() == "Path")
		h_hat = heuristics_type::PATH;
	else if (vm["inadmissibleH"].as<string>() == "Local")
		h_hat = heuristics_type::LOCAL;
	else if (vm["inadmissibleH"].as<string>() == "Conflict")
		h_hat = heuristics_type::CONFLICT;
	else
	{
		cout << "WRONG inadmissible heuristics strategy!" << endl;
		return -1;
	}

	conflict_selection conflict = conflict_selection::EARLIEST;
	node_selection n = node_selection::NODE_CONFLICTPAIRS;


	srand((int)time(0));

	///////////////////////////////////////////////////////////////////////////
	// load the instance
	 for(int t=0;t<1;t++)
	  {
	Instance instance(vm["map"].as<string>(), vm["agents"].as<string>(),
		vm["agentNum"].as<int>());
 
	srand(0);
	int runs = 1 + vm["restart"].as<int>();
	cout<<"runs:"<<runs<<"##";
	//////////////////////////////////////////////////////////////////////
    // initialize the solver
   

	if (vm["lowLevelSolver"].as<bool>())
    {   

	
        ECBS ecbs(instance, vm["sipp"].as<bool>(), vm["screen"].as<int>());
        ecbs.setPrioritizeConflicts(vm["prioritizingConflicts"].as<bool>());
        ecbs.setDisjointSplitting(vm["disjointSplitting"].as<bool>());
        ecbs.setBypass(vm["bypass"].as<bool>());
        ecbs.setRectangleReasoning(vm["rectangleReasoning"].as<bool>());
        ecbs.setCorridorReasoning(vm["corridorReasoning"].as<bool>());
        ecbs.setHeuristicType(h, h_hat);
        ecbs.setTargetReasoning(vm["targetReasoning"].as<bool>());
        ecbs.setMutexReasoning(false);
        ecbs.setConflictSelectionRule(conflict);
        ecbs.setNodeSelectionRule(n);
        ecbs.setSavingStats(vm["stats"].as<bool>());
        ecbs.setHighLevelSolver(s, vm["suboptimality"].as<double>());
        //////////////////////////////////////////////////////////////////////
        // run
        double runtime = 0;
        int lowerbound = 0;
	  // insert goals list
      
	
            
			for (int i = 0; i < 1; i++)
			{
				ecbs.clear();
				//ecbs.solve(vm["cutoffTime"].as<double>() / runs, lowerbound);
				ecbs.solve(vm["cutoffTime"].as<double>() / 1, lowerbound);

				runtime += ecbs.runtime;
				if (ecbs.solution_found)
					break;
				lowerbound = ecbs.getLowerBound();
				ecbs.randomRoot = true;
				cout << "Failed to find solutions in Run " << i << endl;
			}
			ecbs.runtime = runtime;
			if (vm.count("output"))
				ecbs.saveResults(vm["output"].as<string>(), vm["agents"].as<string>());
			if (ecbs.solution_found && vm.count("outputPaths"))
				ecbs.savePaths(vm["outputPaths"].as<string>());
			/*size_t pos = vm["output"].as<string>().rfind('.');      // position of the file extension
			string output_name = vm["output"].as<string>().substr(0, pos);     // get the name without extension
			cbs.saveCT(output_name); // for debug*/
			if (vm["stats"].as<bool>())
				ecbs.saveStats(vm["output"].as<string>(), vm["agents"].as<string>());
			ecbs.clearSearchEngines();
           
        auto start_time1 = Clock::now();
		ifstream file("/home/sslv1/EECBS/paths.txt");
		if (!file.is_open()) {
			cerr << "Error opening file paths.txt" << endl;
			return 1;
		}

		// Stores information for each agent
		vector<Agent> agents;

		string line;
		int x=0;
		while (getline(file, line)) {
		if (line.empty()) {
			continue;
		}

		// Extract agent ID
		istringstream iss(line);
		cout << line << endl;  // Print the full line for debugging (optional)
		string agentStr;
		iss >> agentStr;
		int agentId;
		sscanf(agentStr.c_str(), "Agent %d:", &agentId);

		Agent agent;
		agent.id = x;
		x++;
		//cout << "$$" << agentId << "$$" << endl;  // Print agent ID for debugging (optional)
        std::regex pattern(R"(\d+,\d+)");
		// Read path points
	    std::smatch match;
		std::string remaining = line;

		// Iterate through all matches of the regex pattern
		while (std::regex_search(remaining, match, pattern)) {
			std::string coordinatesStr = match[0];

			// Extract x and y values
			std::stringstream coordinateStream(coordinatesStr);
			int x, y;
			std::getline(coordinateStream, coordinatesStr, ',');
			x = std::stoi(coordinatesStr);
			std::getline(coordinateStream, coordinatesStr, ',');
			y = std::stoi(coordinatesStr);

			// Add the coordinates to the vector
            agent.path.push_back(make_pair(x, y));
			// Update remaining string to search for next coordinates
			remaining = match.suffix().str();
		}

		agents.push_back(agent);
		}

		file.close();

		// Find agent with shortest path
		int shortestPath = INT_MAX;
		int shortestAgentId = -1;
		for (const Agent& agent : agents) {
			int pathLength = agent.path.size();
			if (pathLength < shortestPath) {
			shortestPath = pathLength;
			shortestAgentId = agent.id;
			}
		}

		// Print result for agent with shortest path
		if (shortestAgentId != -1) {
			cout << "Agent with shortest path: Agent " << shortestAgentId << endl;
			cout << "Path length: " << shortestPath << endl;
		} else {
			cout << "No agents found in the file." << endl;
		}
        ifstream inputFile("/home/sslv1/EECBS/random-32-32-20_120.scen");
			if (!inputFile.is_open()) {
			cerr << "Error opening file random-32-32-20.map" << endl;
			return 1;
			}
            std::string tempFileName="random-32-32-20_modified.scen";
			ofstream outputFile(tempFileName);  // Create a new file for modified content
			if (!outputFile.is_open()) {
			cerr << "Error creating output file random-32-32-20_modified.map" << endl;
			inputFile.close();
			return 1;
			}

			string line1;
			int y=0;
			int c=0;
			while (getline(inputFile, line1)) {
            c++;
			// Split the line by whitespace
			istringstream iss(line1);
			//cout<<"!!"<<line1<<"@@";
			vector<string> tokens;
			string token;
			while (getline(iss, token, ' ')) {
			//cout<<"##"<<tokens.size()<<"$$";
			//cout<<"iteration no:"<<y<<"finish";
			y++;
			if(y!=2)
			{
			 if(y>=3){
			 //cout<<"$$"<<std::to_string(agents[y-3].path[shortestPath-1].first)[0]<<"##";
			 token[32]=std::to_string(agents[y-3].path[shortestPath-1].first)[0];
			 token[33]=std::to_string(agents[y-3].path[shortestPath-1].first)[1];
			 token[29]=std::to_string(agents[y-3].path[shortestPath-1].second)[0];
			 token[30]=std::to_string(agents[y-3].path[shortestPath-1].second)[1];
            if(y-3==shortestAgentId)
			{ 

			   cout<<"35:"<<token[35];
			   cout<<"36:"<<token[36];
			   cout<<"37:"<<token[37];
			   cout<<"38:"<<token[38];
			   if(t==0){
				token[35]='1';
			    token[36]='1';
				token[38]='1';
			    token[39]='1';}
			    
				else if(t==1){
				token[35]='1';
			    token[36]='0';
				token[38]='1';
			    token[39]='0';

				}
			}

			 }
             
			 
			 //cout<<"$$"<<y<<"##";
			//cout<<"##"<<agent.path[shortestPath-1].first<<"$$";
			}
			tokens.push_back(token);
			}

			// Modify 5th and 6th columns (assuming indexing starts from 0)
			

			// Rebuild the line with modified tokens
			string modifiedLine;
			//cout<<"modifiedLine";
			for (int i = 0; i < tokens.size(); i++) {
		    cout<<"##"<<modifiedLine<<"$$";
			 modifiedLine += tokens[i] + (i != tokens.size() - 1 ? "" : " ");
			//cout<<"##"<<modifiedLine<<"i:"<<i<<"$$";

			}
            
			// Write the modified line to the output file
			
			 //cout<<"##"<<c<<"$$";
			 outputFile << modifiedLine << endl;
			}
            cout<<"OUTSIDE";
			inputFile.close();
			outputFile.close();
			if (rename(tempFileName.c_str(), "random-32-32-20_120.scen") != 0) {
				cerr << "Error renaming temporary file" << endl;
				return 1;
			}

			std::ifstream inputFile1("random-32-32-20_120.scen");

			if (!inputFile1.is_open()) {
				std::cerr << "Error opening file" << std::endl;
				return 1;
			}

			int lineCount = 0;
			std::string line2;

			// Read the file line by line
			while (std::getline(inputFile1, line2)) {
				lineCount++;
			}

			inputFile1.close();

			// Print the number of lines
			std::cout << "The file has " << lineCount << " lines." << std::endl;

			cout << "Modified content written to random-32-32-20_modified.map" << endl;
		// Prepare and write robot locations to CSV file
		ofstream csvFile("robot_locations.csv");
		if (csvFile.is_open()) {
			csvFile << "Agent,X,Y" << endl;
			for (const Agent& agent : agents) {
				cout<<"##"<<agent.path[shortestPath-1].first<<"$$";
				csvFile << "Agent " << agent.id << "," << agent.path[shortestPath-1].first << "," << agent.path[shortestPath-1].second << endl;
			}
			csvFile.close();
			cout << "Robot locations written to robot_locations.csv" << endl;
		} else {
			cerr << "Error creating robot_locations.csv" << endl;
		}
		auto end_time1 = Clock::now();
    auto duration1=std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time1).count();
    cout << "Time taken for reading: " << duration1 << " milliseconds\n";

	  }
	
    }
	
	/*
    else
    {
        CBS cbs(instance, vm["sipp"].as<bool>(), vm["screen"].as<int>());
        cbs.setPrioritizeConflicts(vm["prioritizingConflicts"].as<bool>());
        cbs.setDisjointSplitting(vm["disjointSplitting"].as<bool>());
        cbs.setBypass(vm["bypass"].as<bool>());
        cbs.setRectangleReasoning(vm["rectangleReasoning"].as<bool>());
        cbs.setCorridorReasoning(vm["corridorReasoning"].as<bool>());
        cbs.setHeuristicType(h, h_hat);
        cbs.setTargetReasoning(vm["targetReasoning"].as<bool>());
        cbs.setMutexReasoning(false);
        cbs.setConflictSelectionRule(conflict);
        cbs.setNodeSelectionRule(n);
        cbs.setSavingStats(vm["stats"].as<bool>());
        cbs.setHighLevelSolver(s, vm["suboptimality"].as<double>());
        //////////////////////////////////////////////////////////////////////
        // run
        double runtime = 0;
        int lowerbound = 0;
        for (int i = 0; i < 1; i++)
        {
            cbs.clear();
            //cbs.solve(vm["cutoffTime"].as<double>() / runs, lowerbound);
			cbs.solve(vm["cutoffTime"].as<double>() / 1, lowerbound);

            runtime += cbs.runtime;
            if (cbs.solution_found)
                break;
            lowerbound = cbs.getLowerBound();
            cbs.randomRoot = true;
            cout << "Failed to find solutions in Run " << i << endl;
        }
        cbs.runtime = runtime;
        if (vm.count("output"))
            cbs.saveResults(vm["output"].as<string>(), vm["agents"].as<string>());
        if (cbs.solution_found && vm.count("outputPaths"))
            cbs.savePaths(vm["outputPaths"].as<string>());
        if (vm["stats"].as<bool>())
            cbs.saveStats(vm["output"].as<string>(), vm["agents"].as<string>());
        cbs.clearSearchEngines();
    }
    */

	return 0;

}