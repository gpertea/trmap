//============================================================================
// Name        : mapped_nc.cpp
// Author      : 
// Version     :
// Copyright   : right
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <unordered_map>

#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>

#include "GBase.h"
#include "gff.h"
#include "GArgs.h"
#include "GIntervalTree.h"

using std::cout;
using std::endl;

/*
struct ObjInterval: public GInterval {
	GffObj* obj;
	//virtual int GetLowPoint() const{return obj->start;}
	//virtual int GetHighPoint() const{return obj->end;}
	uint getStart() const {return obj->start;}
	uint getEnd() const {return obj->end;}
	ObjInterval(GffObj* obj): obj(obj){}
};
*/
struct qInterval {
	std::string seg;
	int start;
	int end;
	std::string attrs; //everything else
};

//bool cmp_tabix(const GffObj& a, const GffObj& b){
//	return (a.start==b.start)?(a.end==b.end?a.obj->getID()<b.obj->getID():(a.end<b.end)):(a.start<b.start);
//}

std::vector<qInterval> readQueries(std::istream& input) { //RVO should make this OK
	std::vector<qInterval> queries;
	std::string line;
	while (std::getline(input, line)){
		if (line.length()!=0 && line.front()!='#'){ //neither blank nor comment, could use continue
			std::stringstream ss;
			ss.str(line);
			std::string seg;
			int start;
			int end;
			std::string attrs;
			ss>>seg>>start>>end;
			std::getline(ss, attrs);
			queries.push_back({seg,start,end,attrs});
		}
	}
	return queries;
}

int main(int argc, char * const argv[]) {
//	const std::string usage=std::string("Usage: ")+argv[0]+"\n";
	const std::string usage = std::string("Positional arguments:\n")+
			"<input_gff>	reference file in GFF format\n"+
			"<query_BED>	query file in BED format, or \"-\" to take from stdin\n";
			//"Options:\n"+
			//"-T		use interval trees";

	GArgs args(argc, argv, "ho:");
	args.printError(usage.c_str(), true);
	if (args.getOpt('h')) {
		cout << usage;
		exit(EXIT_SUCCESS);
	}

	//bool doNCLorITT = !args.getOpt('T');
	//bool doNCLorITT = true;
	std::unordered_map<std::string, GIntervalTree> map_trees;
	//std::unordered_map<std::string, NCList<ObjInterval> > map_nc;

	const char* o_file = args.getOpt('o') ? args.getOpt('o') : "mapped_ov.tab";

//	if(args.isError()) {
//		std::cerr << "Error in arg #" << args.isError() << "\n";
//		exit(EXIT_FAILURE);
//	}

	if (args.startNonOpt()!=2) {
		std::cerr << "Only " << args.startNonOpt() << " arguments provided (expected 2)" << "\n";
		//print usage here?
		exit(EXIT_FAILURE);
	}
	const char* gff_file = args.nextNonOpt();
	const char* q_file = args.nextNonOpt();

	FILE* f=fopen(gff_file, "r");

	//always good to check if the file is actually there and can be read
	if (f==NULL) GError("Error: could not open reference annotation file (%s)!\n", gff_file);

	GffReader myR(f, true, true);
	myR.readAll(false, true, true);
	//fclose(f); // don't do it, myR will do it when going out of scope

//	cout << myR.gflst[0]->getGSeqName() << "\n";

	//Interval tree
	for (int i=0; i<myR.gflst.Count(); i++) {
		GffObj* transcript=myR.gflst[i];
		map_trees[transcript->getGSeqName()].Insert(transcript);
	}

//		cout << map_trees.size() << "\n";
//		for (std::unordered_map<std::string, IntervalTree>::const_iterator it = map_trees.begin(); it!=map_trees.end(); ++it){
//			cout << it->first << "\n";
//		}

	std::vector<qInterval> queries;
	if(*q_file=='-') {
		queries = readQueries(std::cin);
	} else {
		std::ifstream query_file(q_file);
		if (query_file.is_open()){
			queries = readQueries(query_file);
		} else {
			std::cerr << "query file didn't open";
			exit(EXIT_FAILURE);
		}
//		query_file.close(); //will be done automagically
	}

	FILE* oFile2=fopen(o_file, "w");

	//std::chrono::time_point<std::chrono::high_resolution_clock> pre_ov = std::chrono::high_resolution_clock::now();
	for (std::vector<qInterval>::const_iterator it=queries.begin(); it!=queries.end(); ++it){
			//		std::vector<ObjInterval> overlaps;
			TemplateStack<GSeg*> * enu = map_trees.at(it->seg).Enumerate(it->start, it->end);
			if(enu->Size()!=0){
				fprintf(oFile2, "##Qry|%s: %i-%i %s\n", it->seg.c_str(), it->start, it->end, it->attrs.c_str());
				//			oFile << "##Qry|" << it->seg << ":" << it->start << "-" << it->end << it->attrs << "\n";
				//			oFile.flush();
				for (int i=0; i<enu->Size(); ++i) {
					//static_cast<ObjInterval*>((*enu)[i])->obj->printGxf(oFile2);
					((GffObj*)((*enu)[i]))->printGxf(oFile2);
				}
			}
			delete enu;
	}
	//std::chrono::duration<double> ov_time = std::chrono::high_resolution_clock::now()-pre_ov;
	fclose(oFile2);

	//std::cerr << (doNCLorITT?"NCList time: ":"tree time: ") << ov_time.count() << "\n";
	return 0;
}
