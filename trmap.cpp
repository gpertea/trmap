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
//#include <chrono>

#include "GBase.h"
#include "GVec.hh"
#include "gff.h"
#include "GArgs.h"
#include "GIntervalTree.h"

using std::cout;
using std::endl;

struct qInterval {
	std::string name;
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
			std::string name;
			int start;
			int end;
			std::string attrs;
			ss>>name>>start>>end;
			--start;
			std::getline(ss, attrs);
			queries.push_back({name,start,end,attrs});
		}
	}
	return queries;
}

char getOvlClass(GffObj* t, GffObj* r) {
 char c=0;

 return c;
}



struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: + strand, 2 : - strand
};

int main(int argc, char * const argv[]) {
//	const std::string usage=std::string("Usage: ")+argv[0]+"\n";
	const std::string usage = std::string("Positional arguments:\n")+
			"<ref_gff>    reference file name in GFF/BED format\n"+
			"<query_gff>  query file name in GFF/BED format or \"-\" to take from stdin\n";
			//"Options:\n"+
			//"-T		use interval trees";

	GArgs args(argc, argv, "ho:");
	args.printError(usage.c_str(), true);
	if (args.getOpt('h')) {
		cout << usage;
		exit(EXIT_SUCCESS);
	}

	std::unordered_map<std::string, GSTree> map_trees;

	const char* o_file = args.getOpt('o') ? args.getOpt('o') : "-";

	if (args.startNonOpt()!=2) {
		std::cerr << usage << "\nOnly " << args.startNonOpt() << " arguments provided (expected 2)\n";
		exit(1);
	}
	const char* ref_file = args.nextNonOpt();
	const char* q_file = args.nextNonOpt();

	FILE* fr=fopen(ref_file, "r");

	//always good to check if the file is actually there and can be read
	if (fr==NULL) GError("Error: could not open reference annotation file (%s)!\n", ref_file);
	const char* fext=getFileExt(ref_file);
	GffReader myR(fr, true, true);
	if (Gstricmp(fext, "bed")==0) myR.isBED();
	/*
	myR.readAll(false, true, true);
	for (int i=0; i<myR.gflst.Count(); i++) {
		GffObj* t=myR.gflst[i];
		map_trees[t->getGSeqName()].Insert(t);
	}
    */
	GffObj* t=NULL;
	GPVec<GffObj> toFree(true);
	while ((t=myR.readNext())!=NULL) {
		if (t->strand=='+')
		 map_trees[t->getGSeqName()].it[1].Insert(t);
		else if (t->strand=='-')
			map_trees[t->getGSeqName()].it[2].Insert(t);
		else map_trees[t->getGSeqName()].it[0].Insert(t);
		toFree.Add(t);
	}
	FILE* outFH=NULL;
	if (strcmp(o_file, "-")==0) outFH=stdout;
	                       else outFH=fopen(o_file, "w");
	FILE* fq=NULL;
	fext=NULL;
	if (strcmp(q_file,"-")==0) fq=stdin;
	else {
		fq=fopen(q_file, "r");
		if (fq==NULL)
			GError("Error: could not open query file (%s)!\n", q_file);
		fext=getFileExt(q_file);
	}
	GffReader myQ(fq, true, true);
	if (fext && Gstricmp(fext, "bed")==0) myQ.isBED();
	//myQ.readAll(false, true, true);
	//for (int i=0; i<myQ.gflst.Count(); i++) {
	//	GffObj* t=myQ.gflst[i];
	t=NULL;
	while ((t=myQ.readNext())!=NULL) {
		if (map_trees.count(t->getGSeqName())==0) continue;
		GVec<int> sidx;
		int v=0;
		sidx.Add(v); //always search the '.' strand
		if (t->strand=='+') { v=1; sidx.Add(v); }
		else if (t->strand=='-') { v=2; sidx.Add(v); }
		else { v=1; sidx.Add(v); v=2; sidx.Add(v); }
		for (int k=0;k<sidx.Count();++k) {
			TemplateStack<GSeg*> * enu = map_trees.at(t->getGSeqName()).it[sidx[k]].Enumerate(t->start, t->end);
			if(enu->Size()!=0) {
				fprintf(outFH, ">%s %s:%d-%d %c ", t->getID(), t->getGSeqName(), t->start, t->end, t->strand);
				t->printExonList(outFH);
				fprintf(outFH, "\n");
				for (int i=0; i<enu->Size(); ++i) {
					//static_cast<ObjInterval*>((*enu)[i])->obj->printGxf(oFile2);
					GffObj* r=(GffObj*)((*enu)[i]);
					r->printGTab(outFH);
				}
			}
			delete enu;
		}
		delete t;
	}

	fclose(outFH);
	return 0;
}
