/**
 * @brief Heuristic manager
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2022 Domenico Salvagnin
 */

#ifndef HEURMGR_H
#define HEURMGR_H

#include "mip.h"
#include "worker.h"
#include <map>
#include <string>

struct RunType
{
public:
	bool propagate;
	bool repair;
	bool backtrackOnInfeas;
};


static std::map<std::string, RunType> runTypes {
	{"dfs",      {true,  false, true}},
	{"dfsrep",   {true,  true,  true}},
	{"dive",     {false, true,  false}},
	{"diveprop", {true,  true,  false}},
	{"walkmip",  {false,  false,  false}}
};


/* Dependencies */
const int DEP_NONE      =  0;
const int DEP_COREPOINT =  1;
const int DEP_VERTEX    =  2;
const int DEP_ZEROOBJ   =  4;

struct Strategy
{
public:
	std::string varSelect;
	std::string valueSelect;
	int depends;
};


static std::map<std::string, Strategy> strategies {
	{"random",    {"typec",   "random",            DEP_NONE}},
	{"random2",   {"random",  "random",            DEP_NONE}},
	{"goodobj",   {"type",    "good_obj",          DEP_NONE}},
	{"badobj",    {"type",    "bad_obj",           DEP_NONE}},
	{"goodobjcl", {"ccover_goodobj",  "up",        DEP_NONE}},
	{"badobjcl",  {"ccover_badobj",   "up",        DEP_NONE}},
	{"locks",     {"LR",      "loosedyn",          DEP_NONE}},
	{"locks2",    {"locks",   "loosedyn",          DEP_NONE}},
	{"cliques",   {"cliques", "up",                DEP_COREPOINT | DEP_ZEROOBJ}},
	{"zerocore",  {"typec",   "zeroobj_corepoint", DEP_COREPOINT | DEP_ZEROOBJ}},
	{"core",      {"typec",   "corepoint",         DEP_COREPOINT}},
	{"cliques2",  {"cliques2","up",                DEP_VERTEX | DEP_ZEROOBJ}},
	{"zerolp",    {"typec",   "zeroobj_lp",        DEP_VERTEX | DEP_ZEROOBJ}},
	{"lp",        {"typec",   "lp",                DEP_VERTEX}}
};


/* Run a given heuristic */
void runHeur(
	const std::string& run,
	const std::string& strat,
	MIPData& data,
	WorkerDataManager& wManager,
	Params params,
	int trial = 0);

#endif /* HEURMGR_H */
