/* 
 * ControlMenu.cc, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "ControlMenu.hh"
#include <utility>
#include <cassert>

bool StreamInteractor::menutils_CheckStackSize(unsigned int n) {
	assert(mydeque && mystack);
	menutils_StackSize(this);
	unsigned int nstack = popInt();
	if(nstack < n) {
		mydeque->push_front(NameSelector::barf_control+" Insufficient Arguments ["+itos(n-nstack)+"]");
		return false;
	}
	return true;
}

//---- InputRequester

InputRequester InputRequester::exitMenu("Exit Menu", &menutils_Exit);

InputRequester::InputRequester(string d, void (*f)(StreamInteractor*), StreamInteractor* fObj):
NamedInteractor(d), myFuncObject(fObj), myFunc(f) { }

void InputRequester::doIt() {

	assert(mydeque && mystack);
	
	for(unsigned int i=0; i<inputFilters.size(); i++) {
		if(inputFilters[i]) {
			// call input filter if set
			inputFilters[i]->mydeque = mydeque;
			inputFilters[i]->mystack = mystack;
			inputFilters[i]->doIt();
		} else if(mydeque->size()) {
			// pull next item from deque, if available
			mystack->push(popStringD());
		} else {
			// prompt stdin for entry
			char indata[1024];
			if(argDescrips[i].size())
				printf("\n// %s\n",argDescrips[i].c_str());
			printf("%s",argNames[i].c_str());
			if(defaultArgs[i].size())
				printf(" [%s]",defaultArgs[i].c_str());
			printf(" > ");
			assert(fgets(indata,1024,stdin));
			string argin = strip(indata);
			if(!argin.size())
				mystack->push(defaultArgs[i]);
			else
				mystack->push(argin);
		}
	}	
	
	// function call using collected results
	if(myFunc) {
		if(myFuncObject) {
			myFuncObject->mydeque = mydeque;
			myFuncObject->mystack = mystack;
			(*myFunc)(myFuncObject);
		} else (*myFunc)(this);
	}
}

void InputRequester::addArg(const std::string& s, const std::string& dflt, const std::string& descrip, NamedInteractor* filter) {
	argNames.push_back(s);
	argDescrips.push_back(descrip);
	defaultArgs.push_back(dflt);
	inputFilters.push_back(filter);
}

void InputRequester::setArgOpts(unsigned int i, string s, string dflt, NamedInteractor* filter) { 
	assert(i<argNames.size()); 
	argNames[i] = s;
	defaultArgs[i] = dflt;
	inputFilters[i] = filter;
}

string InputRequester::getArgname(unsigned int i) const { 
	assert(i<argNames.size()); 
	return argNames[i];
}

string InputRequester::getDescription() {
	string s = name;
	if(argNames.size() > 0) {
		s += " (";
		for(unsigned int i=0; i<argNames.size(); i++) {
			if(inputFilters[i] && !getArgname(i).size()) {
				s += inputFilters[i]->getDescription();
			} else {
				s += getArgname(i);
				if(defaultArgs[i].size())
					s += " = "+defaultArgs[i];
			}
			if(i<argNames.size()-1)
				s += +", ";
		}
		s += + ")";
	}
	return s;
}


//---- NameSelector


bool nameselector_default_softmatch(const std::string& a, const std::string& b) {
	return lower(a) == lower(b.substr(0,a.size()));
}

string NameSelector::barf_control = "\033_BARF";
string NameSelector::exit_control = "\033_EXIT";

NameSelector::NameSelector(string t, string promptval, bool persist): InputRequester(t), catchAll(NULL), isPersistent(persist) {
	InputRequester::addArg(promptval);
	setSoftmatch();
}

void NameSelector::addChoice(string d, string nm, Selector_Option_Flags o, string mname, StreamInteractor* action) {
	if(!nm.size())
		nm = itos(choiceNames.size()+1);
	assert(nameMap.find(nm) == nameMap.end());
	if(!mname.size())
		mname = nm;
	nameMap.insert(std::make_pair(nm,choiceNames.size()));
	
	choiceNames.push_back(nm);
	choiceDescrips.push_back(d);
	choiceOut.push_back(mname);
	oflags.push_back(o);
	actions.push_back(action);
}

void NameSelector::addSynonym(string arg0, string syn) {
	map<std::string,unsigned int>::iterator it = nameMap.find(arg0);
	assert(it != nameMap.end());
	NameSelector::addChoice(choiceDescrips[it->second],syn,
			  Selector_Option_Flags(oflags[it->second] | SELECTOR_SYNONYM | SELECTOR_HIDDEN),
			  choiceOut[it->second],actions[it->second]);
}

void NameSelector::displayOptions() {
	printf("%s:\n---------------------------\n",name.c_str());
	for(unsigned int i=0; i<choiceNames.size(); i++)
		if(!(oflags[i] & SELECTOR_HIDDEN))
			printf("%s\t%s\n",choiceNames[i].c_str(),choiceDescrips[i].c_str());
}		

string NameSelector::getDescription() {
	string s = name;
	if(defaultArgs[0].size())
		s += " = "+defaultArgs[0];
	return s; 
}

void NameSelector::doIt() {
	bool forceBreak = false;
	assert(mystack && mydeque);
	do {
		
		// display options if choice not already made
		if(!mydeque->size()) {
			displayOptions();
			printf("---------------------------\n");
		}
		
		// input request loop
		while(1) {
			// prompt for input value onto stack
			InputRequester::doIt();
			
			// get argument off stack
			string myArg = popString();
			
			// break for control sequences
			if(startsWith(myArg, exit_control) || startsWith(myArg, barf_control) ) {
				if(startsWith(myArg, barf_control))
					mydeque->push_front(myArg);
				forceBreak = true;
				break;
			}
			
			// keep trying on empty input
			if(!myArg.size())
				continue;
			
			// find matching selection name
			map<std::string,unsigned int>::iterator it = nameMap.find(myArg);
			
			// if no direct match, apply soft matching when available
			if(it == nameMap.end() && softmatch) {
				
				// find soft matches
				vector<map<std::string,unsigned int>::iterator> matches;
				for(map<std::string,unsigned int>::iterator i = nameMap.begin(); i != nameMap.end(); i++)
					if(!(oflags[i->second] & SELECTOR_DISABLED) && softmatch(myArg,i->first))
						matches.push_back(i);
				
				// exactly one soft match is just right
				if(matches.size() == 1) {
					myArg = matches.back()->first;
					it = matches.back();
				}
				// too many ambiguous soft matches
				else if(matches.size() > 1) {
					printf("Error: ambiguous selection from:\n");
					for(unsigned int i=0; i<matches.size(); i++)
						printf("\t%s\n",matches[i]->first.c_str());
					continue;
				}
			}
			
			if(it == nameMap.end() || (oflags[it->second] & SELECTOR_DISABLED)) {
				if(catchAll) {
					mystack->push(myArg);
					catchAll->mystack = mystack;
					catchAll->mydeque = mydeque;
					catchAll->doIt();
					break;
				} else {
					printf("Error: unknown selection '%s'\n",myArg.c_str());
				}
			} else {
				if(actions[it->second]) {
					actions[it->second]->mystack = mystack;
					actions[it->second]->mydeque = mydeque;
					actions[it->second]->doIt();
				} else {
					mystack->push(choiceOut[it->second]);
				}
				break;
			}
		}
		
	} while(isPersistent && !forceBreak);
}


void OptionsMenu::addChoice(NamedInteractor* M, string nm, Selector_Option_Flags o) {
	NameSelector::addChoice(M->getDescription(),nm,o,"",M);
}


//----------------------------

void menutils_PrintQue(StreamInteractor* S) {
	printf("[ ");
	for(std::deque<string>::iterator it = S->mydeque->begin(); it != S->mydeque->end(); it++)
		printf("'%s' ",it->c_str());
	printf("]\n");
}

void menutils_PrintStack(StreamInteractor* S) {
	std::stack<string> tmp;
	while(S->mystack->size())
		tmp.push(S->popString());
	printf("[ ");
	while(tmp.size()) {
		S->mystack->push(S->popString());
		printf("'%s' ",S->mystack->top().c_str());
	}
	printf("]\n");
}

void menutils_StackSize(StreamInteractor* S) {
	S->mystack->push(itos(S->mystack->size()));
}

void menutils_Drop(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	S->mystack->pop();
}

void menutils_Dup(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	S->mystack->push(S->mystack->top());
}

void menutils_DropN(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	unsigned int n = S->popInt();
	if(!S->menutils_CheckStackSize(n))
		return;
	while(n--)
		S->mystack->pop();
}

void menutils_ClearStack(StreamInteractor* S) {
	while(!S->mystack->empty())
		S->mystack->pop();
}

void menutils_Swap(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(2))
		return;
	string a = S->popString();
	string b = S->popString();
	S->mystack->push(a);
	S->mystack->push(b);
}

void menutils_Rot(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	unsigned int n = S->popInt();
	if(!n || !S->menutils_CheckStackSize(n))
		return;
	for(unsigned int i=0; i<n-1; i++)
		S->mydeque->push_front(S->popString());
	string s = S->popString();
	for(unsigned int i=0; i<n-1; i++)
		S->mystack->push(S->popStringD());
	S->mystack->push(s);
}

void menutils_Select(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(3))
		return;
	string c = S->popString();
	string b = S->popString();
	string a = S->popString();
	if(c == "true" || atof(c.c_str()))
		S->mystack->push(a);
	else
		S->mystack->push(b);
}

void menutils_Exec(StreamInteractor* S) {
	if(!S->menutils_CheckStackSize(1))
		return;
	vector<string> v = split(S->popString());
	while(v.size()) {
		S->mydeque->push_front(v.back());
		v.pop_back();
	}
}

void menutils_Barf(StreamInteractor* S) {
	S->mydeque->push_front(NameSelector::barf_control);
}

void menutils_Exit(StreamInteractor* S) {
	S->mydeque->push_front(NameSelector::exit_control);
}
