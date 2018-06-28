#define TreeMakerTreeClass_cxx
#include "TreeMaker/TreeMaker/interface/TreeMakerTreeClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TreeMakerTreeClass::FillFriendTree_JERC() {
	if (fChain == 0 || fFriend == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    	Long64_t ientry = LoadTree(jentry);
    	if (ientry < 0) break;
    	fChain->GetEntry(jentry);
      	if (Cut(ientry) < 0) continue;

		// get all variable values
		GenerateNewIndex();
    	for(auto & variable : variables){
    		if(variable->GetNameInTree().find("Factor")!=std::string::npos || variable->GetNameInTree().find("jecUnc")!=std::string::npos) {
    			if(variable->GetNameInTree()=="Jets_jecFactor")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*Jets_jecFactor);
	    		else if(variable->GetNameInTree()=="Jets_jecUnc")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*Jets_jecUnc);
	    		else if(variable->GetNameInTree()=="Jets_jerFactor")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*Jets_jerFactor);
	    		else if(variable->GetNameInTree()=="Jets_jerFactorUp")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*Jets_jerFactorUp);
	    		else if(variable->GetNameInTree()=="Jets_jerFactorDown")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*Jets_jerFactorDown);
	    		else if(variable->GetNameInTree()=="JetsJECdown_jerFactor")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*JetsJECdown_jerFactor);
	    		else if(variable->GetNameInTree()=="JetsJECup_jerFactor")
	    			dynamic_cast<TreeObject<std::vector<double> > *>(variable)->SetValue(*JetsJECup_jerFactor);
	    	}
	        else if(variable->GetNameInTree().find("origIndex")!=std::string::npos) {
	        	if(variable->GetNameInTree()=="Jets_origIndex")
		        	dynamic_cast<TreeObject<std::vector<int> > *>(variable)->SetValue(*Jets_origIndex);
		        else if(variable->GetNameInTree()=="JetsJECdown_origIndex")
		        	dynamic_cast<TreeObject<std::vector<int> > *>(variable)->SetValue(*JetsJECdown_origIndex);
		        else if(variable->GetNameInTree()=="JetsJECup_origIndex")
		        	dynamic_cast<TreeObject<std::vector<int> > *>(variable)->SetValue(*JetsJECup_origIndex);
		        else if(variable->GetNameInTree()=="JetsJERdown_origIndex")
		        	dynamic_cast<TreeObject<std::vector<int> > *>(variable)->SetValue(*JetsJERdown_origIndex);
		        else if(variable->GetNameInTree()=="JetsJERup_origIndex")
		        	dynamic_cast<TreeObject<std::vector<int> > *>(variable)->SetValue(*JetsJERup_origIndex);
		    }
	        else {
	        	if(variable->GetNameInTree()=="GenJets")
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*GenJets);
		        else if(variable->GetNameInTree()=="Jets")
		        	//dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*Jets);
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(RemakeJets());
		        else if(variable->GetNameInTree()=="JetsJECdown")
		        	//dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*JetsJECdown);
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(RemakeJetsJECdown());
		        else if(variable->GetNameInTree()=="JetsJECup")
		        	//dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*JetsJECup);
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(RemakeJetsJECup());
		        else if(variable->GetNameInTree()=="JetsJERdown")
		        	//dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*JetsJERdown);
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(RemakeJetsJERdown());
		        else if(variable->GetNameInTree()=="JetsJERup")
		        	//dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(*JetsJERup);
		        	dynamic_cast<TreeObject<std::vector<TLorentzVector> > *>(variable)->SetValue(RemakeJetsJERup());
	        }
    	}
      	fFriend->Fill();
    }
}

void TreeMakerTreeClass::GenerateNewIndex() {
	newIndex.reserve(Jets_origIndex->size());
	for(unsigned j = 0; j < Jets_origIndex->size(); ++j){
    	//reverse the index vector
    	newIndex[Jets_origIndex->at(j)] = j;
	}
}

void TreeMakerTreeClass::GetActiveBranches() {
	//Prototype: make_pair("<branch_name>","<draw_command>[:<additional_draw_commands>]")
	//Example: make_pair("Jets","Jets.Pt():Jets.Eta()"),
	branchesVTLorentzVector = {
		//make_pair("GenJets",     "GenJets.Pt()"),
		make_pair("Jets",        "Jets.Pt()"),
		make_pair("JetsJECdown", "JetsJECdown.Pt()"),
		make_pair("JetsJECup",   "JetsJECup.Pt()"),
		make_pair("JetsJERdown", "JetsJERdown.Pt()"),
		make_pair("JetsJERup",   "JetsJERup.Pt()")
	};
	branchesVInt = {
		make_pair("Jets_origIndex",        ""),
		make_pair("JetsJECdown_origIndex", ""),
		make_pair("JetsJECup_origIndex",   ""),
		make_pair("JetsJERdown_origIndex", ""),
		make_pair("JetsJERup_origIndex",   "")
	};
	branchesVDouble = {
		make_pair("Jets_jecFactor",        ""),
		make_pair("Jets_jecUnc",           ""),
		make_pair("Jets_jerFactor",        ""),
		make_pair("Jets_jerFactorUp",      ""),
		make_pair("Jets_jerFactorDown",    ""),
		make_pair("JetsJECdown_jerFactor", ""),
		make_pair("JetsJECup_jerFactor",   "")
	};
	branchesCombined.reserve( branchesVTLorentzVector.size() + branchesVInt.size() + branchesVDouble.size() ); // preallocate memory
	branchesCombined.insert( branchesCombined.end(), branchesVTLorentzVector.begin(), branchesVTLorentzVector.end() );
	branchesCombined.insert( branchesCombined.end(), branchesVInt.begin(), branchesVInt.end() );
	branchesCombined.insert( branchesCombined.end(), branchesVDouble.begin(), branchesVDouble.end() );
}

void TreeMakerTreeClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TreeMakerTreeClass.C
//      root> TreeMakerTreeClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    	Long64_t ientry = LoadTree(jentry);
    	if (ientry < 0) break;
    	nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    }
}

void TreeMakerTreeClass::MakeEmptyFriendTree(std::string name) {
	VarTypeNames = {
		"VectorInt","VectorDouble","VectorTLorentzVector",
    };
    VarTypes = {
        t_vint,t_vdouble,t_vlorentz
    };

	std::set<std::string>          nameSet;
    std::stringstream              skipMessage;
    std::stringstream              message;
	for(unsigned v = 0; v < VarTypeNames.size(); ++v){
		VBranchType VarNames;
		switch(VarTypes[v]) {
			case  TreeTypes::t_vint     : VarNames = branchesVInt; break;
			case  TreeTypes::t_vdouble  : VarNames = branchesVDouble; break;
			case  TreeTypes::t_vlorentz : VarNames = branchesVTLorentzVector; break;
			default                     : continue;
		}
        message << VarTypeNames.at(v) << ":" << "\n";
        for(const auto & VarNamePair : VarNames){
        	const auto & VarName = VarNamePair.first;
            //check for an exact repeat of an existing name
            if(nameSet.find(VarName)!=nameSet.end()){
                skipMessage << VarName << "\n";
                continue;
            }
            else nameSet.emplace(VarName);
            //check for the right type
            TreeObjectBase* tmp = nullptr;
            switch(VarTypes[v]){
                case TreeTypes::t_vint     : tmp = new TreeObject<std::vector<int> >(VarName); break;
                case TreeTypes::t_vdouble  : tmp = new TreeObject<std::vector<double> >(VarName); break;
                case TreeTypes::t_vlorentz : tmp = new TreeObject<std::vector<TLorentzVector> >(VarName); break;
                default                    : continue;
            }
            //if a known type was found, initialize and store the object
            if(tmp) {
            	tmp->Initialize(nameCache,message);
                variables.push_back(tmp);
            }
        }
    }
    //print info
    if(!skipMessage.str().empty()) std::cout << "MakeEmptyFriendTree: Skipping repeated branches:\n" << skipMessage.str();
    std::cout << "MakeEmptyFriendTree:\n" << message.str();

    TFile* friendFile = TFile::Open("friendFile.root","RECREATE");
	fFriend = new TTree(name.c_str(), name.c_str());
	fFriend->SetAutoSave(10000000000);
    fFriend->SetAutoFlush(1000000);
    //add branches to tree
    for(auto & variable : variables){
        variable->SetTree(fFriend);
        variable->AddBranch();
    }
}

std::pair<std::string,int> TreeMakerTreeClass::MakeScanString(VBranchType branches, std::string friendTreeName, bool compare, bool test) {
	std::string varexp = "";
	int items = 0;

	for(auto branch : branches) {
		if(branch.second.empty()) continue;

		std::vector<std::string> vars;
		if(branch.second.find(":")!=std::string::npos) {
			std::string remainder = branch.second;
			while(!remainder.empty()) {
				int pos = (remainder.find(":")==std::string::npos) ? -1 : remainder.find(":");
				if(pos==-1) {
					vars.push_back(remainder);
					remainder = "";
				}
				else {
					vars.push_back(remainder.substr(0,pos));
					remainder = remainder.substr(pos+1,remainder.size()-pos);
				}
			}
		}
		else {
			vars.push_back(branch.second);
		}

		for(auto var : vars) {
			items++;

			// Put in joiner if item not blank (by construction) and it isn't the first item (varexp isn't empty)
			// This specific condition is needed for sublists of variables
			if(!varexp.empty())
				varexp += (test) ? "+" : ":";

			if(!friendTreeName.empty()) {
				if(compare) {
					varexp += (var + "/" + friendTreeName + "." + var);
				}
				else {
					varexp += (var + ":" + friendTreeName + "." + var);
				}
			}
			else {
				varexp += var;
			}	
		}

		// Last entry condition
		//if(branch.first != branches.back().first)
		//	varexp += (test) ? "+" : ":";
	}
	varexp = (test) ? "Sum$(" + varexp + ")" : varexp;
	return make_pair(varexp,items);
}

// Jets(JEC,JER)
VTLorentzVector TreeMakerTreeClass::RemakeJets() {
	VTLorentzVector JetsFriend(Jets->size());
	for(unsigned j = 0; j < Jets_origIndex->size(); ++j){
    	int i = newIndex[Jets_origIndex->at(j)];
    	JetsFriend[j] = Jets->at(i)*(1./(Jets_jecFactor->at(i)*Jets_jerFactor->at(i)))*Jets_jecFactor->at(j)*Jets_jerFactor->at(j);
	}
	return JetsFriend;
}

// Jets(JECup,JER)
VTLorentzVector TreeMakerTreeClass::RemakeJetsJECup() {
	VTLorentzVector JetsJECupFriend(Jets->size());
	for(unsigned j = 0; j < JetsJECup_origIndex->size(); ++j){
    	//JetsJECup_origIndex is sorted in the final order after JEC uncertainty variation
    	//go up to common ancestor, then down to central smeared collection
    	int i = newIndex[JetsJECup_origIndex->at(j)];
    	//undo central smearing, apply JEC unc, redo smearing w/ new smearing factor
    	JetsJECupFriend[j] = Jets->at(i)*(1./Jets_jerFactor->at(i))*(1+Jets_jecUnc->at(i))*JetsJECup_jerFactor->at(j);
	}
	return JetsJECupFriend;
}

// Jets(JECdown,JER)
VTLorentzVector TreeMakerTreeClass::RemakeJetsJECdown() {
	VTLorentzVector JetsJECdownFriend(Jets->size());
	for(unsigned j = 0; j < JetsJECdown_origIndex->size(); ++j){
    	int i = newIndex[JetsJECdown_origIndex->at(j)];
    	JetsJECdownFriend[j] = Jets->at(i)*(1./Jets_jerFactor->at(i))*(1-Jets_jecUnc->at(i))*JetsJECdown_jerFactor->at(j);
	}
	return JetsJECdownFriend;
}

// Jets(JEC,JERup)
VTLorentzVector TreeMakerTreeClass::RemakeJetsJERup() {
	VTLorentzVector JetsJERupFriend(Jets->size());
	for(unsigned j = 0; j < JetsJERup_origIndex->size(); ++j){
    	int i = newIndex[JetsJERup_origIndex->at(j)];
    	JetsJERupFriend[j] = Jets->at(i)*(1./Jets_jerFactor->at(i))*Jets_jerFactorUp->at(i);
	}
	return JetsJERupFriend;
}

// Jets(JEC,JERdown)
VTLorentzVector TreeMakerTreeClass::RemakeJetsJERdown() {
	VTLorentzVector JetsJERdownFriend(Jets->size());
	for(unsigned j = 0; j < JetsJERdown_origIndex->size(); ++j){
    	int i = newIndex[JetsJERdown_origIndex->at(j)];
    	JetsJERdownFriend[j] = Jets->at(i)*(1./Jets_jerFactor->at(i))*Jets_jerFactorDown->at(i);
	}
	return JetsJERdownFriend;
}

/*
Run using:
.L ../../TreeMaker/src/TreeMakerTreeClass.C+
TreeMakerTreeClass tmtc
tmtc.TestJERCBranches()

Replicating to some extent:
PreSelection->Scan("GenJets.Pt():Jets.Pt():Jets.Eta():Jets_origIndex:Jets_jecFactor:Jets_jecUnc:Jets_jerFactor:Jets_jerFactorUp:Jets_jerFactorDown:JetsJECdown.Pt():JetsJECdown_jerFactor:JetsJECdown_origIndex:JetsJECup.Pt():JetsJECup_jerFactor:JetsJECup_origIndex:JetsJERdown.Pt():JetsJERdown_origIndex:JetsJERup.Pt():JetsJERup_origIndex","","colsize=22")
*/
void TreeMakerTreeClass::TestJERCBranches(int scanDepth, bool verbose) {
    GetActiveBranches();
    fChain->SetBranchStatus("*",0);
    for(auto branch : branchesCombined) {
    	fChain->SetBranchStatus(branch.first.c_str(),1);
    }


    MakeEmptyFriendTree("JERCFriend");
    FillFriendTree_JERC();
    if(verbose) {
    	std::cout << "Original tree structure:" << std::endl;
    	fChain->Show(0);
    	std::cout << "Friend tree structure:" << std::endl;
    	fFriend->Show(0);
	}
    fChain->AddFriend(fFriend);

    std::string scanString = MakeScanString(branchesCombined,"JERCFriend",true,false).first;
    std::cout << "\nScan varexp: " << scanString << std::endl;
    fChain->SetScanField(scanDepth); 
    fChain->Scan(scanString.c_str(),"","colsize=10 precision=6");

    std::cout << "\nTestJERCBranches::Testing all entries for inconsistencies ... " << std::endl;
    std::pair<std::string,int> cutStringPair = MakeScanString(branchesCombined,"JERCFriend",true,true);
    scanString = cutStringPair.first;
    std::stringstream cutString; cutString << "abs(" << scanString << " - (" << cutStringPair.first.replace(0,3,"Length") << "*" << cutStringPair.second << ")) > 0.0001";
    std::cout << "Scan varexp: " << scanString << std::endl;
    std::cout << "Scan cut: " << cutString.str() << std::endl;
    fChain->SetScanField(0);
	Long64_t inconsistencies = fChain->Scan((scanString+":"+cutString.str()).c_str(),cutString.str().c_str(),"colsize=22 precision=6");
	if(inconsistencies==0) {
		std::cout << "\nAll values in the friend tree match those in the original tree!" << std::endl;
	}
	else {
		std::cout << "\nThere were " << inconsistencies << " events where one or more values in the friend tree didn't match the original tree ;(" << std::endl;
	}
}
