#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace std;

int main() {
	vector<float> single_logFC;	//log2 fold change for sgRNA+sgCon
	string entry;
	int gene_count=77;	
//	int num_seq=0;
		
	ifstream singleZ;
	singleZ.open("Z-single.txt");
	
	for (int x=0; x<gene_count;x++)
	{
		getline(singleZ,entry);
		float FC=stof(entry);
		single_logFC.push_back(FC);
	}	
	singleZ.close();


	ifstream fold_changes;
	fold_changes.open("Z-avg.txt");	
	string output="";
	double FC;
	//Reads barcode list
	if (fold_changes.is_open())
	{
		int a=0;
		int b=0;
		
		while (getline(fold_changes,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			stringstream stream(entry);
			
			stream>>a;
			stream>>b;
			output+=entry+"\t"+to_string(single_logFC[a]+single_logFC[b])+"\n";
		}
	}
	fold_changes.close();
	
	ofstream results;
	results.open("Z-observed-expected.txt");
	results<<"gene A\tgene B\tobserved Z\texpected Z\n";
	results<<output;
	results.close();

	return 0;
}
