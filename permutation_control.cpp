#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int main() {
	vector<float> logFC;	//log2 fold change list
	string entry;
	int num_seq=0;
		
	ifstream fold_changes;
	fold_changes.open("FC-raw.txt");
	
	//Reads barcode list
	if (fold_changes.is_open())
	{
		while (getline(fold_changes,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			if (entry=="NA")
				logFC.push_back(100.0);
			else
			{
			float FC=stof(entry);
			logFC.push_back(FC);
			}
			num_seq++;
		}
	}
	fold_changes.close();
	
	int sgRNA_count=sqrt(num_seq);	//infers number of sgRNAs by square root of combinations
		
	ofstream perm_control;
	perm_control.open("permutation_control.txt");
		
	for (unsigned int x=0; x<sgRNA_count;x++)
	{
		for (unsigned int y=x; y<sgRNA_count;y++)
			perm_control<<x<<"\t"<<y<<"\t"<<logFC[y+x*sgRNA_count]<<"\t"<<logFC[x+y*sgRNA_count]<<"\n";
	}
	perm_control.close();
	
	return 0;
}
