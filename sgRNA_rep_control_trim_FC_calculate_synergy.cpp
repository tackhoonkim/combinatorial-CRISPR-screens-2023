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
	int sgRNA_count=233;	
//	int num_seq=0;
	for (int x=0; x<sgRNA_count*sgRNA_count;x++)
		logFC.push_back(0);
		
	ifstream fold_changes;
	fold_changes.open("FC-perm added.txt");
	
	//Reads barcode list
	if (fold_changes.is_open())
	{
		int a=0;
		int b=0;
		while (getline(fold_changes,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			float FC=stof(entry);
			logFC[a+b*sgRNA_count]+=FC;
			a++;
			if (a>=sgRNA_count)
			{
				b++;
				a=b;
			}
		}
		for (unsigned int x=0; x<sgRNA_count;x++)
		{
			for (unsigned int y=x; y<sgRNA_count;y++)
				logFC[x+y*sgRNA_count]=logFC[y+x*sgRNA_count];
		}
	}
	fold_changes.close();
	
	int gene_count=77;	//hard coding gene count (76 genes, 3 controls)
	
	ofstream sgRNA_rep_control;
	sgRNA_rep_control.open("sgRNA_rep_control.txt");
		
	for (unsigned int x=0; x<gene_count;x++)
	{
		for (unsigned int y=x; y<gene_count;y++)
			sgRNA_rep_control<<x<<"\t"<<y<<"\t"<<logFC[3*y+3*x*sgRNA_count]<<"\t"			//sgRNA #1-sgRNA #1
										<<logFC[3*y+1+(3*x+1)*sgRNA_count]<<"\t"			//sgRNA #2-sgRNA #2
										<<logFC[3*y+2+(3*x+2)*sgRNA_count]<<"\n";			//sgRNA #3-sgRNA #3
	}
	sgRNA_rep_control.close();

	
	vector<float> single_sgRNA_logFC;	//sgRNA+sgCon average calculating.
	for (unsigned int x=0;x<sgRNA_count;x++)
	{
		single_sgRNA_logFC.push_back(0);
		int count=0;
		for (int y=228; y<sgRNA_count;y++)
		{
			if (logFC[sgRNA_count*x+y]<100)		// if the count was not excluded because of low coverage: excluded ones were given pseudo-FC of 100
			{
				single_sgRNA_logFC[x]+=logFC[sgRNA_count*x+y];
				count++;
			}
		}
		single_sgRNA_logFC[x]/=count;
	}

	ofstream expected;
	expected.open("expected FC-trim.txt");
	for (unsigned int x=0; x<sgRNA_count;x++)
	{
		for (unsigned int y=x; y<sgRNA_count;y++)
			expected<<x<<"\t"<<y<<"\t"<<single_sgRNA_logFC[x]+single_sgRNA_logFC[y]<<"\n";
	}
	expected.close();



	ofstream synergy;
	synergy.open("synergy score-trim.txt");
	for (unsigned int x=0; x<sgRNA_count;x++)
	{
		for (unsigned int y=x; y<sgRNA_count;y++)
		{
			if (logFC[y+x*sgRNA_count]==100)
				synergy<<x<<"\t"<<y<<"\t"<<"NA"<<"\n";
			else
				synergy<<x<<"\t"<<y<<"\t"<<logFC[y+x*sgRNA_count]-single_sgRNA_logFC[x]-single_sgRNA_logFC[y]<<"\n";
		}
	
	}
	synergy.close();
	
		
	return 0;
}
