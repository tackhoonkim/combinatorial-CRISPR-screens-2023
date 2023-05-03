#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace std;

int main() {
	int gene_num=77;	//77 genes- tyrosine kinases
	float rep1scores[gene_num*gene_num];	//list of Z scores for D20/D0 rep1
	float rep2scores[gene_num*gene_num]; //list of Z scores for D20/D0 rep2
	int reps[gene_num*gene_num];	//number of sgRNAs per gene
	for (int i=0; i<gene_num*gene_num; i++)
	{
		rep1scores[i]=0.0;
		rep2scores[i]=0.0;
		reps[i]=0;
	}
		
	ifstream input;
	input.open("Z.txt");	//Z score list

	string entry;
	int x1, x2;
	float y;
	while (getline(input,entry))	//parse gene number, Z score
	{
		stringstream stream(entry);
		stream>>x1;
		stream>>x2;
		reps[x1*gene_num+x2]++;
		stream>>y;
		rep1scores[x1*gene_num+x2]+=y;
		stream>>y;
		rep2scores[x1*gene_num+x2]+=y;
	}
	input.close();
		
	for (int i=0; i<gene_num; i++)
	{
		for (int j=i; j<gene_num; j++)
		{
			rep1scores[i*gene_num+j]/=reps[i*gene_num+j];
			rep2scores[i*gene_num+j]/=reps[i*gene_num+j];	//averaging the Z scores
		}
	}

	ofstream output;
	output.open("Z-avg.txt");
	output<<"gene A\tgene B\trep1\trep2\n";
	for (int i=0; i<gene_num; i++)
	{
		for (int j=i; j<gene_num; j++)
			output<<i<<"\t"<<j<<"\t"<<rep1scores[i*gene_num+j]<<"\t"
									<<rep2scores[i*gene_num+j]<<"\n";
	}
	output.close();

	return 0;
}
