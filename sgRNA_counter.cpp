#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>


using namespace std;

int match(string read, vector<string> dict)
{
	//read: NGS read, dict: list of sgRNAs, key: constant seq for finding sgRNA location
	//match length: the length of sgRNA to do alignment
	read.erase(read.find_last_not_of(" \r\n\t")+1);
			
	for (int i=0; i<dict.size();i++)
	{
		if (read==dict[i])
			return i;
		
	}
	return -1;
}

int main() {
	vector<string> dictionary;	//sgRNA library ref seq
	string entry;
	string read_sequence1;
	string read_sequence2;
	int num_seq=0;
	
	int num_reads=0;
	int perfect_matches=0;
	int non_perfect_matches=0;
	
	ifstream sgRNAs;
	sgRNAs.open("sgRNAs-rev compl.txt");
	
	//Reads barcode list
	if (sgRNAs.is_open())
	{
		while (getline(sgRNAs,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			dictionary.push_back(entry);
			num_seq++;
		}
	}
	//declares barcode count matrix
	int counter_matrix[num_seq*num_seq];
	for (unsigned i=0; i<num_seq*num_seq;i++)
		counter_matrix[i]=1;	//pseudocount to enable logarithmic transformation of dropouts
	
	//declares positive control (sgRPL) count matrix
	int poscon[num_seq];
	for (unsigned int i=0; i<num_seq;i++)
		poscon[i]=1;
	
	//Reads Fastq file
	ifstream fastq1;
	fastq1.open("NGS-1.fastq");
	ifstream fastq2;
	fastq2.open("NGS-2.fastq");
	int num=0;
	
	while (getline(fastq1,read_sequence1))
	{
		getline(fastq2,read_sequence2);
		num++;
		if(num%4==2){
			num_reads++;
			int sg1=match(read_sequence1,dictionary);
			int sg2=match(read_sequence2,dictionary);
			if (sg1>=0 && sg2>=0)
			{
				counter_matrix[sg2+sg1*num_seq]++;
				perfect_matches++;
			}else
				non_perfect_matches++;
		}
	}
	fastq1.close();
	fastq2.close();
	
	int guides_with_reads=0;
	int guides_no_reads=0;
	ofstream count_table;
	count_table.open("library_count.txt");
	for (unsigned int x=0; x<num_seq;x++)
	{
		for (unsigned int y=0; y<num_seq;y++)
		{
			count_table<<x<<"\t"<<y<<"\t"<<dictionary[x]<<"\t"<<dictionary[y]<<"\t"<<counter_matrix[y+x*num_seq]<<"\n";
			if (counter_matrix[y+x*num_seq]>0){
				guides_with_reads++;
			}else{
				guides_no_reads++;
			}
		}
	}
	count_table.close();
	
	
	float percent_matched=round((float)perfect_matches/(perfect_matches+non_perfect_matches)*1000)/10.0;
	float percent_no_reads=round((float)guides_no_reads/(guides_no_reads+guides_with_reads)*1000)/10.0;
	
	sort(counter_matrix,(counter_matrix+num_seq*num_seq));
	
	int top_10=counter_matrix[(int)round((float)num_seq*num_seq*9/10)];
	int bottom_10=counter_matrix[(int)round((float)num_seq*num_seq/10)];
	
	float skew_ratio=-1.0;
	if (top_10*bottom_10>0){
		skew_ratio=(float)top_10/bottom_10;
	}
	
	ofstream summary;
	summary.open("statistics.txt");
	summary<<"Number of perfect guide matches: "<<perfect_matches<<"\n";
	summary<<"Number of nonperfect guide matches: "<<non_perfect_matches<<"\n";
	summary<<"Number of reads processed: "<<num_reads<<"\n";
	summary<<"Percenteage of guides that matched perfectly: "<<percent_matched<<"\n";
	summary<<"Percentage of undetected guides: "<<percent_no_reads<<"\n";
	summary<<"Skew ratio of top 10% to bottom 10%: "<<skew_ratio<<"\n";
	summary.close();

	return 0;
}
