#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int main() {
	vector<int> count;	//sgRNA count list
	string entry;
	int num_seq=0;
		
	ifstream count_read;
	count_read.open("library_count.txt");
	
	//Reads barcode list
	if (count_read.is_open())
	{
		while (getline(count_read,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			int cnt=stoi(entry);	//subtract pseudocount
			count.push_back(cnt);
			num_seq++;
		}
	}
	count_read.close();
	
	int sgRNA_count=sqrt(num_seq);	//infers number of sgRNAs by square root of combinations
		
	ofstream perm_sum;
	perm_sum.open("permutation_added.txt");
		
	for (unsigned int x=0; x<sgRNA_count;x++)
	{
		for (unsigned int y=x; y<sgRNA_count;y++)
		{
			if (x==y)
				perm_sum<<x<<"\t"<<y<<"\t"<<count[y+x*sgRNA_count]<<"\n";
			else		
				perm_sum<<x<<"\t"<<y<<"\t"<<count[y+x*sgRNA_count]+count[x+y*sgRNA_count]-1<<"\n";	//subtract one duplicate pseudocount
		}
	}
	perm_sum.close();
	
	return 0;
}
