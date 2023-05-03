#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

int main() {
	vector<string> names;	//gene names for sgRNA numbers
	string entry;
		
	ifstream gene_names;
	gene_names.open("gene names.txt");
	
	//Reads barcode list
	if (gene_names.is_open())
	{
		while (getline(gene_names,entry))
		{
			entry.erase(entry.find_last_not_of(" \r\n\t")+1);	//trims new line
			names.push_back(entry);
		}
	}
	gene_names.close();

	ifstream numbers1;
	numbers1.open("sgRNA numbers1.txt");
	ifstream numbers2;
	numbers2.open("sgRNA numbers2.txt");
	string entry1;
	string entry2;

	ofstream output;
	output.open("gene names-converted.txt");
	if (numbers1.is_open() && numbers2.is_open())
	{
		while (getline(numbers1,entry1))
		{
			getline(numbers2,entry2);	
			entry1.erase(entry1.find_last_not_of(" \r\n\t")+1);	//trims new line
			entry2.erase(entry2.find_last_not_of(" \r\n\t")+1);	//trims new line
			output<<names[stoi(entry1)]<<","<<names[stoi(entry2)]<<"\n";
		}
	}
	numbers1.close();
	numbers2.close();
	output.close();

	return 0;
}
