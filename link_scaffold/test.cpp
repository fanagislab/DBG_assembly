#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>


using namespace std;

inline int ctgStr2Id(string ctg_id_str)
{	return atoi( (ctg_id_str.substr(4,ctg_id_str.size()-4)).c_str() );
}

int get_pair_id( int id)
{   return ( id % 2 ==0 ) ? (id - 1) : (id + 1);
}

int main(int argc, char *argv[])
{	
	vector<int> scaff_ids;
	scaff_ids.push_back(-201);
	
	cout << scaff_ids[0] << endl;
	
}

