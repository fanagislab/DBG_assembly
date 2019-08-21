//author: Fan Wei, email: fanw@genomics.org.cn, date: 2011-1-20

//注意：本文件内的函数可以供外界使用，而变量则只能供本文件内函数使用
//注意：内联函数要慎用，经常出错，仅特别简单的函数可以作为内联函数，务必逐个测试通过

//关于坐标系：本程序采用从1开始的坐标系来记录位置，涉及continuous region start and end
//position, 需要修改的低频区的start和end位置。

#include "seqKmer.h"
#include "global.h"
#include "correct.h"
#include "gzstream.h"


//get exactly continuous low or high frequence kmers on reads, this is a cluster function
void get_cont_kmerfreq_region(string &read, uint8_t *freq, vector <FreqContReg> &continuous_regions)
{	
	int total_kmers = read.size() - KmerSize + 1;
	int i = 0;
	while (i < total_kmers)
	{	
		// pass one low kmer frequence region
		int low_start = i;
		while ( i < total_kmers )
		{	string kseq = read.substr(i,KmerSize);
			uint64_t kbit = seq2bit(kseq);
			int kbit_freq = get_freq(freq, kbit);
		
			if (kbit_freq == 0)  //低频kmer的freq已被设置为0
 			{	i++;
			}else{
				break;
			}
		} 
		int low_end = i;

		if (low_end > low_start)
		{	FreqContReg element;
			element.start = low_start + 1;
			element.end = low_end;
			element.status = 0;
			continuous_regions.push_back(element);
		}
		
		// get one high kmer frequence region
		int high_start = i;
		while ( i < total_kmers )
		{	string kseq = read.substr(i,KmerSize);
			uint64_t kbit = seq2bit(kseq);
			int kbit_freq = get_freq(freq, kbit);
			
			if (kbit_freq == 1)  //高频kmer的freq已被设置为1
			{	i++;
			}else{
				break;
			}
		}
		int high_end = i;

		if (high_end > high_start)
		{	FreqContReg element;
			element.start =  high_start + 1;
			element.end = high_end;
			element.status = 1;
			continuous_regions.push_back(element);
		}
	}

}


//The fast error correction method, 找到一种方案后就停止继续搜索。
//因为此处对应的错误理论上应发生在unique genomic region, 所以只有一种正确的选择
int correct_one_base(string &read, uint8_t *freq,  int error_pos, int check_start, int check_end)
{	int is_corrected = 0;
	char error_base = read[error_pos-1];
	int check_num = check_end - check_start + 1;
	char correct_base = error_base;

	for (int i=0; i<4; i++)
	{	if (error_base != bases[i])
		{	read[error_pos-1] = bases[i];
			int high_kmer_num = 0;
			for (int j=check_start-1; j<check_end; j++)
			{	string kseq = read.substr(j,KmerSize);
				uint64_t kbit = seq2bit(kseq);
				int kbit_freq = get_freq(freq, kbit);
				if (kbit_freq == 1)
				{	high_kmer_num ++;
				}else{
					break;
				}
			}
			if (high_kmer_num == check_num)
			{	correct_base = bases[i];
				break; //第一次找到后便结束工作
			}
		}
	}

	if (correct_base != error_base)
	{	is_corrected = 1;
	}else{
		read[error_pos-1] = error_base;
	}
	return is_corrected;
}



//连接相邻的高频区域，同时去除小的高频区域(length < HighFreqRegLenCutoff)
void get_high_freq_region(vector <FreqContReg> &freqContRegs, vector <FreqContReg> &highFreqContRegs, int &num_in_highFreqRegs, int &num_out_highFreqRegs)
{	
	num_in_highFreqRegs = 0;
	num_out_highFreqRegs = 0;
	int total_blocks = freqContRegs.size();
	int i = 0;

	while (i < total_blocks)
	{
		while ( i < total_blocks && freqContRegs[i].status == 0 )
		{       i++;
		} // pass the low_freq regions
		
		int highFreq_start = i;
		while ( i < total_blocks && freqContRegs[i].status == 1 )
		{       i++;
				num_in_highFreqRegs ++;
		} // get one continuous high freq region
		int highFreq_end = i;   

		if (highFreq_end > highFreq_start && freqContRegs[highFreq_end-1].end - freqContRegs[highFreq_start].start + 1 >= HighFreqRegLenCutoff)
		{       FreqContReg element;
				element.start = freqContRegs[highFreq_start].start;
				element.end = freqContRegs[highFreq_end-1].end;
				element.status = 1;
				highFreqContRegs.push_back( element );
				num_out_highFreqRegs ++;
		}
	}
	
}


//To correct the error base, modify on the read directly
void correct_one_read(string &read_head, string &read, uint8_t *freq, int &correct_one_base_score, int &correct_multi_base_score, int &is_deleted, int &trim_left_end_len, int &trim_right_end_len)
{	
	
	vector <FreqContReg> freqContRegs;
	int read_len =  read.size();
	
	int accum_change = 0;
	is_deleted = 0;
	correct_one_base_score = 0;
	correct_multi_base_score = 0;
	
	//这两个变量只用在对尾巴进行纠错
	int right_last_change_pos = read_len + 1;
	trim_right_end_len = 0;

	//这两个变量只用在对头部进行纠错
	int left_last_change_pos = 0;
	trim_left_end_len = 0;

	//在整条reads范围对连续的高低频进行分区，获得每一个高频区域和低频区域的起始和终止坐标
	get_cont_kmerfreq_region(read,freq,freqContRegs);


	//简单纠错方法：仅用于纠正两个高频区域中间的低频区域
	int num_cRegs = freqContRegs.size();
	for (int i=1; i<num_cRegs-1; i++)
	{	
		if (freqContRegs[i].status != 0) {continue;}
		
		//判断：如果没有可改的量了，则停止继续往下纠错
		if (accum_change >= Max_change_in_one_read)
		{	break;
		}

		int region_size = freqContRegs[i].end - freqContRegs[i].start + 1;
		int is_corrected = 0;

		if (i > 0 && i < num_cRegs - 1 && region_size == KmerSize)
		{	is_corrected = correct_one_base(read, freq, freqContRegs[i].end, freqContRegs[i].start, freqContRegs[i].end);
		}

		if (is_corrected)
		{	correct_one_base_score ++;
			freqContRegs[i].status = 1;
			accum_change ++;
		}
	}
	
	//成功纠错后相应的低频区域变为高频区域，连接相邻的高频区域，去除小的高频区域(length < HighFreqRegLenCutoff)，获得满足长度要求的高频区域
	vector <FreqContReg> highFreqContRegs;
	int num_in_highRegs = 0;
	int num_out_highRegs = 0;
	get_high_freq_region(freqContRegs, highFreqContRegs, num_in_highRegs, num_out_highRegs);
	int num_hRegs = highFreqContRegs.size(); 
	
	//修理大高频区的边缘(因为高频区也可能含有测序错误)，两边各削去1/3个HighFreqRegLenCutoff长度，这样做将使高频区的边界更加准确，使后面的纠错更容易成功
	int edge_cut_len = HighFreqRegLenCutoff / 3;
	int Kmer_num_on_read = read_len - KmerSize + 1;
	for (int i=0; i<highFreqContRegs.size(); i++)
	{	if (highFreqContRegs[i].start != 1)
		{	highFreqContRegs[i].start += edge_cut_len;
		}
		if (highFreqContRegs[i].end != Kmer_num_on_read)
		{	highFreqContRegs[i].end -= edge_cut_len;
		}
	}

	//如果没有找到任何高频区，则退出
	if (num_hRegs == 0)
	{	is_deleted = 1;
		return;
	}

	//分支定界方法：修改两个高频区域中间的错误base
	//因为高频区的边界可能不准确，也就是包含测序错误，所以有的时候，从左边不能纠对，却能从右边纠对，故当一个方向失败时尝试另一个方向纠错
	vector<int> failCorrectReg_ids;
	if (num_hRegs >= 2)
	{	for (int i=0; i<num_hRegs-1; i++)
		{	
			if (accum_change >= Max_change_in_one_read)
			{	//把后面的region都存入failCorrectReg_ids
				for (int k=i; k<num_hRegs-1; k++)
				{	failCorrectReg_ids.push_back(k);
				}
				break;
			}

			//先从左往右来做
			int highFreq_end = highFreqContRegs[i].end + KmerSize - 1;
			int lowFreq_end = highFreqContRegs[i+1].start - 1 + KmerSize - 1;
			int this_trim_right_end_len = 0;
			int this_right_last_change_pos = -1;

			int num_corrected = correct_multi_bases_rightward(read, freq, highFreq_end + 1, lowFreq_end, this_trim_right_end_len, 0, Max_change_in_one_read - accum_change, this_right_last_change_pos);
			
			if (this_trim_right_end_len == 0 && num_corrected > 0)
			{	correct_multi_base_score += num_corrected;
				accum_change += num_corrected;
			}
			
			//如果上一步没有纠对，再从右往左来纠错
			if (this_trim_right_end_len > 0 || num_corrected == 0)
			{	int highFreq_start = highFreqContRegs[i+1].start;
				int lowFreq_start = highFreqContRegs[i].end + 1;
				int this_trim_left_end_len = 0;
				int this_left_last_change_pos = -1;

				int num_corrected = correct_multi_bases_leftward(read, freq, highFreq_start - 1, lowFreq_start, this_trim_left_end_len, 0, Max_change_in_one_read - accum_change, this_left_last_change_pos);

				if (this_trim_left_end_len == 0 && num_corrected > 0)
				{	correct_multi_base_score += num_corrected;
					accum_change += num_corrected;
				}else
				{	failCorrectReg_ids.push_back(i);
				}
			}
		}
	}
	
	//成功纠错后相应的低频区域变为高频区域，把相邻的高频区域连起来，然后再找到最大的一个高频区
	//此处的核心目的是给出在前面计算基础上得到一个最大高频区，以便最后进行头和尾的纠错和trim操作	
	FreqContReg maxHighFreqReg;
	get_max_highFreq_region(highFreqContRegs, maxHighFreqReg, failCorrectReg_ids); //该步骤之后，highFreqContRegs仅保存一个max highfreq region
	

	//如果还有可改的量，则先修改头部
	//分支定界方法：修改开头上的错误base	
	int highFreq_start = maxHighFreqReg.start;
	if (highFreq_start > 1)
	{	if (accum_change < Max_change_in_one_read)
		{	int num_corrected = correct_multi_bases_leftward(read, freq, highFreq_start - 1, 1,	trim_left_end_len, 1, Max_change_in_one_read - accum_change, left_last_change_pos);
			if (num_corrected > 0)
			{	correct_multi_base_score += num_corrected;
				accum_change += num_corrected;
			}else
			{	trim_left_end_len = highFreq_start - 1;
				left_last_change_pos = 0;
			}
		}
		else
		{	trim_left_end_len = highFreq_start - 1;
			left_last_change_pos = 0;
		}
		
	}


	//如果还有可改的量，则最后修改尾部
	//分支定界方法：修改尾巴上的错误base
	int highFreq_end = maxHighFreqReg.end + KmerSize - 1;
	if (highFreq_end < read_len)
	{	if (accum_change < Max_change_in_one_read)
		{	int num_corrected = correct_multi_bases_rightward(read, freq, highFreq_end + 1, read_len, trim_right_end_len, 1, Max_change_in_one_read - accum_change, right_last_change_pos);
			if (num_corrected > 0)
			{	correct_multi_base_score += num_corrected;
				accum_change += num_corrected;
			}else
			{	trim_right_end_len = read_len - highFreq_end;
				right_last_change_pos = read_len + 1;
			}
		}
		else
		{	trim_right_end_len = read_len - highFreq_end;
			right_last_change_pos = read_len + 1;
		}
		
	}

	//在reads头部和尾部增加trim量, 提高最终结果的准确率	
	//如果trim_end_len不为零，则直接增加Further_trim_len；
	//否则，trim_end_len等于0，就看last_change_pos是否距离末端的距离在Further_trim_len之内
	if (trim_left_end_len > 0 || (left_last_change_pos > 0 && left_last_change_pos <= Further_trim_len))
	{	trim_left_end_len += Further_trim_len; 
		if (trim_left_end_len > read_len)
		{	trim_left_end_len = read_len;
		}
	}
	if (trim_right_end_len > 0 || (right_last_change_pos < read_len + 1 && right_last_change_pos >= read_len - Further_trim_len + 1))
	{	trim_right_end_len += Further_trim_len;
		if (trim_right_end_len > read_len)
		{	trim_right_end_len = read_len;
		}
	}
	
	//如果trimmed reads最终长度小于Min_trimmed_read_len，则删去
	if (read_len - trim_left_end_len - trim_right_end_len < Min_trimmed_read_len)
	{	is_deleted = 1;
	}
	
}

//get the max highFreq region, 先连接，然后再找最大
void get_max_highFreq_region(vector <FreqContReg> &highFreqContRegs, FreqContReg &maxHighFreqReg, vector<int> &failCorrectReg_ids)
{	
	vector <FreqContReg> CombinedhighFreqContRegs;
	int cur_start_pos = highFreqContRegs[0].start;
	failCorrectReg_ids.push_back(highFreqContRegs.size()-1);
	
	//combine the neighering high-freq-cont-regions, which intermediate have been corrected
	for (int i=0; i<failCorrectReg_ids.size(); i++)
	{	int stop_id = failCorrectReg_ids[i];
		int cur_end_pos = highFreqContRegs[stop_id].end;
		FreqContReg element;
		element.start = cur_start_pos;
		element.end = cur_end_pos;
		element.status = 1;
		CombinedhighFreqContRegs.push_back(element);

		if (stop_id != highFreqContRegs.size()-1)
		{	cur_start_pos = highFreqContRegs[stop_id+1].start;
		}
		
	}
	
	//find the max region from combined high-freq-cont-regions
	int max_len = 0;
	int max_id = 0;
	for (int i=0; i<CombinedhighFreqContRegs.size(); i++)
	{	int this_len = CombinedhighFreqContRegs[i].end - CombinedhighFreqContRegs[i].start + 1 ;
		if (this_len > max_len)
		{	max_len = this_len;
			max_id = i;
		}
	}
	
	maxHighFreqReg.start = CombinedhighFreqContRegs[max_id].start;
	maxHighFreqReg.end = CombinedhighFreqContRegs[max_id].end;
	maxHighFreqReg.status = 1;
}
	

//检查方向为：高频区kmer尾(左)->低频区kmer尾(右)
//只有当len_need_trim为0，并且num_corrected大于0的时候，才说明进行了成功彻底的纠错
//当因为可改量不够或者错误碱基太多无法彻底完成全部纠错任务时，对于无法纠错区域给出需要trim掉的长度信息
int correct_multi_bases_rightward(string &read, uint8_t *freq, int check_start, int check_end, int &len_need_trim, int is_modify_trimmed_reads, int max_allowed_change, int &last_change_pos)
{	
	int num_corrected = 0;
	
	string start_point_str = read.substr(check_start - KmerSize, KmerSize - 1); 
	uint64_t start_point_bit = seq2bit(start_point_str);
	
	if (max_allowed_change > 2)
	{	max_allowed_change = 2 ;  //每一次的可改量不能超过两个，也就是说，每一个连续低频区里面，最多可以发生两个测序错误
	}

	vector <TreeNode> node_vec;  //store all the high-frequency nodes in the treee
	vector <uint32_t> node_cur;  //the the postions in node_vec for the nodes in current parsing cycle
	TreeNode node;               //the root node of the tree
	node.change = 0;
	node_vec.push_back(node);
	node_cur.push_back(0);
	uint32_t node_vec_pos = 0;

	//build the tree data structure
	int cycle_cur = check_start;
	while (cycle_cur <= check_end)
	{	
		vector <uint32_t> node_tmp;
		
		for (uint32_t i = 0; i < node_cur.size(); i ++)
		{	
			uint32_t parent_pos = node_cur[i];
			for (uint32_t j =0; j < 4; j++)
			{	
				uint64_t kbit = 0;
				get_kmer_rightward(kbit, j, node_vec, parent_pos, start_point_bit);
				int kbit_freq = get_freq(freq, kbit);
				int is_equal_base = (bases[j] == read[cycle_cur - 1]) ? 1 : 0;
				int this_change = is_equal_base ? node_vec[parent_pos].change : node_vec[parent_pos].change + 1;
				if (kbit_freq == 1 && this_change <= max_allowed_change)
				{	TreeNode node;
					node.pointer = parent_pos;
					node.base = j;
					node.change = this_change;
					node.same = is_equal_base;
					node_vec.push_back(node);
					node_vec_pos ++;
					node_tmp.push_back(node_vec_pos);
					
				}
			}
		}
		//
		if (node_tmp.size() >= 1 && node_vec_pos < Max_node_in_BB_tree)
		{	node_cur = node_tmp;
		}else{
			if (node_vec_pos >= Max_node_in_BB_tree) 
			{ cerr << "node_vec_pos exceed Max_node_in_BB_tree" << endl; 
			}
			break;
		}
		
		//cerr << "cycle" << cycle_cur << endl;
		cycle_cur ++;
	}
	
	//cout << "##  " << node_vec_pos << endl;

	//get the mini change, node_cur.size()肯定是>=1的，但是里面可能只有一个根节点
	int cur_first = node_cur[0];
	int min_change = node_vec[cur_first].change;
	int min_change_pos = cur_first;
	int min_change_path = 0;
	for (int i = 0; i < node_cur.size(); i++)
	{	int cur_pos = node_cur[i];
		if (node_vec[cur_pos].change < min_change)
		{	min_change = node_vec[cur_pos].change;
			min_change_pos = cur_pos;
			min_change_path = 1;
		}
		else if(node_vec[cur_pos].change == min_change)
		{	min_change_path ++;
		}
	}

	
	len_need_trim = check_end - cycle_cur + 1; //完全结束时为0

	// trace back, correct the error base on reads directly
	if ( min_change_path == 1 &&  (len_need_trim ==0 || (len_need_trim >0 && is_modify_trimmed_reads)) )
	{	num_corrected = min_change;
		int read_len = read.size();
		uint32_t pos = min_change_pos;
		int read_pos = cycle_cur - 1;
		while (pos > 0)
		{	
			if (! node_vec[pos].same)
			{	read[read_pos - 1] = bases[node_vec[pos].base];
				if (last_change_pos == read_len + 1)
				{	last_change_pos = read_pos;
				}
			}
			pos = node_vec[pos].pointer;
			read_pos --;
		}
	}
	

	return  num_corrected;  //当返回值为0时，表明有多条路径
}

//get the kmer bit by trace back to the right forward
//加法操作的速度与位与位或操作速度相当
void get_kmer_rightward(uint64_t &kbit, uint8_t cur_base, vector <TreeNode> &node_vec, uint32_t pos, uint64_t start_point_bit)
{
	int i = 0;
	kbit =  bitLeft[cur_base];
	i ++;
	while (pos > 0 && i < KmerSize)
	{	
		kbit = (kbit >> 2) | bitLeft[node_vec[pos].base];
		pos = node_vec[pos].pointer;
		i ++;
	}
	
	while (i < KmerSize)
	{	kbit = (kbit >> 2) | bitLeft[start_point_bit & 0x3];
		start_point_bit >>= 2;
		i ++;
	}
	
	kbit >>= (64 - KmerSize * 2);
}


//检查方向为：高频区kmer头(右)->低频区kmer头(左)
//只有当len_need_trim为0，并且num_corrected大于0的时候，才说明进行了成功彻底的纠错
//当因为可改量不够或者错误碱基太多无法彻底完成全部纠错任务时，对于无法纠错区域给出需要trim掉的长度信息
int correct_multi_bases_leftward(string &read, uint8_t *freq, int check_start, int check_end, int &len_need_trim, int is_modify_trimmed_reads, int max_allowed_change, int &last_change_pos)
{	
	int num_corrected = 0;
	
	string start_point_str = read.substr(check_start, KmerSize - 1); 
	uint64_t start_point_bit = seq2bit(start_point_str);

	if (max_allowed_change > 2)
	{	max_allowed_change = 2 ;  //每一次的可改量不能超过两个，也就是说，每一个连续低频区里面，最多可以发生两个测序错误
	}

	vector <TreeNode> node_vec;   //store all the high-frequency nodes in the treee
	vector <uint32_t> node_cur;   //the the postions in node_vec for the nodes in current parsing cycle
	TreeNode node;                //the root node of the tree
	node.change = 0;
	node_vec.push_back(node);
	node_cur.push_back(0);
	uint32_t node_vec_pos = 0;

	//build the tree data structure
	int cycle_cur = check_start;
	while (cycle_cur >= check_end)
	{	
		vector <uint32_t> node_tmp;
		
		for (uint32_t i = 0; i < node_cur.size(); i ++)
		{	
			uint32_t parent_pos = node_cur[i];
			for (uint32_t j =0; j < 4; j++)
			{	
				uint64_t kbit = 0;
				get_kmer_leftward(kbit, j, node_vec, parent_pos, start_point_bit);
				int kbit_freq = get_freq(freq, kbit);
				int is_equal_base = (bases[j] == read[cycle_cur - 1]) ? 1 : 0;
				int this_change = is_equal_base ? node_vec[parent_pos].change : node_vec[parent_pos].change + 1;
				if (kbit_freq == 1 && this_change <= max_allowed_change)
				{	TreeNode node;
					node.pointer = parent_pos;
					node.base = j;
					node.change = this_change;
					node.same = is_equal_base;
					node_vec.push_back(node);
					node_vec_pos ++;
					node_tmp.push_back(node_vec_pos);
				}
			}
		}
		
		if (node_tmp.size() >= 1 && node_vec_pos < Max_node_in_BB_tree)
		{	node_cur = node_tmp;
		}else{
			if (node_vec_pos >= Max_node_in_BB_tree) 
			{ cerr << "node_vec_pos exceed Max_node_in_BB_tree" << endl; 
			}
			break;
		}
		
		//cerr << "cycle" << cycle_cur << endl;
		cycle_cur --;
	}
	
	//cout << "##  " << node_vec_pos << endl;

	//get the mini change, node_cur.size()肯定是>=1的，但是里面可能只有一个根节点
	int cur_first = node_cur[0];
	int min_change = node_vec[cur_first].change;
	int min_change_pos = cur_first;
	int min_change_path = 0;
	for (int i = 0; i < node_cur.size(); i++)
	{	int cur_pos = node_cur[i];
		if (node_vec[cur_pos].change < min_change)
		{	min_change = node_vec[cur_pos].change;
			min_change_pos = cur_pos;
			min_change_path = 1;
		}
		else if(node_vec[cur_pos].change == min_change)
		{	min_change_path ++;
		}
	}

	len_need_trim = cycle_cur - check_end + 1;  //完全结束时为0

	// trace back, get the min_change path and sequence
	if ( min_change_path == 1 &&  (len_need_trim ==0 || (len_need_trim >0 && is_modify_trimmed_reads)) )
	{	
		num_corrected = min_change;
		uint32_t pos = min_change_pos;
		int read_pos = cycle_cur + 1;
		while (pos > 0)
		{	if (! node_vec[pos].same)
			{	read[read_pos - 1] = bases[node_vec[pos].base];
				if (last_change_pos == 0)
				{	last_change_pos = read_pos;
				}
			}
			pos = node_vec[pos].pointer;
			read_pos ++;
		}
	}

	return  num_corrected;  //当返回值为0时，表明有多条路径
}

//get the kmer bit by trace back to the left forward
//加法操作的速度与位与位或操作速度相当
void get_kmer_leftward (uint64_t &kbit, uint8_t cur_base, vector <TreeNode> &node_vec, uint32_t pos, uint64_t start_point_bit)
{
	int i = 0;
	kbit =  cur_base;
	i ++;
	while (pos > 0 && i < KmerSize)
	{	
		kbit = (kbit << 2) | node_vec[pos].base;
		pos = node_vec[pos].pointer;
		i ++;
	}
	
	if (i < KmerSize)
	{	kbit = ( kbit << ((KmerSize-i)*2) ) | ( start_point_bit >> ((i-1)*2) );
	}
	
}


//parse one fasta format file
void parse_one_reads_fa_file(string &raw_reads_file)
{	
	uint64_t num_raw_reads = 0;
	uint64_t num_raw_bases = 0;
	uint64_t num_res_reads = 0;
	uint64_t num_res_bases = 0;
	uint64_t num_trimmed_reads = 0;
	uint64_t num_trimmed_bases = 0;
	uint64_t num_deleted_reads = 0;
	uint64_t total_one_base_correct_score = 0;
	uint64_t total_multi_base_correct_score = 0;
	uint64_t total_all_base_correct_score = 0;
	
	igzstream infile ( raw_reads_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file " << raw_reads_file << endl;
	}
	
	string reads_file_cor = raw_reads_file + ".cor";
	ogzstream corFile (reads_file_cor.c_str());
	if ( ! corFile )
	{	cerr << "fail to open output file " << reads_file_cor << endl;
	}
	
	//parse each read sequence, fa format, 容错读法
	string read, read_head;
	while ( getline( infile, read_head, '\n' ) )
	{	
		if (read_head[0] == '>') 
		{	getline( infile, read, '\n' );
		}
		else
		{	continue;	
		}

		num_raw_reads ++;
		num_raw_bases += read.size();

		int one_base_correct_score = 0;
		int multi_base_correct_score = 0;
		int all_base_correct_score = 0;
		int is_deleted = 0;
		int trim_left_end_len = 0;
		int trim_right_end_len = 0;
		int final_read_len = 0;
		
		correct_one_read(read_head, read, KmerFreq, one_base_correct_score, multi_base_correct_score, is_deleted, trim_left_end_len, trim_right_end_len);

		all_base_correct_score = one_base_correct_score + multi_base_correct_score;
		final_read_len = read.size() - trim_left_end_len - trim_right_end_len;

		if ( ! is_deleted )
		{	
			total_one_base_correct_score += one_base_correct_score;
			total_multi_base_correct_score += multi_base_correct_score;
			if (trim_left_end_len > 0 || trim_right_end_len > 0)
			{	
				read = read.substr(trim_left_end_len, final_read_len);
				num_trimmed_reads ++;
				num_trimmed_bases += trim_left_end_len + trim_right_end_len;
			}
			corFile << read_head << " score: " << all_base_correct_score << "  left_trim: " << trim_left_end_len << "\n" << read << "\n";
			num_res_reads ++;
			num_res_bases += final_read_len;
		}
		else
		{	
			corFile << read_head << " score: " << all_base_correct_score << "  left_trim: " << trim_left_end_len << "\n\n";
			num_deleted_reads ++;
		}
	}
	
	cerr << "Finished to parse reads file: " << raw_reads_file << endl;

	total_all_base_correct_score = total_one_base_correct_score + total_multi_base_correct_score;
	
	//output statistic result
	string reads_file_cor_stat = raw_reads_file + ".cor.stat";
	ofstream staFile (reads_file_cor_stat.c_str());
	if ( ! staFile )
	{	cerr << "fail to open output file " << reads_file_cor_stat << endl;
	}
	
	double data_filter_ratio = (double)(num_raw_bases - num_res_bases) / (double)num_raw_bases;
	double corrected_error_ratio = (double)total_all_base_correct_score / (double)num_res_bases;

	staFile << "num_raw_reads " << num_raw_reads << endl;
	staFile << "num_raw_bases " << num_raw_bases << endl;
	staFile << "num_result_reads " << num_res_reads << endl;
	staFile << "num_result_bases " << num_res_bases << endl;
	
	staFile << "\nnum_trimmed_reads " << num_trimmed_reads << endl;
	staFile << "num_trimmed_bases " << num_trimmed_bases << endl;
	staFile << "num_deleted_reads " << num_deleted_reads << endl;
	
	staFile << "\nnum_corrected_bases_by_Fast_method " << total_one_base_correct_score << endl;
	staFile << "num_corrected_bases_by_BBtree_method " << total_multi_base_correct_score << endl;
	staFile << "num_corrected_bases_by_two_methods " << total_all_base_correct_score << endl;
	
	staFile << "\nlow_quality_bases_filter_ratio " << data_filter_ratio << endl;
	staFile << "estimated_raw_base_error_ratio " << corrected_error_ratio << endl;

}


//parse one fastaq format file
void parse_one_reads_fq_file(string &raw_reads_file)
{	
	uint64_t num_raw_reads = 0;
	uint64_t num_raw_bases = 0;
	uint64_t num_res_reads = 0;
	uint64_t num_res_bases = 0;
	uint64_t num_trimmed_reads = 0;
	uint64_t num_trimmed_bases = 0;
	uint64_t num_deleted_reads = 0;
	uint64_t total_one_base_correct_score = 0;
	uint64_t total_multi_base_correct_score = 0;
	uint64_t total_all_base_correct_score = 0;
	
	igzstream infile ( raw_reads_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file " << raw_reads_file << endl;
	}
	
	string reads_file_cor = raw_reads_file + ".cor";
	ogzstream corFile (reads_file_cor.c_str());
	if ( ! corFile )
	{	cerr << "fail to open output file " << reads_file_cor << endl;
	}
	
	//parse each read sequence, fq format, 容错读法
	string read, read_head, empty_line;
	while ( getline( infile, read_head, '\n' ) )
	{	
		if (read_head[0] == '@') 
		{	read_head[0] = '>';
			getline( infile, read, '\n' );
			getline( infile, empty_line, '\n' );
			getline( infile, empty_line, '\n' );
		}

		num_raw_reads ++;
		num_raw_bases += read.size();

		int one_base_correct_score = 0;
		int multi_base_correct_score = 0;
		int all_base_correct_score = 0;
		int is_deleted = 0;
		int trim_left_end_len = 0;
		int trim_right_end_len = 0;
		int final_read_len = 0;
		
		
		correct_one_read(read_head, read, KmerFreq, one_base_correct_score, multi_base_correct_score, is_deleted, trim_left_end_len, trim_right_end_len);

		all_base_correct_score = one_base_correct_score + multi_base_correct_score;
		final_read_len = read.size() - trim_left_end_len - trim_right_end_len;

		if ( ! is_deleted )
		{	
			total_one_base_correct_score += one_base_correct_score;
			total_multi_base_correct_score += multi_base_correct_score;
			if (trim_left_end_len > 0 || trim_right_end_len > 0)
			{	
				read = read.substr(trim_left_end_len, final_read_len);
				num_trimmed_reads ++;
				num_trimmed_bases += trim_left_end_len + trim_right_end_len;
			}
			corFile << read_head << " score: " << all_base_correct_score << "  left_trim: " << trim_left_end_len << "\n" << read << "\n";
			num_res_reads ++;
			num_res_bases += final_read_len;
		}
		else
		{	
			corFile << read_head << " score: " << all_base_correct_score << "  left_trim: " << trim_left_end_len << "\n\n";
			num_deleted_reads ++;
		}
	}
	
	cerr << "Finished to parse reads file: " << raw_reads_file << endl;

	total_all_base_correct_score = total_one_base_correct_score + total_multi_base_correct_score;
	
	//output statistic result
	string reads_file_cor_stat = raw_reads_file + ".cor.stat";
	ofstream staFile (reads_file_cor_stat.c_str());
	if ( ! staFile )
	{	cerr << "fail to open output file " << reads_file_cor_stat << endl;
	}
	
	double data_filter_ratio = (double)(num_raw_bases - num_res_bases) / (double)num_raw_bases;
	double corrected_error_ratio = (double)total_all_base_correct_score / (double)num_res_bases;

	staFile << "num_raw_reads " << num_raw_reads << endl;
	staFile << "num_raw_bases " << num_raw_bases << endl;
	staFile << "num_result_reads " << num_res_reads << endl;
	staFile << "num_result_bases " << num_res_bases << endl;
	
	staFile << "\nnum_trimmed_reads " << num_trimmed_reads << endl;
	staFile << "num_trimmed_bases " << num_trimmed_bases << endl;
	staFile << "num_deleted_reads " << num_deleted_reads << endl;
	
	staFile << "\nnum_corrected_bases_by_Fast_method " << total_one_base_correct_score << endl;
	staFile << "num_corrected_bases_by_BBtree_method " << total_multi_base_correct_score << endl;
	staFile << "num_corrected_bases_by_two_methods " << total_all_base_correct_score << endl;
	
	staFile << "\nlow_quality_bases_filter_ratio " << data_filter_ratio << endl;
	staFile << "estimated_raw_base_error_ratio " << corrected_error_ratio << endl;

}

//merge the read1 and read2 files into pair and single file
void merge_two_corr_files(string &read1_file, string &read2_file)
{
	igzstream READ1 ( read1_file.c_str() );
	if ( ! READ1 )
	{	cerr << "fail to open input file " << read1_file << endl;
	}

	igzstream READ2 ( read2_file.c_str() );
	if ( ! READ2 )
	{	cerr << "fail to open input file " << read2_file << endl;
	}
	
	string pair_file = read1_file + ".pair.fa.gz";
	ogzstream PAIR (pair_file.c_str());
	if ( ! PAIR )
	{	cerr << "fail to open output file " << pair_file << endl;
	}

	string single_file = read1_file + ".single.fa.gz";
	ogzstream SINGLE (single_file.c_str());
	if ( ! SINGLE )
	{	cerr << "fail to open output file " << single_file << endl;
	}

	uint64_t pair_reads_num = 0;
	uint64_t pair_bases_num = 0;
	uint64_t single_reads_num = 0;
	uint64_t single_bases_num = 0;

	//suppose the read number and line number are equal in read1_file and read2_file
	string read1, read_head1, read2, read_head2;
	while ( getline( READ1, read_head1, '\n' ) )
	{	
		getline( READ1, read1, '\n' );
		
		//get one reads from the read2_file
		getline( READ2, read_head2, '\n' );
		getline( READ2, read2, '\n' );
		
		if (read1.size() && read2.size())
		{
			PAIR << read_head1 << "\n" << read1 << "\n" << read_head2 << "\n" << read2 << "\n"; 
			pair_reads_num += 2;
			pair_bases_num += read1.size() + read2.size();
		}
		else
		{	
			if ( read1.size() )
			{	SINGLE << read_head1 << "\n" << read1 << "\n";
				single_reads_num ++;
				single_bases_num += read1.size();
			}
			if ( read2.size() )
			{	SINGLE << read_head2 << "\n" << read2 << "\n";
				single_reads_num ++;
				single_bases_num += read2.size();
			}
		}
	}

	//output the stat file
	string pair_single_stat_file = read1_file + ".pair.single.stat";
	ofstream PSS (pair_single_stat_file.c_str());
	if ( ! PSS )
	{	cerr << "fail to open output file " << pair_single_stat_file << endl;
	}
	PSS << "pair reads:   " << pair_reads_num << endl;
	PSS << "pair bases:   " << pair_bases_num << endl;
	PSS << "single reads: " << single_reads_num << endl;
	PSS << "single bases: " << single_bases_num << endl;

}
