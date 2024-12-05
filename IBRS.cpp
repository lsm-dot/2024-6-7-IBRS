#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
//********* choose just one of these pairs **********
//#define MR_PAIRING_CP      // AES-80 security   
//#define AES_SECURITY 80

#define MR_PAIRING_MNT	// AES-80 security
#define AES_SECURITY 80

//#define MR_PAIRING_BN    // AES-128 or AES-192 security
//#define AES_SECURITY 128
//#define AES_SECURITY 192

//#define MR_PAIRING_KSS    // AES-192 security
//#define AES_SECURITY 192

//#define MR_PAIRING_BLS    // AES-256 security
//#define AES_SECURITY 256
//*********************************************

#include "pairing_3.h"
#include "ecn.h"
using namespace std;


extern "C" {
#include "miracl.h"
#include "mirdef.h"
}
PFC pfc(AES_SECURITY);  // initialise pairing-friendly curve
struct KeyParams {
	Big q;
	Big s_trac;
	Big s;
	G2 PK_trac;
	G1 PK1;
	G2 PK2;
	G1 P;
	G2 Q;
};

// KeyGen1 function
void KeyGen1(KeyParams& params, char*& VID, G1& PID, G1& PSK) {
	pfc.hash_and_map(PID, VID);
	PSK = pfc.mult(PID, params.s);
	/*return PID, PSK;*/
}

// KeyGen2 function
void KeyGen2(KeyParams& params, char*& VID, G2& RID, G2& RSK) {

	pfc.hash_and_map(RID, VID);

	RSK = pfc.mult(RID, params.s);
	//return make_pair(RID, RSK);
}

void sign(KeyParams& params, vector<G1>& Ls, G1& PSK, char* VIDk, int& k, char* m, vector<G1>& U_array, G1& V, GT& tag) {
	int len_L = Ls.size();
	U_array.resize(len_L);
	G1 sum_Ls, RIDk;
	sum_Ls.g.clear();
	int flag = 0;
	//cout << sum_Ls.g << endl;

	pfc.hash_and_map(RIDk, VIDk);
	tag = pfc.pairing(params.PK_trac, RIDk);
	//Big temp_m = pfc.hash_to_group((char*)m);
	Big Ls_result = 0;
	pfc.start_hash();
	for (int i = 0; i < len_L; i++)
	{
		pfc.add_to_hash(Ls[i]);
	}
	Ls_result = pfc.finish_hash_to_group();
	Big hi;
	G1 Ui, mul;
	for (int i = 0; i < len_L; i++) {
		if (i == k)
		{
			continue;
		}
		
		pfc.random(Ui);
		U_array[i] = Ui;

		pfc.start_hash();
		pfc.add_to_hash(m);
		pfc.add_to_hash(tag);
		pfc.add_to_hash(Ls_result);
		pfc.add_to_hash(Ui);
		hi = pfc.finish_hash_to_group();
		
		mul = pfc.mult(Ls[i], hi);
		if (!flag)
		{
			sum_Ls = U_array[i] + mul;
			flag = 1;
		}
		else {
			sum_Ls = sum_Ls + U_array[i] + mul;
		}

	}

	Big r;
	pfc.random(r);
	G1 Uk;
	sum_Ls = -sum_Ls;
	Uk = pfc.mult(Ls[k], r) + sum_Ls;
	U_array[k] = Uk;
	pfc.start_hash();
	pfc.add_to_hash(m);
	pfc.add_to_hash(tag);
	pfc.add_to_hash(Ls_result);
	pfc.add_to_hash(Uk);
	Big hk = pfc.finish_hash_to_group();
	//Big hk = H(m, tag, Ls, Uk);
	V = pfc.mult(PSK, hk + r);
}

bool single_verify(KeyParams& params, char* m, vector<G1>& Ls, vector<G1>& U_array, G1& V, GT& tag) {
	/*const vector<G1>& U_array = (sigma);
	const G1& V = get<1>(sigma);
	const GT& tag = get<2>(sigma);*/
	int flag = 0;
	int len_L = Ls.size();
	G1 sum_Ls;
	sum_Ls.g.clear();
	//sum_Ls = 0;
	//Big temp_m = pfc.hash_to_group((char*)m);
	Big Ls_result = 0;

	GT righteq = pfc.pairing(params.Q, V);
	pfc.start_hash();
	for (int i = 0; i < len_L; i++)
	{
		pfc.add_to_hash(Ls[i]);
	}
	Ls_result = pfc.finish_hash_to_group();
	for (int i = 0; i < len_L; i++) {
		pfc.start_hash();
		pfc.add_to_hash(m);
		pfc.add_to_hash(tag);
		pfc.add_to_hash(Ls_result);
		pfc.add_to_hash(U_array[i]);
		Big hi = pfc.finish_hash_to_group();
		//Big hi = H(m, tag, Ls, U_array[i]);
		G1 mul = pfc.mult(Ls[i], hi);;
		//mul *= hi;
		
		if (!flag)
		{
			sum_Ls = U_array[i] + mul;
			flag = 1;
		}
		else {
			sum_Ls = sum_Ls + U_array[i] + mul;
		}
		
	}

	GT lefteq = pfc.pairing(params.PK2, sum_Ls);
	

	if (lefteq != righteq) {
		cout << "single verify error!" << endl;
		return false;
	}
	return true;
}

// Function to calculate trace
char* trace(KeyParams& params, GT& tag, vector<G1>& Ls, vector<char*>& VID) {
	GT tag1 = pfc.power(tag, inverse(params.s_trac, params.q)); // Compute tag1 = tag^(-s_trac)

	for (size_t i = 0; i < Ls.size(); ++i) {
		G1 Hi;
		pfc.hash_and_map(Hi, VID[i]); // Hash VID[i] to G1

		GT pairing_result = pfc.pairing(params.Q, Hi); // Compute pairing result with Q

		if (pairing_result == tag1) {
			return VID[i]; // Return the matching VID
		}
	}
}

void test_trial(KeyParams& params) {
	vector<double> sign_time_list;
	vector<double> trace_time_list;
	vector<double> verify_time_list;
	vector<double> entire_time_list;

	cout << "power\tsign\tverify\ttrace\tentire(ms)" << endl;

	for (int power = 0; power < 6; power++) {
		double sign_time_sum = 0;
		double trace_time_sum = 0;
		double verify_time_sum = 0;
		double entire_time_sum = 0;
		int power_of_2 = power + 1;
		int PK_num = (int)pow(2, power_of_2);
		int time_trail = 10;

		vector<G1> fake_Pid(PK_num);
		vector<G1> fake_Ssk(PK_num);
		vector<char*> VID(PK_num);
		
		
		for (int j = 0; j < PK_num; j++) {
			VID[j] = (char*)("hello world!");
			
			KeyGen1(params, VID[j], fake_Pid[j], fake_Ssk[j]);
		}
		
		time_t seed;
		time(&seed);
		irand((long)seed);

		for (int ii = 0; ii < time_trail; ii++) {
			clock_t start_time = clock();
			int random_position = rand(params.q) % PK_num;
			//int random_position = 0;
			G1 pid;
			G1 psk;
			KeyGen1(params, VID[random_position], pid, psk);
			fake_Pid[random_position] = pid;

			vector<G1> U_array(PK_num);
			G1 V;
			GT tag;
			char *message = (char*)"I am a girl";
			clock_t sign_start_time = clock();
			sign(params, fake_Pid, psk, VID[random_position], random_position, message, U_array, V, tag);
			clock_t sign_end_time = clock();
			sign_time_sum += ((double)sign_end_time - (double)sign_start_time) / CLOCKS_PER_SEC;

			clock_t trace_start_time = clock();
			char* VID_trace = trace(params, tag, fake_Pid, VID);
			clock_t trace_end_time = clock();
			trace_time_sum += ((double)trace_end_time - (double)trace_start_time) / CLOCKS_PER_SEC;

			clock_t verify_start_time = clock();
			if (!single_verify(params, message, fake_Pid, U_array, V, tag)) {
				cout << "failed" << endl;
			}
			clock_t verify_end_time = clock();
			verify_time_sum += ((double)verify_end_time - (double)verify_start_time) / CLOCKS_PER_SEC;

			/*clock_t trace_start_time = clock();
			char* VIDk = trace(params, pfc, tag, fake_Pid, &VID[0]);
			clock_t trace_end_time = clock();
			trace_time_sum += ((double)trace_end_time - (double)trace_start_time) / CLOCKS_PER_SEC;*/
		
			clock_t end_time = clock();
			entire_time_sum += ((double)end_time - (double)start_time) / CLOCKS_PER_SEC;
		}

		cout << power + 1 << "\t" << (sign_time_sum / time_trail) * 1000 << "\t" << (verify_time_sum / time_trail) * 1000
			<< "\t" << (trace_time_sum / time_trail) * 1000 << "\t" << (entire_time_sum / time_trail) * 1000 << endl;

		sign_time_list.push_back((sign_time_sum / time_trail) * 1000);
		trace_time_list.push_back((trace_time_sum / time_trail) * 1000);
		verify_time_list.push_back((verify_time_sum / time_trail) * 1000);
		entire_time_list.push_back((entire_time_sum / time_trail) * 1000);
	}
	// Get current date
	std::time_t t = std::time(nullptr);
	std::tm now;
	localtime_s(&now, &t);

	// Format date as YYYY-MM-DD
	std::ostringstream date_stream;
	date_stream << (now.tm_year + 1900) << '-'
		<< (now.tm_mon + 1) << '-'
		<< now.tm_mday << '-' << now.tm_hour << '-' << now.tm_min << '-';
	std::string date_str = date_stream.str();

	// Create filename with date
	std::string filename = date_str + "Identity-Based Ring Signature Time Analysis.txt";

	// 指定文件路径
	std::string directory = "D:\\VS2019\\repo\\2024-6-7-IBRS\\time_test\\";
	//std::string filename = "output.txt";
	std::string filepath = directory + filename;

	// Writing to file
	std::ofstream text_file;
	text_file.open(filepath);
	text_file << "2^n\tSign\tVerify\tTrace\n";
	for (size_t i = 0; i < sign_time_list.size(); i++) {
		text_file << i + 1 << "\t" << sign_time_list[i] << "\t" << verify_time_list[i] << "\t" << trace_time_list[i] << endl;
	}
	text_file.close();
}

int main()
{
	

	Big q = pfc.order();
	Big s_trac, s;
	G1 P, PK1; 
	G2 Q, PK_trac, PK2;
	time_t seed;

	time(&seed);
	irand((long)seed);

	// setup
	pfc.random(P);
	pfc.random(Q);

	// 生成LEA secret key and public key
	pfc.random(s_trac);
	PK_trac = pfc.mult(Q, s_trac);

	// generate master secret key and public key
	pfc.random(s);

	PK1 = pfc.mult(P, s);
	PK2 = pfc.mult(Q, s);

	KeyParams params; // 定义公共参数
	//params.pfc = pfc;
	params.q = q;
	params.P = P;
	params.Q = Q;
	params.s_trac = s_trac;
	params.s = s;
	params.PK_trac = PK_trac;
	params.PK1 = PK1;
	params.PK2 = PK2;

	test_trial(params);

	return 0;
}