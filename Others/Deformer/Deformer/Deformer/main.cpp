#include <iostream>
#include <fstream>
#include <string>


using namespace std;

static std::string INPUTFILE{ "in.txt" };   //输入文件名称
static std::string OUTPUTFILE{ "out.txt" }; //输出文件名称 

int main(int argc ,char* argv)
{
	//////////////////////////////////////////////////////////////////////////
	//文件读取
	std::ifstream infile(INPUTFILE);
	
	if (!infile.is_open())
	{
		cout << "未成功打开输入文件" << endl;
		return 0;
	}

	string line;
	while (getline(infile, line))
	{
		cout << line << endl;
		
		///操作字符串
		


	}
	infile.close();
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//实现算法在这里调用

	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//计算完成输出结果
	ofstream outfile(OUTPUTFILE);
	if (!outfile.is_open())
	{
		cout << "未成功打开输出文件" << endl;
		return 0;
	}

	outfile << 1 << "  " << 2.03 << endl;
	outfile << "end";

	outfile.close();

	//////

	return 1;
}