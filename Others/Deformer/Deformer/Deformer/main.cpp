#include <iostream>
#include <fstream>
#include <string>


using namespace std;

static std::string INPUTFILE{ "in.txt" };   //�����ļ�����
static std::string OUTPUTFILE{ "out.txt" }; //����ļ����� 

int main(int argc ,char* argv)
{
	//////////////////////////////////////////////////////////////////////////
	//�ļ���ȡ
	std::ifstream infile(INPUTFILE);
	
	if (!infile.is_open())
	{
		cout << "δ�ɹ��������ļ�" << endl;
		return 0;
	}

	string line;
	while (getline(infile, line))
	{
		cout << line << endl;
		
		///�����ַ���
		


	}
	infile.close();
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//ʵ���㷨���������

	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//�������������
	ofstream outfile(OUTPUTFILE);
	if (!outfile.is_open())
	{
		cout << "δ�ɹ�������ļ�" << endl;
		return 0;
	}

	outfile << 1 << "  " << 2.03 << endl;
	outfile << "end";

	outfile.close();

	//////

	return 1;
}