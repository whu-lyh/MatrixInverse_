#include<iostream>
#include<fstream>
#include<sstream>//istringstream 必须包含这个头文件
#include<cstdio>
#include<iomanip>
#include<vector>
#include<string>
using namespace std;
#define isZero     0.00000001
#define maxn       50
struct coef_matrix {
	int row;
	vector<vector<double>> matrix;
	vector<vector<double>> matrix_p;
};
 coef_matrix read()//读入系数矩阵
 {
	cout << "please input the coefficient matrix: "<<endl;
	ifstream infile("F:\\VS2015\\gauss_solve\\DATA.txt");
	if (!infile){
		cout << "Cannot open point cloud file." << endl;
	}
	int row = 0,i=0,j=0;
	string line;
	getline(infile, line);
	row = stoi(line);//系数的阶数
	vector<vector<double>> a(row);//总的系数矩阵，可以捕不设置二维向量的维度
	for (i = 0; i != row; ++i) {
		getline(infile, line); //getline命令之前用过一次之后就会在后面接着运行
		istringstream record(line);
		string a1;
		while (record >> a1){
			a[i].push_back(stod(a1));
		}
	}
	cout << "^^^^^^^^^Input coefficient matrix:" << endl;
	for (unsigned i = 0; i != row; ++i) {
		for (unsigned j = 0; j != row; ++j)
			cout <<right<<fixed<<setw(9)<< a[i][j] << " ";
		cout << endl;
	}
	coef_matrix matrix_a;//结构体实例化
	matrix_a.matrix = a;
	matrix_a.row = row;
	infile.close();//关闭文件输入流
	return matrix_a;//返回系数阵的结构体
}
 //采用选主元的三角分解方法
 coef_matrix doolittele(coef_matrix m_a)//传入系数矩阵
 {
	 int row = m_a.row;
	 vector<double> aa(row);
	 vector<vector<double>> a = m_a.matrix;
	 vector<vector<double>> temp(row, vector<double>(row, 0));//用于寻找最大主元
	 vector<vector<double>> p(row, vector<double>(row, 0));//交换矩阵
	 unsigned i = 0, r = 0, k = 0;
	 for (i = 0; i != row; ++i)
		 for (r = 0; r != row; ++r) {
			 temp[i][r] = 0;
			 if (r == i) {
				 p[i][r] = 1;
			 }
		 }
	 //考虑到用向量可能不是很方便，基于练习的话，果断改成数组的形式，反正也知道了系数矩阵的维度。
	 //double a[row][row];不能用非常量表达式...
	 double sum = 0, temp_p = 0;
	 vector<int> max(row);
	 for (r = 1; r != row; ++r)//r代表列,这里需要注意
		 for (i = r; i != row; ++i) {//i代表行
			 for (k = 0; k != r - 1; ++k) {//r-1还是什么，还需要进一步搞清楚
				 sum += (a[i][k] /= a[k][r]);
			 }
			 a[i][r] -= sum;//si
			 sum = 0;
		 }
	 for (r = 0; r != row; ++r) {//r代表列,这里需要注意
		 max[r] = r;//默认对角线是最大的
		 for (i = r; i != row; ++i) {//i代表行
			 if (abs(a[i][r]) > abs(a[r][r])) {//寻找最大的主元，并交换这两行的元素
				 max[r] = i;//发现更大的就跟新一次
				 for (k = 0; k != row; ++k) {//和最大的那个元素对换
					 temp[i][k] = a[r][k];
					 a[r][k] = a[i][k];
					 a[i][k] = temp[i][k];
				 }
			 }
		 }
	 }//??????????
	 for (i = 0; i != row; ++i)
		 for (r = i; r != row; ++r) {
				 if (max[r] != i) {
					temp[i][r]=p[max[r]][r];//交换矩阵
					p[max[r]][r]=p[i][r];
					p[i][r] = temp[i][r];
				}
		 }
//还应该确认一下det是不是=0？
	for (r = 0; r != row-1; ++r)//用L,U冲散A矩阵。
		for (i = r+1; i != row; ++i) {
			a[i][r] /= a[r][r];//下三角矩阵
			if (r!=0) {
				for (k = 0; k != r - 1; ++k) {
					sum += (a[r][k] /= a[k][i]);
				}
			}
			a[r][i] -= sum;//上三角矩阵
			sum = 0;
		}
	coef_matrix doolittle_coef_m;
	doolittle_coef_m.row = row;
	doolittle_coef_m.matrix = a;
	doolittle_coef_m.matrix_p = p;
	return doolittle_coef_m;
}
//求逆
coef_matrix inverse(const coef_matrix a) {//感觉需要L,U,P分开存储分开操作啊
	int row = a.row;
	coef_matrix inverse_m, inverse_mU, inverse_mV, inverse_mL;
	inverse_m.row = row;
	//inverse_mU.row = row;
	inverse_mU.matrix = a.matrix;
	//inverse_mV.row = row;
	//inverse_mL.row = row;
	vector<double> vv(row,0);
	for (unsigned i = 0; i != row; ++i) {
		inverse_m.matrix.push_back(vv);
		inverse_mV.matrix.push_back(vv);
		inverse_mL.matrix.push_back(vv);
		inverse_mU.matrix_p.push_back(vv);
	}
	for (unsigned r = 0; r != row; ++r) {
		inverse_mV.matrix[r][r] = 1 / a.matrix[r][r];
		inverse_mL.matrix[r][r] = 1;
		for (unsigned i = 0; i != row; ++i) {//对角线元素,r还是表示列
			if (i > r) {
				inverse_mU.matrix[i][r] = 0;
			}
		}
	}
	double sum = 0;
	for (unsigned  i = row - 2,r = row-1; r != row, i != 0; ++r, --i)//U上三角矩阵的逆阵
		//for (unsigned i = r-1 ; i !=row; ++i) 
	{
			for (unsigned k = i+1; k != r; ++k)//这个k也表示列
				sum += inverse_mU.matrix[i][k] * inverse_mU.matrix[k][r];
			inverse_mV.matrix[i][r] = -sum/inverse_mV.matrix[i][i];
		}
	for (unsigned r = 0; r != row-1; ++r) //L下三角矩阵的逆阵
		for (unsigned i = r+1; i != row; ++i) {
			inverse_mL.matrix[i][r] = -1 / a.matrix[i][r];
		}
	inverse_m.matrix_p =a.matrix_p ;
	for (unsigned i = 0; i != row; ++i) 
		for (unsigned r = 0; r != row; ++r) 
			for (unsigned k = 0; k != row; ++k) 
				inverse_mU.matrix_p[i][r] += inverse_mV.matrix[i][k]*inverse_mL.matrix[k][r];
	for (unsigned i = 0; i != row; ++i)
		for (unsigned r= 0; r != row; ++r)
			for (unsigned k = 0; k != row; ++k)
				inverse_m.matrix[i][r] += inverse_mU.matrix_p[i][k]*inverse_m.matrix_p[k][r];
	return inverse_m;
}

coef_matrix feedback(const coef_matrix a, const coef_matrix b) {
	coef_matrix feed_m;
	int row = a.row;
	vector<double> vv(row, 0);
	for (unsigned i = 0; i != row; ++i) 
		feed_m.matrix.push_back(vv);
	for (unsigned r = 0; r != row; ++r)
		for(unsigned i=0;i!=row;++i)
			for (unsigned k = 0; k != row; ++k) {
				feed_m.matrix[i][r] +=a.matrix[i][k] * b.matrix[k][r];
			}
	feed_m.row = row;
	return feed_m;
}

coef_matrix print(const coef_matrix a,const int b)//打印矩阵
{
	cout << "\n";
	int row = a.row;
	switch (b){
		case 1: {
			for (unsigned i = 0; i != row; ++i) {
				for (unsigned j = 0; j != row; ++j)
					cout << right << fixed << setw(9) << a.matrix[i][j] << " ";
				cout << endl;
			}
			break;
		}
		case 2: {
			for (unsigned i = 0; i != row; ++i) {
				for (unsigned j = 0; j != row; ++j)
					cout << right << fixed << setw(9) << a.matrix_p[i][j] << " ";
				cout << endl;
			}
			break;//break不可少
		}
		default:
			cout << "row: "<<row << endl;
	}
	return a;
}

void main()
{
	coef_matrix matrix_a;
	coef_matrix matrix_LU;
	coef_matrix matrix_inv;
	coef_matrix matrix_f;
	matrix_a=read();
	matrix_LU=doolittele(matrix_a);
	cout << "LU_matrix:";
	print(matrix_LU,1);//调用时一定要加小括号
	cout << "sort matrix:";
	print(matrix_LU,2);
	matrix_inv = inverse(matrix_LU);
	cout << "Inv_matrix:";
	print(matrix_inv,1);
	matrix_f=feedback(matrix_a, matrix_inv);
	cout << "Comfirm: ";
	print(matrix_f,1);
	system("pause");
}
