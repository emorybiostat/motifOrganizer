// clustering gene expression profiles 11/07/05
// programmer: steve qin
// 02/16/06 add command line arguments
// modified 04/11/06, add
// 1 - standardization
// 2 - co-occurance matrix 

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<list>
#include<string>
#include<math.h>
#include <time.h>
#include "ranlib.h"
using namespace std;

//#define INITIAL_CLUST 15//20
#define diffe 0//use difference =1, use original =0
//#define CHAINS 10
//#define ROUND 20
#define FACTOR_0 1
#define pie 3.1415926
//#define invert false//true//false 
#define MAX_LENGTH 150
//#define row 10
//#define column 5
//#define group 2
//#define priorbeta0 -0.08//7.7//5//0
//#define priora 2.08//17//9.5//1.111//3
//#define priorb 0.086//240//127.5//2.111//5

class clusterClass
{
private:
public:
	list <int> geneID;
	vector <double> ysum;
	vector <double> ysum2;
	//double *ysum;
	//double *ysum2;
/*	clusterClass(int j, int length)
	{
		geneID.push_back(j);
		ysum = new double[length];
		ysum2 = new double[length];
	}*/
};

/*class geneClass   
{
private:
public:
	string name;
	int trueClass;
	//list <clusterClass>::iterator member;
	//double postProb;
	//int shift;
};
*/
void readData(const char dataFileName[], vector <double> &data, int *nCol);

void standard(vector <double> &data, const int nCol, const int nRow);

double gammaln(double xx);

void operation(const int geneID, const int length, const int maxShift, const double priorbeta0, const double priora, double priorb);
// const int column,

void operationinvert(const int geneID, const int length, const int maxShift, const double priorbeta0, 
		const double priora, double priorb,const int invert);
// const int column,

double bayesratio(const int geneID, const int shift, const int length, const double *ysum, 
				  const double *ysum2, const int nsize, const double priorbeta0, const double priora, double priorb);// const int column,

double bayesratioinvert(const int geneID, const int shift, const int length, const double *ysum, 
				  const double *ysum2, const int nsize, const double priorbeta0, const double priora, double priorb);// const int column,

//void simudata();

int singleCompareInvert(const vector <double> ratioVec, bool &outsign,const int invert);

int singleCompare(const vector <double> ratioVec);

int compareInvert(const vector <double> ratioVec, const int maxShift, int &shift, bool &outsign,const int invert);

int compare(const vector <double> ratioVec, const int maxShift, int &shift);

double unran(int *na, int *nb, int *nc);

double clust(const int row, const int length, const int maxShift,list <clusterClass> &finalClust,
			 int *&finalShift, bool *&finalSign, const double priorbeta0, 
			 const double priora, double priorb, const int ROUND,const int invert, int **&occur);//column

double loglikelihood(const int length, const double priorbeta0, const double priora, double priorb);//column

double finalRatio(list <clusterClass>::iterator allClustsIter, list <clusterClass> finalClust,
	       int *finalShift, bool *finalSign, const int length, 
		const double priorbeta0, const double priora, double priorb,const int invert);

void output(const char outPutFileName[], const int row, const int column, const int length, const double loglike,
			list <clusterClass> finalClust, int *finalShift, bool *finalSign, double *postProba, 
			const double priorbeta0, const double priora, double priorb, const char memberFileName[],
                        const char lengthFileName[], const char resultFileName[],const int invert,
                        const double probcutoff,const int CHAINS, const int ROUND, const int misscount, int **occur);

void postProb(const int row, const int length, const int maxShift, list <clusterClass> finalClust,
			  int *finalShift, bool *finalSign, double *&postProba, const double priorbeta0, 
			const double priora, double priorb,const int invert);

void oldpostProb(const int row, const int length, const int maxShift, list <clusterClass> finalClust,
			  int *finalShift, double *&postProba, const double priorbeta0, const double priora, double priorb);

void randomInitialClust(const int row, const int length, const int maxShift, const int updatedClustNumber,const int invert);

void cooccur(const int row, int **&occur);

list <clusterClass> allClusts;
//vector <geneClass> allGenes;
list <clusterClass>::iterator * member;
bool *sign;
bool *missing;
double *average;
vector <string> geneName;
//308 vector <int> trueClass;//308

int * shift;
double ** expre;

int main(int argc, char **argv)
{
	int j,k;
	double loglike;
	vector <double> data;
	int column, row, length;
	list <int> tempID;
	clusterClass tempClust;
	list < clusterClass >::iterator allClustsIter;
	list < int >::iterator allGeneIDIter;
	long seed1, seed2;
	//int maxShift = 0;//1;
	list <clusterClass> finalClust,midClust;
	int *finalShift, *midShift;
	bool *finalSign, *midSign;
	double *postProba;	
	double sum, sum2;
	int count =0,misscount=0;
	double mean, var;
	double maxlike =0;
	list < int >::iterator geneIDIter;
	int updatedClustNumber;
	double priorbeta0;
	double priora;
	double priorb;
	bool allmiss;
        char inputfile[MAX_LENGTH]="inp.txt";
	char outputfile[MAX_LENGTH]="out.txt";
        char resultfile[MAX_LENGTH]="result.txt";
        char lengthfile[MAX_LENGTH]="length.txt";
        char memberfile[MAX_LENGTH]="member.txt";
		char standfile[MAX_LENGTH] = "stand.txt";
	int CHAINS =1;
	int ROUND = 10;
//	int BURN_IN = ROUND/2;
	bool invert = 1;
	int maxShift = 0;
	double probcutoff = 0;
	double strengthcutoff =0;
	int INITIAL_CLUST = 15;
	int **occur;
//	ofstream normalFile;
/*
	normalFile.open(standfile);
	if(!normalFile)
	{
		cout << "ERROR: Unable to open file: " << standfile << endl;
	    exit(20);
	}//end of if
*/
	srand((unsigned) time(NULL));
	seed1=rand()-2;//1;
	seed2=rand()+1;//2;
	setall(seed1,seed2);

// simulate data
//	simudata();
//	exit(0);
        if(argc<8)
        {
          printf("7 options need to be specified:\n\tinput file,\n\toutput cluster result,\n\tnumber of chains,\n\tnumber of cycles,\n\tinvert switch (0-no, 1-yes),\n\tmax shift,\n\tprob_cutoff.\n");
            exit(0);
        }  
        for(j=0;j<MAX_LENGTH;j++)
        {
            inputfile[j]=argv[1][j];
            outputfile[j]=argv[2][j];
//            lengthfile[j]=argv[3][j];
//            memberfile[j]=argv[4][j];
        }       
	CHAINS = atoi(argv[3]);
	ROUND = atoi(argv[4]);
	invert = atoi(argv[5]);
	maxShift = atoi(argv[6]);
	probcutoff = atof(argv[7]);
//	strengthcutoff = atof(argv[8]);
//	BURN_IN = ROUND/2;
// read in data
	cout <<"read in data."<<endl;
	readData(inputfile,data, &column);//yeast.txt//gal205.all.txt//tra.txt
	row = data.size()/column;
	//length = column - maxShift;
// standardize data
	//standard(data, column, row);
/*	for(j=0;j<row;j++)
	{
		for(k=0;k<column;k++)
		{
			normalFile << data[j*column + k]<<"	";
		}//end of k
		normalFile <<endl;
	}//end of j
	normalFile <<endl;
	normalFile.close();
*/
	expre = new double*[row];
	member = new list <clusterClass>::iterator[row];
	shift = new int[row];
//	initial shift position is in the middle
	for(j=0;j<row;j++)
		shift[j] = 0;//maxShift;
	if(invert)
	{
		sign = new bool[row];
		for(j=0;j<row;j++)
			sign[j] = true;
	}	
// initialize indicator for inversion
	if(diffe<1)
	{
		count =0;
		for(j=0;j<row;j++)
		{
			allmiss = true;
			for(k=0;k<column;k++)
			{
				if(data[j*column+k]>-900)
				{
					allmiss = false;
					break;
				}//end of if
			}//end of k
			if(allmiss)
			{
				cout << "the "<<j<<"th gene is completely missing."<<endl;
				count ++;
				continue;
			}
			expre[j] = new double[column];
			for(k=0;k<column;k++)
				expre[j][k] = data[j*column+k];
		}//end of j
	}//end of if
	else
	{
		count =0;
		for(j=0;j<row;j++)
		{
			allmiss = true;
			for(k=0;k<column-1;k++)
			{
				if((data[j*column+k]>-900)&&(data[j*column+k+1]>-900))
				{
					allmiss = false;
					break;
				}
			}//end of k
			if(allmiss)
			{
				cout << "the "<<j<<"th gene is completely missing when calculate difference."<<endl;
				count ++;
				continue;
			}
			expre[j] = new double[column-1];
			for(k=0;k<column-1;k++)
				expre[j][k] = data[j*column+k+1]-data[j*column+k];
		}//end of j
	}//end of else
	row = row -count;
	INITIAL_CLUST = (int) sqrt((double)row);
// calculate average for each gene, prepare to try inverted expression profiles
	if(invert)
	{
		sum =0;
		count =0;
		for(j=0;j<row;j++)
		{
			for(k=0;k<column;k++)
			{
				if(data[j*column+k]>-900)
				{
					sum = sum + data[j*column+k];
					count ++;
				}//end of if
			}//end of k
		}//end of j
		average = new double[row];
		for(j=0;j<row;j++)
		{
			average[j] = sum/(double) count;//0;
		}//end of j
	}//end of if
	if(diffe<1)
		length = column - maxShift;//2*maxShift
	else
		length = column-1 - maxShift;//2*maxShift
// calculate mean and variance, use them to choose parameters for priors
	sum = 0;
	sum2 = 0;
	count =0;
	missing = new bool[row];
	for(j=0;j<(row*column);j++)
	{
		//double a = data[j];
		double b = data[j];
		if(data[j]<-900)
		{
			count++;
			missing[j/column] = true;
		}
		else
		{
			missing[j/column] = false;
			sum = sum + data[j];
			sum2 = sum2 + data[j]*data[j];
		}//end of else
	}//end of j
	misscount = count;
	mean = sum/(double) (row*column - count);
	var = (sum2 - sum*mean)/((double) row*column - count -1);
//	cout << "mean= "<<mean <<" var= "<<var<<endl;//-02/16/06
	//exit(0);
	priorbeta0 = mean;
	//priora = 1;//2 + var;//0.5;
	if((invert)||(maxShift>0))//added 01/02/06
	{
		priora = 0.5;//2+var;
		priorb = var;//var*(1+var);
	}
	else
	{
		priora = 1;
		priorb = 2*var;//var+var*var;//var*0.5
	}
	//-02/16/06 cout <<priorbeta0<<" "<<priora<<" "<<priorb<<endl;
	//exit(0);
	data.clear();
/*	for(j=0;j<row;j++)
	{
		for(k=0;k<column;k++)
			cout <<expre[j][k]<<" ";
		cout <<endl;
	}
*/	//exit(0);
//	initialize clusters
// ********************** special initial for simulation.
/*	for(j=0;j<length;j++)
	{
		tempClust.ysum.push_back(0);
		tempClust.ysum2.push_back(0);
	}
	for(j=0;j<20;j++)
	{
		tempClust.geneID.push_back(j);
		for(k=0;k<length;k++)//column
		{
			if(expre[j][k]> -900)// test if it is missing, added 12/11/05//k + maxShift
			{		
				tempClust.ysum[k] = tempClust.ysum[k] + expre[j][k];
				tempClust.ysum2[k] = tempClust.ysum2[k] + expre[j][k]*expre[j][k];
			}//end of if
		}//end of k
	}//end of j
	allClusts.push_back(tempClust);
	tempClust.geneID.clear();
	tempClust.ysum.clear();
	tempClust.ysum2.clear();
	for(j=0;j<length;j++)
	{
		tempClust.ysum.push_back(0);
		tempClust.ysum2.push_back(0);
	}
	for(j=20;j<40;j++)
	{
		tempClust.geneID.push_back(j);
		for(k=0;k<length;k++)//column
		{
			if(expre[j][k]> -900)// test if it is missing, added 12/11/05//k + maxShift
			{		
				tempClust.ysum[k] = tempClust.ysum[k] + expre[j][k];
				tempClust.ysum2[k] = tempClust.ysum2[k] + expre[j][k]*expre[j][k];
			}//end of if
		}//end of k
	}//end of j
	allClusts.push_back(tempClust);
	tempClust.geneID.clear();
	tempClust.ysum.clear();
	tempClust.ysum2.clear();
// assign members
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(allGeneIDIter = (*allClustsIter).geneID.begin();allGeneIDIter != 
			(*allClustsIter).geneID.end(); allGeneIDIter++)
			member[(*allGeneIDIter)] = 	allClustsIter;	
	}
*/
// **************special initial above stop here.
// sequentially added each of the first (INITIAL_CLUST) genes to a cluster. 
	cout <<"start clustering."<<endl;
	for(j=0;j<INITIAL_CLUST;j++)
	{
		tempClust.geneID.push_back(j);
//		tempClust.ysum = new double[length];//column
//		tempClust.ysum2 = new double[length];//column
//		tempClust = clusterClass(j, length);
		for(k=0;k<length;k++)//column
		{
			if(expre[j][k]> -900)// test if it is missing, added 12/11/05//k + maxShift
			{
				//tempClust.ysum[k] = expre[j][k + maxShift];//k
				//tempClust.ysum2[k] = expre[j][k + maxShift]*expre[j][k + maxShift];//k
				tempClust.ysum.push_back(expre[j][k]);//k + maxShift]);
				tempClust.ysum2.push_back(expre[j][k]*expre[j][k]);//k + maxShift
			}//end of if
			else
			{
				tempClust.ysum.push_back(0);
				tempClust.ysum2.push_back(0); 
			}//end of else
		}//end of k
		allClusts.push_back(tempClust);
		tempClust.geneID.clear();
		tempClust.ysum.clear();
		tempClust.ysum2.clear();
	}
// sequentially add the remaining genes one by one to these clusters.
	for(j=INITIAL_CLUST;j<row;j++)
	{
		if ((j % INITIAL_CLUST) == 0)// set allClusts pointer to the beginning
			allClustsIter = allClusts.begin();
		else//add genes to existing clusters
			allClustsIter ++;
		(*allClustsIter).geneID.push_back(j);
		for(k=0;k<length;k++)
		{
			if(expre[j][k]> -900)// test if it is missing, added 12/11/05 //[k + maxShift]
			{
				(*allClustsIter).ysum[k] = (*allClustsIter).ysum[k] + expre[j][k];//k + maxShift];//k
				(*allClustsIter).ysum2[k] = (*allClustsIter).ysum2[k] + expre[j][k]*expre[j][k];//k + maxShift
				double a = (*allClustsIter).ysum[0];
			}//end of if
		}//end of k
	}//end of j
// assign members
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(allGeneIDIter = (*allClustsIter).geneID.begin();allGeneIDIter != 
			(*allClustsIter).geneID.end(); allGeneIDIter++)
			member[(*allGeneIDIter)] = 	allClustsIter;	
	}
/*
// output initial cluster information
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(allGeneIDIter = (*allClustsIter).geneID.begin();allGeneIDIter != 
			(*allClustsIter).geneID.end(); allGeneIDIter++)
		{
			cout << (*allGeneIDIter)<< " ";	
		}
		cout <<endl;
		cout <<(*allClustsIter).ysum[0]<< " "<<(*allClustsIter).ysum2[0]<<endl;
	}
	exit(0);
*/
// initialize co-occurance matrix
	occur = new int*[row];
	for(j=0;j<row;j++)
	{
		occur[j] = new int[row];
		for(k=0;k<row;k++)
			occur[j][k] = 0;
	}//end of j
// initial runs to determine likely number of clusters.
	maxlike = clust(row,length,maxShift,finalClust,finalShift,finalSign, priorbeta0,priora,priorb,ROUND,invert,occur);//column
	updatedClustNumber = finalClust.size();
	//-cout <<"preliminary run finished, there are "<<updatedClustNumber<<" clusters\n";
	//exit(0);
// fine search for mostly likely clustering using multiple chains
	for(j=0;j<CHAINS;j++)
	{
		cout <<"CHAIN " <<j+1<<endl;
		randomInitialClust(row,length,maxShift,INITIAL_CLUST,invert);//updatedClustNumber+j-CHAINS/2);
		loglike = clust(row,length,maxShift,midClust,midShift,midSign, priorbeta0,priora,priorb,ROUND,invert,occur);
		if(maxlike<loglike)
		{
			maxlike = loglike;
			finalShift = new int[row];
			finalClust.clear();
			finalClust.insert(finalClust.end(),midClust.begin(),midClust.end());
//			midClust.clear();
			for(k=0;k<row;k++)
				finalShift[k] = midShift[k];
			if(invert)
			{
				finalSign = new bool[row];
				for(k=0;k<row;k++)
					finalSign[k] =	midSign[k];
			}//end of if
		}//end of if
		midClust.clear();
		delete [] midShift;
		if(invert)
			delete [] midSign;
	}//end of j
	postProb(row, length,maxShift,finalClust,finalShift,finalSign,postProba,priorbeta0,priora,priorb,invert);
// output co-occurance matrix
/*	for(j=0;j<10;j++)
	{
		for(k=0;k<j;k++)
		{
			cout << occur[k][j]<<" ";
		}
		cout <<endl;
	}
	cout <<endl;
	//exit(0);
*/
	//oldpostProb(row, length,maxShift,finalClust,finalShift,postProba,priorbeta0,priora,priorb);
//	cout <<finalClust.size()<<" "<<finalShift[0]<<endl;
//	cout << postProba[0]<<" "<<postProba[1]<<endl;
	output(outputfile,row,column,length, maxlike,finalClust,finalShift,finalSign,
		 postProba,priorbeta0,priora,priorb,memberfile,lengthfile,resultfile,invert,
		 probcutoff,CHAINS, ROUND, misscount,occur);//column
	allClusts.clear();
	delete [] finalShift;
	for(j=0;j<row;j++)
		delete [] expre[j];
	delete [] expre;
	delete postProba;
	finalClust.clear();
	if(invert)
	{
		delete [] average;
		delete [] finalSign;
	}//end of if
	return(0);
}//end of main

void randomInitialClust(const int row, const int length, const int maxShift, const int updatedClustNumber,const int invert)
{
	int j,k,m;
	int na,nb,nc;
	double ran;
	int initialClust;
	int size;
	vector <int> *geneList;
	int order;
	clusterClass tempClust;
	list < clusterClass >::iterator allClustsIter;
	list < int >::iterator allGeneIDIter;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=50;

//	initial shift position is in the middle
	geneList = new vector <int>[updatedClustNumber]; 
	member = new list <clusterClass>::iterator[row];
	shift = new int[row];
//	initial shift position is in the middle
	for(j=0;j<row;j++)
		shift[j] = 0;//maxShift;
	if(invert)
	{
		sign = new bool[row];
		for(j=0;j<row;j++)
			sign[j] = true;
	}
	for(j=0;j<row;j++)
	{
		ran = unran(&na, &nb, &nc);
		initialClust = (int) (ran*(double) updatedClustNumber);
		geneList[initialClust].push_back(j);
	}//end of j
	for(j=0;j<updatedClustNumber;j++)
	{
		size = geneList[j].size();
		//tempClust.ysum = new double[length];
		//tempClust.ysum2 = new double[length];
		for(m=0;m<length;m++)
		{
			tempClust.ysum.push_back(0);
			tempClust.ysum2.push_back(0);
		}//end of m	
		for(k=0;k<size;k++)
		{
			order = geneList[j][k];
			tempClust.geneID.push_back(order);
			for(m=0;m<length;m++)
			{
				if(expre[order][m]> -900)// test if it is missing, added 12/11/05 //[m + maxShift]
				{	
					tempClust.ysum[m] = tempClust.ysum[m] + expre[order][m];//[m + maxShift]
					tempClust.ysum2[m] = tempClust.ysum2[m] + expre[order][m]*expre[order][m];//[m + maxShift]
				}//end of if
			}//end of m
		}//end of k
		allClusts.push_back(tempClust);
		tempClust.geneID.clear();
		tempClust.ysum.clear();
		tempClust.ysum2.clear();
//		delete [] tempClust.ysum;
//		delete [] tempClust.ysum2;
	}//end of j
	for(j=0;j<updatedClustNumber;j++)
		geneList[j].clear();
	delete [] geneList;
// assign members
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(allGeneIDIter = (*allClustsIter).geneID.begin();allGeneIDIter != 
			(*allClustsIter).geneID.end(); allGeneIDIter++)
			member[(*allGeneIDIter)] = 	allClustsIter;	
	}

// output initial cluster information
// delete 02/10/06
/*	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(allGeneIDIter = (*allClustsIter).geneID.begin();allGeneIDIter != 
			(*allClustsIter).geneID.end(); allGeneIDIter++)
		{
			cout << (*allGeneIDIter)<< " ";	
		}
		cout <<endl;
		cout <<(*allClustsIter).ysum[0]<< " "<<(*allClustsIter).ysum2[0]<<endl;
	}
*/	//exit(0);
}//end of randomInitialClust

void output(const char outPutFileName[], const int row, const int column, const int length, const double loglike,
			list <clusterClass> finalClust, int *finalShift, bool *finalSign, double *postProba, 
			const double priorbeta0, const double priora, double priorb, const char memberFileName[],
                        const char lengthFileName[], const char resultFileName[],const int invert,
			const double probcutoff,const int CHAINS, const int ROUND, const int misscount, int **occur)
{
	list <clusterClass>::iterator finalClustsIter;
	list < int >::iterator geneIDIter,ageneIDIter, bgeneIDIter;
	ofstream outPutFile;
	ofstream clustFile;
//	ofstream lengthFile;
//	ofstream resultFile;
//	ofstream stableFile;
	int finalClustNumber;
	int count=0,counter=0;
	double value =0;
	double ratio;
	int *classes;
	int agene,bgene;
	double sum =0;
	int nsize=0;
	double stable;

	outPutFile.open(outPutFileName);
	if(!outPutFile)
	{
		cout << "ERROR: Unable to open file: " << outPutFileName << endl;
	    exit(30);
	}//end of if
	clustFile.open(memberFileName);//"member.txt"
	if(!clustFile)
	{
		cout << "ERROR: Unable to open file: member.txt" << endl;
	    exit(31);
	}//end of if
//	lengthFile.open(lengthFileName);//"length.txt"
/*	if(!lengthFile)
	{
		cout << "ERROR: Unable to open file: length.txt" << endl;
	    exit(32);
	}//end of if
	stableFile.open("stable.txt");//stableFileName);
	if(!stableFile)
	{
		cout << "ERROR: Unable to open file: stable.txt" << endl;
	    exit(33);
	}//end of if
*/
	finalClustNumber = finalClust.size();
// find out minimum and maximum cluster width
	outPutFile << "**************************************** \n";
	outPutFile << "*                                      * \n";
	outPutFile << "*       CRC 1.0 clustering Result      * \n";
	outPutFile << "*                                      * \n";
	outPutFile << "**************************************** \n\n";
	
	outPutFile << "Total number of genes: " << row <<endl;
	    cout << "\nTotal number of genes: " << row <<endl;
	outPutFile << "Total number of experiments: " << column <<endl;
	    cout << "\nTotal number of experiments: " << column <<endl;
        outPutFile << "Rate of missing data: " << 100*(double) misscount/(double)(row*column) <<"%"<<endl;
            cout << "\nRate of missing data: " << 100*(double) misscount/(double)(row*column) <<"%"<<endl;
	outPutFile << "\nThere are total of " << finalClustNumber << " clusters"<<endl;
    cout <<"\n There are total of of " << finalClustNumber << " clusters, see output file " << outPutFileName << " for details"<<endl<<endl;
	outPutFile << "alpha = " << FACTOR_0<<"  ";
        outPutFile << "posterior probability threshold = " << probcutoff<<endl;
//        outPutFile << "cluster strength threshold = " << strengthcutoff<<endl;
	outPutFile << "log likelihood = " << loglike <<endl<<endl;
    cout <<"output log likelihood = "<<loglike <<endl<<endl;//-   
    count = 0;
//	resultFile.open(resultFileName);//"result.txt"
/*	if(!resultFile)
	{
		cout << "ERROR: Unable to open file: result.txt" << endl;
	    exit(33);
	}//end of if
*/
	classes = new int[row];
	for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
	{
	     ratio = finalRatio(finalClustsIter,finalClust, finalShift, finalSign,length,priorbeta0,priora,priorb,invert);
//	     if(ratio > strengthcutoff)
//	     {
		outPutFile <<"cluster "<<count+1<<" size = " <<(*finalClustsIter).geneID.size()<<"\nlog Bayes ratio = "<<ratio<<endl;
		cout       <<"cluster "<<count+1<<" size = " <<(*finalClustsIter).geneID.size()<<"\nlog Bayes ratio = "<<ratio<<endl;
		counter = 0;
// get co-occurance matrix
		nsize = (*finalClustsIter).geneID.size();
		if(nsize == 1)
			stable =0;
		else
		{
			sum = 0;
			for(ageneIDIter = (*finalClustsIter).geneID.begin();ageneIDIter != (*finalClustsIter).geneID.end();ageneIDIter++)
			{
				int a = (*((*finalClustsIter).geneID.end()));
				agene = (*ageneIDIter);
				for(bgeneIDIter = (*finalClustsIter).geneID.begin();bgeneIDIter != (*finalClustsIter).geneID.end();bgeneIDIter++)
				{
					bgene = (*bgeneIDIter);
					if(agene > bgene)
						sum = sum + occur[bgene][agene];
					else if(agene < bgene)
						sum = sum + occur[agene][bgene];
				}	
			}
			stable = sum/((double) 2*nsize*(nsize-1)*(CHAINS+1)*(ROUND-(ROUND/2)));
		}
		cout <<"average co-occurrences = "<<stable<<endl;
		outPutFile <<"average co-occurrences =  "<<stable<<endl;
//		stableFile <<stable<<" "<<ratio<<endl;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			classes[(*geneIDIter)] = count +1;
			//int a = (*geneIDIter);
			//double b = postProba[(*geneIDIter)];
			//int c= finalShift[(*geneIDIter)];
			
//			outPutFile <<(*geneIDIter)<<" "<<geneName[(*geneIDIter)]<<" "<<trueClass[(*geneIDIter)]<<" ("<<finalShift[(*geneIDIter)]<<") "<<" ["<<postProba[(*geneIDIter)]<<"] "<<endl;
//			cout       <<(*geneIDIter)<<" "<<geneName[(*geneIDIter)]<<" "<<trueClass[(*geneIDIter)]<<" ("<<finalShift[(*geneIDIter)]<<") "<<" ["<<postProba[(*geneIDIter)]<<"] "<<endl;
//			clustFile <<(*geneIDIter)<<" ";
            if(postProba[(*geneIDIter)]>probcutoff)
			{
				clustFile <<count+1<<"\t"<<(*geneIDIter)+1<<"\t"<<geneName[(*geneIDIter)]<<"\t";
				if((invert)&&(!finalSign[(*geneIDIter)]))
				{
					outPutFile <<(*geneIDIter)+1<<" "<<geneName[(*geneIDIter)]<<" (-"<<finalShift[(*geneIDIter)]<<") "<<"["<<postProba[(*geneIDIter)]<<"] ";
					cout       <<(*geneIDIter)+1<<" "<<geneName[(*geneIDIter)]<<" (-"<<finalShift[(*geneIDIter)]<<") "<<"["<<postProba[(*geneIDIter)]<<"] ";
					if(missing[(*geneIDIter)])
					{
						outPutFile <<"*"<<endl;
						cout <<"*"<<endl;
					}
					else
					{
						outPutFile <<endl;
						cout <<endl;
					}
					if(finalShift[(*geneIDIter)]>0)
						clustFile <<"4"<<endl;
					else
						clustFile <<"3"<<endl;
				}
				else
				{
					outPutFile <<(*geneIDIter)+1<<" "<<geneName[(*geneIDIter)]<<" (+"<<finalShift[(*geneIDIter)]<<") "<<"["<<postProba[(*geneIDIter)]<<"] ";
					cout       <<(*geneIDIter)+1<<" "<<geneName[(*geneIDIter)]<<" (+"<<finalShift[(*geneIDIter)]<<") "<<"["<<postProba[(*geneIDIter)]<<"] ";
					if(missing[(*geneIDIter)])
					{
						outPutFile <<"*"<<endl;
						cout <<"*"<<endl;
					}
					else
					{
						outPutFile <<endl;
						cout <<endl;
					}
					//clustFile <<(*geneIDIter)<<" ";
					if(finalShift[(*geneIDIter)]>0)
						clustFile <<"2"<<endl;
					else
						clustFile <<"1"<<endl;
				}
				counter ++;
            }//end of if
		}
		outPutFile <<"number of genes above threshold = " <<counter<<endl;
                cout       <<"number of genes above threshold = " <<counter<<endl;
		outPutFile <<endl;
		cout <<endl;
		//clustFile <<endl;
//		lengthFile << counter /*(*finalClustsIter).geneID.size()*/<<" ";
		count ++;
/* due to removal of strength cutoff 04/21/06
 	     }//end of if
	     else
	     {
			  outPutFile <<"cluster "<<count<<" size = " <<(*finalClustsIter).geneID.size()<<" Bayes ratio = "<<ratio<<", below cutoff"<<endl<<endl;
                 cout       <<"cluster "<<count<<" size = " <<(*finalClustsIter).geneID.size()<<" Bayes ratio = "<<ratio<<", below cutoff"<<endl<<endl;
	     }
due to removal of strength cutoff 04/21/06 */
	}
//	lengthFile <<endl;
//	for(j=0;j<row;j++)
//		resultFile <<classes[j]<<" ";
//	resultFile <<endl;
	outPutFile.close();
	clustFile.close();
//	lengthFile.close();
//	resultFile.close();
//	stableFile.close();
}//end of output

void postProb(const int row, const int length, const int maxShift, list <clusterClass> finalClust,
	      	int *finalShift, bool *finalSign, double *&postProba, const double priorbeta0, 
		const double priora, double priorb,const int invert)
{
	int k;
	list <clusterClass>::iterator finalClustsIter;
	list <clusterClass>::iterator allClustsIter;
	list < int >::iterator geneIDIter;
	list <clusterClass>::iterator genePlace;
	double *ysum, *ysum2;
	double *sum1,*sum2;
	int start, nsize,nsizenew;
	double ratio;
	double sum=0, win, numerator =0;

	ysum = new double[length];//column
	ysum2 = new double[length];//column
	sum1 = new double[length];//column
	sum2 = new double[length];//column
	postProba = new double[row];
/*	for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
	{
		//cout <<"*****"<<(*finalClustsIter).ysum[0]<<endl;
		double tempsum =0;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			cout << (*geneIDIter)<<" "<<finalShift[(*geneIDIter)]<<" ";
			if(invert)
				cout << finalSign[(*geneIDIter)]<<endl;
			tempsum = tempsum+ expre[(*geneIDIter)][finalShift[(*geneIDIter)]];
		}
		//cout <<"&&&&&" <<tempsum<<endl;
	}
*/	//exit(0);
	for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
	{
		for(k=0;k<length;k++)
		{
			sum1[k] = (*finalClustsIter).ysum[k];
			sum2[k] = (*finalClustsIter).ysum2[k];
		}//end of k
		nsize = (*finalClustsIter).geneID.size();//+1	
// test only		
/*		double temsum =0;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			cout << "order= "<<(*geneIDIter)<<endl;
			int a = finalShift[(*geneIDIter)];
			temsum = temsum+ expre[(*geneIDIter)][finalShift[(*geneIDIter)]];
//added exclusively for the simulated data 12/24/05		
			if((invert)&&(sign[(*geneIDIter)]))
				temsum = temsum - expre[(*geneIDIter)][finalShift[(*geneIDIter)]];
// stop here.
		}
		if((temsum - sum1[0])>0.001)
		{
			cout <<"error, temsum ="<<temsum<<" " <<sum1[0]<<endl;
			exit(5);
		}//end of if
*/// test only
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			start = finalShift[(*geneIDIter)];
			if((invert)&&(!finalSign[(*geneIDIter)]))
			{
				for(k=0;k<length;k++)//column
				{	
					if(expre[(*geneIDIter)][k + start]> -900)// test if it is missing, added 12/11/05
					{
						int a= *geneIDIter;
						ysum[k] = sum1[k] - 2*average[(*geneIDIter)] + expre[(*geneIDIter)][k + start];
						ysum2[k] = sum2[k] - (2*average[(*geneIDIter)] - expre[(*geneIDIter)][k + start])*(2*average[(*geneIDIter)] 
							- expre[(*geneIDIter)][k + start]);
					}//end of if
					else
					{
						ysum[k] = sum1[k];
						ysum2[k] = sum2[k];
					}//end of else
				}//end of k
				numerator = exp(bayesratioinvert((*geneIDIter),start,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
			}//end of if
			else
			{
				for(k=0;k<length;k++)//column
				{	
					if(expre[(*geneIDIter)][k + start]> -900)// test if it is missing, added 12/11/05
					{
						ysum[k] = sum1[k] - expre[(*geneIDIter)][k + start];
						ysum2[k] = sum2[k] - expre[(*geneIDIter)][k + start]*expre[(*geneIDIter)][k + start];
					}//end of if
					else
					{
						ysum[k] = sum1[k];
						ysum2[k] = sum2[k];
					}//end of else
				}//end of k
				numerator = exp(bayesratio((*geneIDIter),start,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
			}//end of else	
			sum=0;
			for(allClustsIter = finalClust.begin();allClustsIter != finalClust.end(); allClustsIter++)
			{
				if(allClustsIter != finalClustsIter)//finalClustsIter is not the cluster that contain (*geneIDIter)
				{
					for(k=0;k<length;k++)
					{
						ysum[k] = (*allClustsIter).ysum[k];
						ysum2[k] = (*allClustsIter).ysum2[k];
					}//end of k
// find the largest *******************12/15/05
					nsizenew = (*allClustsIter).geneID.size() + 1;	
					for(k=0;k<(maxShift + 1);k++)//(2*maxShift + 1)
					{	
						ratio = exp(bayesratio((*geneIDIter),k,length,ysum,ysum2,nsizenew,priorbeta0,priora,priorb));
						if((k==0)||(win<ratio))
							win = ratio;
					}//end of k		
					if(invert)
					{
						for(k=0;k<(maxShift + 1);k++)//(2*maxShift + 1)
						{	
							ratio = exp(bayesratioinvert((*geneIDIter),k,length,ysum,ysum2,nsizenew,priorbeta0,priora,priorb));
							if(win<ratio)
								win = ratio;
						}//end of k
					}//end of if
					sum = sum + win;
				}//end of if
			}//allClustsIter
			//int a = (*geneIDIter);
			postProba[(*geneIDIter)] = numerator/(numerator + sum);
			//-02/16/06 cout << numerator << " " <<sum<<endl;
		}//geneIDIter
	}//finalClustsIter
	delete [] sum1;
	delete [] sum2;
	delete [] ysum;
	delete [] ysum2;		
}//end of postProb

void oldpostProb(const int row, const int length, const int maxShift, list <clusterClass> finalClust,
			  int *finalShift, double *&postProba, const double priorbeta0, const double priora, double priorb)
{
	int k;
	list <clusterClass>::iterator finalClustsIter;
	list <clusterClass>::iterator allClustsIter;
	list < int >::iterator geneIDIter;
	list <clusterClass>::iterator genePlace;
	double *ysum, *ysum2;
	double *sum1,*sum2;
	int start, nsize,nsizenew;
	double ratio;
	double sum=0,insum=0, numerator =0;

	ysum = new double[length];//column
	ysum2 = new double[length];//column
	sum1 = new double[length];//column
	sum2 = new double[length];//column
	postProba = new double[row];
/*	for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
	{
		cout <<"*****"<<(*finalClustsIter).ysum[0]<<endl;
		double tempsum =0;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			//cout << (*geneIDIter)<<endl;
			tempsum = tempsum+ expre[(*geneIDIter)][finalShift[(*geneIDIter)]];
		}
		cout <<"&&&&&" <<tempsum<<endl;
	}
*/	//exit(0);
	for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
	{
		for(k=0;k<length;k++)
		{
			sum1[k] = (*finalClustsIter).ysum[k];
			sum2[k] = (*finalClustsIter).ysum2[k];
		}//end of k
		nsize = (*finalClustsIter).geneID.size();//+1	
		
		double temsum =0;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			//cout << "order= "<<(*geneIDIter)<<endl;
			int a = finalShift[(*geneIDIter)];
			temsum = temsum+ expre[(*geneIDIter)][finalShift[(*geneIDIter)]];
		}
		if((temsum - sum1[0])>0.001)
		{
			cout <<"error, temsum ="<<temsum<<" " <<sum1[0]<<endl;
			exit(5);
		}//end of if

		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			start = finalShift[(*geneIDIter)];
			for(k=0;k<length;k++)//column
			{	
				if(expre[(*geneIDIter)][k + start]> -900)// test if it is missing, added 12/11/05
				{
					ysum[k] = sum1[k] - expre[(*geneIDIter)][k + start];
					ysum2[k] = sum2[k] - expre[(*geneIDIter)][k + start]*expre[(*geneIDIter)][k + start];
				}//end of if
				else
				{
					ysum[k] = sum1[k];
					ysum2[k] = sum2[k];
				}//end of else
			}//end of k
			insum =0;
			for(k=0;k<(maxShift +1);k++)//(2*maxShift + 1)
			{	
				ratio = exp(bayesratio((*geneIDIter),k,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
				//cout <<k<<"  "<<ratio;
				insum = insum + ratio;
			}//end of k
			sum=0;
			for(allClustsIter = finalClust.begin();allClustsIter != finalClust.end(); allClustsIter++)
			{
				if(allClustsIter != finalClustsIter)//finalClustsIter is the member of (*geneIDIter)
				{
					for(k=0;k<length;k++)
					{
						ysum[k] = (*allClustsIter).ysum[k];
						ysum2[k] = (*allClustsIter).ysum2[k];
					}//end of k

		double temsum =0;
		list < int >::iterator IDIter;
		for(IDIter = (*allClustsIter).geneID.begin();IDIter != (*allClustsIter).geneID.end();IDIter++)
		{
			//cout << "order 02 = "<<(*IDIter)<<endl;
			int b = finalShift[(*IDIter)];
			temsum = temsum+ expre[(*IDIter)][finalShift[(*IDIter)]];
		}
		if((temsum - ysum[0])>0.001)
		{
			cout <<"error, temsum ="<<temsum<<" " <<sum1[0]<<endl;
			exit(5);
		}//end of if

					nsizenew = (*allClustsIter).geneID.size() + 1;
					for(k=0;k<(maxShift +1);k++)//(2*maxShift + 1)
					{	
						ratio = exp(bayesratio((*geneIDIter),k,length,ysum,ysum2,nsizenew,priorbeta0,priora,priorb));
						//cout <<k<<"  "<<ratio;
						sum = sum + ratio;
					}//end of k		
				}//end of if
			}//allClustsIter
			int a = (*geneIDIter);
			postProba[(*geneIDIter)] = insum/(insum+sum);
			//cout << insum << " " <<sum<<endl;
		}//geneIDIter
	}//finalClustsIter
	/*
	for(j=0;j<row;j++)
	{
		sum =0;
		genePlace = finalGenes[j].member;
		for(finalClustsIter = finalClust.begin();finalClustsIter != finalClust.end(); finalClustsIter++)
		{
			for(k=0;k<length;k++)
			{
				ysum[k] = (*finalClustsIter).ysum[k];
				ysum2[k] = (*finalClustsIter).ysum2[k];
			}//end of k
			nsize = (*finalClustsIter).geneID.size();
			if(genePlace == finalClustsIter)
				self = true;
			else
				self = false;
			if(self)// cluster contain this gene
			{
				start = finalGenes[j].shift;
				for(k=0;k<length;k++)//column
				{		
					ysum[k] = ysum[k] - expre[j][k + start];
					ysum2[k] = ysum2[k] - expre[j][k + start]*expre[j][k + start];
				}//end of k
			}//end of if
			else 
				nsize ++;
			insum =0;
			for(k=0;k<(maxShift + 1);k++)//(2*maxShift + 1)
			{	
				ratio = exp(bayesratio(j,k,length,ysum,ysum2,nsize));
				insum = insum + ratio;
			}//end of j
			sum = sum + insum;
			if(self)
				numerator = insum;
		}
		postProba[j] = numerator/sum;		
	}//end of j
	*/
	delete [] sum1;
	delete [] sum2;
	delete[] ysum;
	delete[] ysum2;		
}//end of oldpostProb

double finalRatio(list <clusterClass>::iterator finalClustsIter, list <clusterClass> finalClust,
		      int *finalShift, bool *finalSign, const int length, const double priorbeta0, 
			const double priora, double priorb,const int invert)
{
	int j,nsize;
	list < int >::iterator geneIDIter;
	double ysum, ysum2,expression;
	double togetherloglike,seperateloglike;
	double loglike, dif, constant;

	nsize = (*finalClustsIter).geneID.size(); 
// test only
/*	togetherloglike = 0;
	seperateloglike = 0;
	for(j=0;j<length;j++)
	{
		ysum = (*finalClustsIter).ysum[j];
		ysum2 = (*finalClustsIter).ysum2[j];
		loglike = gammaln(0.5*(double) nsize+priora) - sqrt((double) nsize + 1)
			- (0.5*(double)nsize + priora)*log(priorb + 0.5*(ysum2 + priorbeta0*priorbeta0 
			- (ysum + priorbeta0)*(ysum + priorbeta0)/((double) nsize+1)));
		togetherloglike = togetherloglike + loglike;
		for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
		{
			expression = expre[(*geneIDIter)][j + finalShift[(*geneIDIter)]];
			loglike = - (0.5 + priora)*log(priorb + 0.5*(expression - priorbeta0)*(expression - priorbeta0));
			seperateloglike = seperateloglike + loglike;
		}	
	}//end of j
*/
	if(nsize == 1)
	{
		return 0;
	}
	else
	{
		togetherloglike = 0;
		seperateloglike = 0;
		for(j=0;j<length;j++)
		{
			ysum = (*finalClustsIter).ysum[j];
			ysum2 = (*finalClustsIter).ysum2[j];
			loglike = gammaln(0.5*(double) nsize+priora) - 0.5*log((double) nsize + 1)
				- (0.5*(double)nsize + priora)*log(priorb + 0.5*(ysum2 + priorbeta0*priorbeta0 
				- (ysum + priorbeta0)*(ysum + priorbeta0)/((double) nsize+1)));
			togetherloglike = togetherloglike + loglike;
			for(geneIDIter = (*finalClustsIter).geneID.begin();geneIDIter != (*finalClustsIter).geneID.end();geneIDIter++)
			{
				expression = expre[(*geneIDIter)][j + finalShift[(*geneIDIter)]];
				if(expression> -900)// test if it is missing, added 12/11/05
				{
					if((invert)&&(!expre[(*geneIDIter)]))
						expression = 2*average[(*geneIDIter)] - expression;
					loglike = - (0.5 + priora)*log(priorb + 0.5*(expression - priorbeta0)*(expression - priorbeta0));	
				}
				seperateloglike = seperateloglike + loglike;
			}
		}//end of j
		constant = (double) length*nsize*(gammaln(0.5 + priora) - 0.5*log(2.0)) + (double) length*(nsize-1)*(priora*log(priorb)-gammaln(priora));
		dif = togetherloglike - seperateloglike - constant;
	}//end of else
	return(dif/(double)nsize);
}//finalratio

double loglikelihood(const int length, const double priorbeta0, const double priora, double priorb)//column
{
	int j;
	double dsize;
	double ysum, ysum2;
	list < clusterClass >::iterator allClustsIter;
	double loglike,loglikesum;
	int count;

	loglikesum = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		dsize = (double) (*allClustsIter).geneID.size();
		for(j=0;j<length;j++)//column
		{
			ysum = (*allClustsIter).ysum[j];
			ysum2 = (*allClustsIter).ysum2[j];
			loglike = - (0.5*dsize + priora)*log(priorb 
				+ 0.5*(ysum2 + priorbeta0*priorbeta0 - (ysum + priorbeta0)*(ysum + priorbeta0)/(dsize+1)));
			loglikesum = loglikesum + loglike;
		}//end of j
		loglikesum = loglikesum + (double)length*(-0.5*dsize*log(2*pie) + gammaln(0.5*dsize+priora) - 0.5*log(dsize + 1));
	}
// the next term appears in every column of every cluster. 
	count = allClusts.size();
	double a = priora*log(priorb) - gammaln(priora);
	loglikesum = loglikesum + (double) count*length*(priora*log(priorb) - gammaln(priora));
	return(loglikesum);
}//end of loglikelihood

double clust(const int row, const int length, const int maxShift,list <clusterClass> &finalClust,
			 int *&finalShift, bool *&finalSign, const double priorbeta0, 
			 const double priora, double priorb, const int ROUND,const int invert, int **&occur)//column
{
	int j,k,m;
	double large;
	double loglike;
	ofstream outputFile;
	clusterClass tempClust;
	list < clusterClass >::iterator allClustsIter;
	list < int >::iterator geneIDIter;

	outputFile.open("like.txt");
	if(!outputFile)
	{
		cout << "ERROR: Unable to open file: like.txt" << endl;
	    exit(32);
	}//end of if	
// start work
	finalShift = new int[row];
	loglike = loglikelihood(length,priorbeta0,priora,priorb);//
	outputFile <<loglike<<endl;
// ******** remove temporaryly
	//cout <<"log likelihood =" <<loglike<<"*****"<<allClusts.size()<<endl;
// ******** remove temporaryly
	large = loglike;
//initialize finalClust by the current cotents on allClusts
	finalClust.clear();
	finalClust.insert(finalClust.end(),allClusts.begin(),allClusts.end());
	for(m=0;m<row;m++)
		finalShift[m] = shift[m];
	if(invert)
	{
		finalSign = new bool[row];
		for(m=0;m<row;m++)
			finalSign[m] = sign[m];
	}//end of if
	for(j=0;j<ROUND;j++)
	{
// test agreement
/*		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		{
//-			cout <<"$$$$$"<<(*allClustsIter).ysum[0]<<endl;
			double tempsum1 =0;
			cout << (*allClustsIter).geneID.size()<<"  ";
			for(geneIDIter = (*allClustsIter).geneID.begin();geneIDIter != (*allClustsIter).geneID.end();geneIDIter++)
			{
				cout << (*geneIDIter)<<" ";
				if((invert)&&(!sign[(*geneIDIter)]))
					tempsum1 = tempsum1 + 2*average[(*geneIDIter)] - expre[(*geneIDIter)][shift[(*geneIDIter)]];
				else
					tempsum1 = tempsum1 + expre[(*geneIDIter)][shift[(*geneIDIter)]];
			}
			cout << endl;
//-			cout <<"@@@@@" <<tempsum1<<endl;
		}
		cout <<endl;
*/
// test agreement	
		for(k=0;k<row;k++)
		{//************** delete later 
/*			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			{
				cout <<" size = " <<(*allClustsIter).geneID.size()<<endl;
				for(geneIDIter = (*allClustsIter).geneID.begin();geneIDIter != (*allClustsIter).geneID.end();geneIDIter++)
					cout <<(*geneIDIter)<<" ("<<shift[(*geneIDIter)]<<") ";
				cout <<endl<<endl;
			}
*/			//exit(0);
//****************** delete later
			if(invert)
				operationinvert(k,length,maxShift,priorbeta0,priora,priorb,invert);	
			else
				operation(k,length,maxShift,priorbeta0,priora,priorb);
			loglike = loglikelihood(length,priorbeta0,priora,priorb);//column
			outputFile <<loglike<<endl;
// ******** remove temporaryly
			//cout <<"log likelihood =" <<loglike<<"*****"<<allClusts.size()<<endl;
// ******** remove temporaryly
			/*if (loglike > large)
			{
				cout << "finalClust replaced"<< endl;
			}*/
			if (loglike > large)
			{
				large=loglike;
				finalClust.clear();
				finalClust.insert(finalClust.end(),allClusts.begin(),allClusts.end());
				for(m=0;m<row;m++)
					finalShift[m] = shift[m];
				if(invert)
				{
					for(m=0;m<row;m++)
						finalSign[m] =sign[m];
				}//end of if
			}//end of if
		}//end of k
// update co-occurance matrix
		if(j>= (ROUND/2))
			cooccur(row, occur);
// test agreement 
/*		for(allClustsIter = finalClust.begin();allClustsIter != finalClust.end(); allClustsIter++)
		{
			cout <<"%%%%%"<<(*allClustsIter).ysum[0]<<endl;
			double tempsum2 =0;
			for(geneIDIter = (*allClustsIter).geneID.begin();geneIDIter != (*allClustsIter).geneID.end();geneIDIter++)
			{
//				cout << (*geneIDIter)<<" " <<finalShift[(*geneIDIter)]<<" ";
//				if(invert)
//					cout <<finalSign[(*geneIDIter)]<<" ";
				if((invert)&&(!finalSign[(*geneIDIter)]))
					tempsum2 = tempsum2 + 2*average[(*geneIDIter)] - expre[(*geneIDIter)][shift[(*geneIDIter)]];
				else
				{
//					cout <<expre[(*geneIDIter)][shift[(*geneIDIter)]]<<endl;
					tempsum2 = tempsum2 + expre[(*geneIDIter)][shift[(*geneIDIter)]];
				}
			}
			cout <<"#####" <<tempsum2<<endl;
		}
		cout <<endl;
*/
// test agreement	
	}//end of j

// clean up
	delete [] member;
	delete [] shift;
	if(invert)//added 12/25/05
		delete [] sign;
	allClusts.clear();
	outputFile.close();
	return(large);//loglike
}//end of clust

void operationinvert(const int geneID, const int length, const int maxShift, const double priorbeta0, 
			const double priora, double priorb,const int invert)// const int column,
{
	int j;
	int outcome,newshift, oldstart;
	int count = 0;
	double ratio;
	double *ysum, *ysum2;
	int nsize,msize;
	bool self,single;
	vector <double> ratioVec;
	clusterClass tempClust;
	list < clusterClass >::iterator allClustsIter;
	list < int >::iterator geneIDIter;
	list <clusterClass>::iterator genePlace;
	list <clusterClass>::iterator temPlace;
	bool change = false;
	int start;
	double *temexpre;
	bool newsign, oldsign;
	bool outsign;

	oldstart = shift[geneID];
	oldsign = sign[geneID];
	ysum = new double[length];//column
	ysum2 = new double[length];//column
	genePlace = member[geneID];
	single = false;
	//remove the cluster if gene is in a singleton cluster
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if(genePlace == allClustsIter)
		{
			if((*allClustsIter).geneID.size() ==1)
			{
				single = true;
				break;
			}//end of if
		}//end of if
	}
	if(single)
		allClusts.erase(allClustsIter);
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(j=0;j<length;j++)//column
		{
			ysum[j] = (*allClustsIter).ysum[j];
			ysum2[j] = (*allClustsIter).ysum2[j];
		}//end of j
		nsize = (*allClustsIter).geneID.size();
		msize = nsize;
		if(genePlace == allClustsIter)//cluster contain this gene
		{
			self = true;
			start = shift[geneID];
// ***** add 12/15/05 for testing on inversion			
			if((invert)&&(!sign[geneID]))
			{
				for(j=0;j<length;j++)//column
				{
					if(expre[geneID][j + start]> -900)// test if it is missing, added 12/11/05
					{
						ysum[j] = ysum[j] - 2*average[geneID] + expre[geneID][j + start];
						ysum2[j] = ysum2[j] - (2*average[geneID] - expre[geneID][j + start])*(2*average[geneID] - expre[geneID][j + start]);
					}//end of if
				}//end of j	
			}//end of if
			else
			{
				for(j=0;j<length;j++)//column
				{
					if(expre[geneID][j + start]> -900)// test if it is missing, added 12/11/05
					{
						ysum[j] = ysum[j] - expre[geneID][j + start];
						ysum2[j] = ysum2[j] - expre[geneID][j + start]*expre[geneID][j + start];
					}//end of if
				}//end of j
			}//end of else
		}//end of if
		else
		{
			self = false;
			nsize ++;
		}//end of else
/*		
		double sum =0;
		for(geneIDIter = (*allClustsIter).geneID.begin();geneIDIter != (*allClustsIter).geneID.end();geneIDIter++)
		{
			//cout << (*geneIDIter)<<endl;
			sum = sum+ expre[(*geneIDIter)][0];
		}
		if((sum - ysum[0])>0.001)
		{
			cout <<"error"<<endl;
			exit(5);
		}//end of if
//		cout <<endl<<sum<<endl<<sum2<<endl;
*/
		for(j=0;j<(maxShift + 1);j++)//j<(2*maxShift + 1)
		{	
			ratio = exp(bayesratio(geneID,j,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
//			cout <<"ratio= "<<ratio<<endl;
			//ratio = ratio * (double)msize;//dirichlet prior 12/13/05
			ratioVec.push_back(ratio);
		}//end of j
// ***** add 12/15/05 for testing on inversion	
		for(j=0;j<(maxShift + 1);j++)//(2*maxShift + 1)
		{	
			ratio = exp(bayesratioinvert(geneID,j,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
//			cout <<"invratio= "<<ratio<<endl;
			//ratio = ratio * (double)msize;//dirichlet prior 12/13/05
			ratioVec.push_back(ratio);
		}//end of j
	}
	delete[] ysum;
	delete[] ysum2;
	if(maxShift==0)
	{
		outcome = singleCompareInvert(ratioVec,outsign,invert);
		newshift = 0;
		newsign = outsign;
	}
	else
	{
		outcome = compareInvert(ratioVec,maxShift,newshift,outsign,invert);
		newsign = outsign;
	}
	ratioVec.clear();
	if(outcome ==0)//new cluster is forming
	{
		//remove gene from previous cluster
		if(!single)
		{
			//count = 0;
			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			{
				if(genePlace == allClustsIter)//found cluster
				{
					(*allClustsIter).geneID.pop_front();
					//oldstart = shift[geneID];
					//oldsign = sign[geneID];
					temexpre = new double[length];
					if(sign[geneID])
					{
						for(j=0;j<length;j++)
							temexpre[j] = expre[geneID][j + oldstart]; 
					}
					else// if previously, it is inverted.
					{
						for(j=0;j<length;j++)
							temexpre[j] = 2* average[geneID] - expre[geneID][j + oldstart]; 
					}//end of else
					for(j=0;j<length;j++)
					{
						if(expre[geneID][j + oldstart]> -900)// test if it is missing, added 12/11/05
						{
							(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] - temexpre[j];
							(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] - temexpre[j]*temexpre[j];
						}//end of if
					}//end of j
					break;
					delete [] temexpre;
				}//end of if		
			//	else
			//		count ++;
			}
		}//end of if
		//create a new cluster that contain this gene only.
		tempClust.geneID.push_back(geneID);
		//tempClust.ysum = new double[length];//changed 121305
		//tempClust.ysum2 = new double[length];//changed 121305
		for(j=0;j<length;j++)
		{
			if(expre[geneID][j]> -900)// test if it is missing, added 12/11/05//j + maxShift
			{
				tempClust.ysum.push_back(expre[geneID][j]);//j + maxShift
				tempClust.ysum2.push_back(expre[geneID][j]*expre[geneID][j]);//j + maxShift
			}//end of if
			else
			{
				tempClust.ysum.push_back(0);
				tempClust.ysum2.push_back(0);
			}//end of else
		}//end of j
		allClusts.push_back(tempClust);
		tempClust.geneID.clear();
		tempClust.ysum.clear();
		tempClust.ysum2.clear();
		temPlace = allClusts.end();
		temPlace --;
		member[geneID] = temPlace;
		shift[geneID] = 0;//maxShift;
		sign[geneID] = true;
	}//end of if
	else //join one of exisiting clusters
	{
		if(!single)
		{
			change = false;
			count = 0;
			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			{
				if(genePlace == allClustsIter)//found the cluster the gene was in before operation
				{
					(*allClustsIter).geneID.pop_front();
					//oldstart = shift[geneID];
					if((outcome != (count + 1))||(oldstart != newshift)||(oldsign != newsign))//move from one cluster to another
					{
						change = true;
						temexpre = new double[length];
						if(sign[geneID])
						{
							for(j=0;j<length;j++)
								temexpre[j] = expre[geneID][j + oldstart]; 
						}
						else// if previously, it is inverted.
						{
							for(j=0;j<length;j++)
								temexpre[j] = 2* average[geneID] - expre[geneID][j + oldstart]; 
						}//end of else
						for(j=0;j<length;j++)
						{
							//double a = (*allClustsIter).ysum[j];
							if(expre[geneID][j + oldstart]> -900)// test if it is missing, added 12/11/05
							{
								(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] - temexpre[j];
								(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] - temexpre[j]*temexpre[j];
							}//end of if
						}//end of j
						delete [] temexpre;
					}//end of if
					break;
				}//end of if	
				else
					count ++;
			}//end of if
		}//end of if
// add this gene to another cluster
		count = 0;
		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		{
			if(outcome == (count+1))//found cluster the gene will join
			{
				(*allClustsIter).geneID.push_back(geneID);
				//if(genePlace != allClustsIter)//move from one cluster to another
				if((single)||(change))
				{
					temexpre = new double[length];
					if(newsign)
					{
						for(j=0;j<length;j++)
							temexpre[j] = expre[geneID][j + newshift]; 
						sign[geneID] = true;
					}
					else// if it is inverted.
					{
						for(j=0;j<length;j++)
							temexpre[j] = 2* average[geneID] - expre[geneID][j + newshift];
						sign[geneID] = false;
					}//end of else
					for(j=0;j<length;j++)
					{
						if(expre[geneID][j + newshift]> -900)// test if it is missing, added 12/11/05
						{
							(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] + temexpre[j];
							(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] + temexpre[j]*temexpre[j];
						}//end of if
					}//end of j
					delete [] temexpre;
					member[geneID] = allClustsIter;
					shift[geneID] = newshift;
				}//end of if
				break;
			}//end of if	
			else
				count ++;
		}
	}//end of else
}//end of operationinvert

void operation(const int geneID, const int length, const int maxShift, const double priorbeta0, const double priora, double priorb)// const int column,
{
	int j;
	int outcome,newshift, oldstart;
	int count = 0;
	double ratio;
	double *ysum, *ysum2;
	int nsize,msize;
	bool self,single;
	vector <double> ratioVec;
	clusterClass tempClust;
	list < clusterClass >::iterator allClustsIter;
	list < int >::iterator geneIDIter;
	list <clusterClass>::iterator genePlace;
	list <clusterClass>::iterator temPlace;
	bool change = false;
	int start;

	ysum = new double[length];//column
	ysum2 = new double[length];//column
	genePlace = member[geneID];
	single = false;
	//remove the cluster if gene is in a singleton cluster
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if(genePlace == allClustsIter)
		{
			if((*allClustsIter).geneID.size() ==1)
			{
				single = true;
				break;
			}//end of if
		}//end of if
	}
	if(single)
		allClusts.erase(allClustsIter);
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(j=0;j<length;j++)//column
		{
			ysum[j] = (*allClustsIter).ysum[j];
			ysum2[j] = (*allClustsIter).ysum2[j];
		}//end of j
		nsize = (*allClustsIter).geneID.size();
		msize = nsize;
		if(genePlace == allClustsIter)
			self = true;
		else
			self = false;
		if(self)// cluster contain this gene
		{
			start = shift[geneID];
			for(j=0;j<length;j++)//column
			{
				if(expre[geneID][j + start]> -900)// test if it is missing, added 12/11/05
				{
					ysum[j] = ysum[j] - expre[geneID][j + start];
					ysum2[j] = ysum2[j] - expre[geneID][j + start]*expre[geneID][j + start];
				}//end of if
			}//end of j
		}//end of if
		else 
			nsize ++;
//***** test only		
/*		double sum =0;
		double sum2=0;
		for(geneIDIter = (*allClustsIter).geneID.begin();geneIDIter != (*allClustsIter).geneID.end();geneIDIter++)
		{
			cout << (*geneIDIter)<<" ";
			sum = sum+ expre[(*geneIDIter)][start];
			sum2 = sum2 + expre[(*geneIDIter)][start]*expre[(*geneIDIter)][start];
		}
		if((sum - ysum[0])>0.001)
		{
			cout <<"error "<<sum <<" " <<ysum[0]<<endl;
			exit(5);
		}//end of if
*/
//		cout <<endl<<sum<<endl<<sum2<<endl;
//***** test only
		for(j=0;j<(maxShift + 1);j++)//j<(2*maxShift + 1)
		{	
			ratio = exp(bayesratio(geneID,j,length,ysum,ysum2,nsize,priorbeta0,priora,priorb));
//			cout <<"ratio= "<<ratio<<endl;
			//ratio = ratio * (double)msize;//dirichlet prior 12/13/05
			ratioVec.push_back(ratio);
		}//end of j
	}
	delete[] ysum;
	delete[] ysum2;
/*	count = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		cout << endl<<count <<" "<<(*allClustsIter).geneID.size()<<" ";
		for(j=0;j<length;j++)
			cout <<(*allClustsIter).ysum[j] <<" ";
		cout<<endl;
		count ++;
	}
*/
	if(maxShift==0)
	{
		outcome = singleCompare(ratioVec);
		newshift = 0;
	}
	else
		outcome = compare(ratioVec,maxShift,newshift);
	ratioVec.clear();
	if(outcome ==0)//new cluster is forming
	{
		//remove gene from previous cluster
		if(!single)
		{		
			count = 0;
			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			{
				if(genePlace == allClustsIter)//found cluster
				{
					(*allClustsIter).geneID.pop_front();
					oldstart = shift[geneID];
					for(j=0;j<length;j++)
					{
						if(expre[geneID][j + oldstart]> -900)// test if it is missing, added 12/11/05
						{
							(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] - expre[geneID][j + oldstart];
							(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] 
										- expre[geneID][j + oldstart]*expre[geneID][j + oldstart];
						}//end of if
					}//end of j
					break;
				}//end of if		
				else
					count ++;
			}
		}//end of if
		//create a new cluster that contain this gene only.
		tempClust.geneID.push_back(geneID);
		//tempClust.ysum = new double[length];//changed 121305
		//tempClust.ysum2 = new double[length];//changed 121305
		for(j=0;j<length;j++)
		{
			if(expre[geneID][j]> -900)// test if it is missing, added 12/11/05//j + maxShift
			{
				tempClust.ysum.push_back(expre[geneID][j]);//j + maxShift
				tempClust.ysum2.push_back(expre[geneID][j]*expre[geneID][j]);//j + maxShift
			}//end of if
			else
			{
				tempClust.ysum.push_back(0);
				tempClust.ysum2.push_back(0);
			}//end of else
		}//end of j
		allClusts.push_back(tempClust);
		tempClust.geneID.clear();
		tempClust.ysum.clear();
		tempClust.ysum2.clear();
		temPlace = allClusts.end();
		temPlace --;
		member[geneID] = temPlace;
		shift[geneID] = 0;//maxShift;
	}//end of if
	else //join one of exisiting clusters
	{
		if(!single)
		{
			change = false;
			count = 0;
			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			{
				if(genePlace == allClustsIter)//found the cluster the gene was in before operation
				{
					(*allClustsIter).geneID.pop_front();
					oldstart = shift[geneID];
					if((outcome != (count + 1))||(oldstart != newshift))//move from one cluster to another
					{
						change = true;
						for(j=0;j<length;j++)
						{
							//double a = (*allClustsIter).ysum[j];
							if(expre[geneID][j + oldstart]> -900)// test if it is missing, added 12/11/05
							{
								(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] - expre[geneID][j + oldstart];
								(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] 
									- expre[geneID][j + oldstart]*expre[geneID][j + oldstart];
							}//end of if
						}//end of j
					}//end of if
					break;
				}//end of if	
				else
					count ++;
			}
		}//end of if
// add this gene to another cluster
		count = 0;
		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		{
			if(outcome == (count+1))//found cluster the gene will join
			{
				(*allClustsIter).geneID.push_back(geneID);
				//if(genePlace != allClustsIter)//move from one cluster to another
				if((single)||(change))
				{
					for(j=0;j<length;j++)
					{
						if(expre[geneID][j + newshift]> -900)// test if it is missing, added 12/11/05
						{
							(*allClustsIter).ysum[j] = (*allClustsIter).ysum[j] + expre[geneID][j + newshift];
							(*allClustsIter).ysum2[j] = (*allClustsIter).ysum2[j] 
								+ expre[geneID][j + newshift]*expre[geneID][j + newshift];
						}//end of if
					}//end of j
					member[geneID] = allClustsIter;
					shift[geneID] = newshift;
				}//end of if
				break;
			}//end of if	
			else
				count ++;
		}
	}//end of else
// check if ysum match sums of all member genes	
/*
	count = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		cout << endl<<count <<" "<<(*allClustsIter).geneID.size()<<" ";
		for(j=0;j<length;j++)
			cout <<(*allClustsIter).ysum[j] <<" ";
		cout<<endl;
		double sum1 =0;
		for(geneIDIter = (*allClustsIter).geneID.begin(); geneIDIter != (*allClustsIter).geneID.end(); geneIDIter++)
		{
			//cout << (*geneIDIter)<< " ";
			sum1 = sum1 + expre[(*geneIDIter)][shift[(*geneIDIter)]];
		}
		//cout <<endl;
		double sum2 = (*allClustsIter).ysum[0];
		cout <<"sum="<<sum1<<" sum= " <<(*allClustsIter).ysum[0]<<"&&&&&"<<endl;
		//if ((sum1 - (*allClustsIter).ysum)>0.001)
		//{
		//	double dif = sum1 - (*allClustsIter).ysum;
		//	double sum2 = (*allClustsIter).ysum;
		//	cout <<"error sum"<<endl;
		//	exit(0);
		//}//end of if
		count ++;
	}
*/
}//end of operation

int singleCompare(const vector <double> ratioVec)
{
	int j;
	int na,nb,nc;
	int number;
	double sum, partialsum;
	double ran, compare;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=30;

	number = ratioVec.size();
	sum = 0;
	for(j=0;j<number;j++)
		sum = sum + ratioVec[j];	
//	ran = genunf(0,1);
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		partialsum = partialsum + ratioVec[j];
		if(partialsum > compare)
			return(j + 1);
	}//end of j
	cout <<"error in decide singelCompare"<<endl;
	exit(5);
}//end of singleCompare

int singleCompareInvert(const vector <double> ratioVec, bool &outsign,const int invert)
{
	int j;
	int na,nb,nc;
	int number;
	double sum, partialsum;
	double ran, compare;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=30;

	number = ratioVec.size();
	sum = 0;
	for(j=0;j<number;j++)
		sum = sum + ratioVec[j];	
//	ran = genunf(0,1);
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		partialsum = partialsum + ratioVec[j];
		if(partialsum > compare)
		{
			if(j%2 ==0)
				outsign = true;
			else
				outsign = false;
			return(j/2 + 1);
			if(invert)
			{
				cout <<"error"<<endl;
				exit(57);
			}
		}//end of if
	}//end of j
	cout <<"error in decide"<<endl;
	exit(5);
}//end of singleCompareInvert

int oldcompare(const vector <double> ratioVec, const int maxShift, int &shift)
{
	int j,k, decide;
	int na,nb,nc;
	int number,shiftSize, clustSize;
	double sum, partialsum, insum, accumu;
	double ran, compare;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=250;
	
	number = ratioVec.size();
	shiftSize = maxShift;//2*maxShift + 1;
	clustSize = number / shiftSize;
	if(number % shiftSize !=0)
	{
		cout <<"ratio vector size problem."<<endl;
		exit(23);
	}
	sum = FACTOR_0;
	for(j=0;j<number;j++)
		sum = sum + ratioVec[j];	
//	ran = genunf(0,1);
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		shift = 0;//maxShift;
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<clustSize;j++)
	{
		insum = 0;
		for(k=0;k<shiftSize;k++)
		{
			insum = insum + ratioVec[j*shiftSize+k];	
		}
		partialsum = partialsum + insum;
		if(partialsum > compare)
		{
			decide = j;
			break;
		}//end of if
	}//end of j
	ran = unran(&na,&nb,&nc);
	compare = insum * ran;
	accumu = 0;
	for(k=0;k<shiftSize;k++)
	{
		accumu = accumu + ratioVec[decide*shiftSize+k];
		if(accumu > compare)
		{
			shift = k;
			return(decide + 1);
		}//end of if
	}//end of k
	cout <<"error in decide"<<endl;
	exit(50);
}//end of oldcompare

int compare(const vector <double> ratioVec, const int maxShift, int &shift)
{
	int j, decide;
	int na,nb,nc;
	int number, clustSize;
	double sum, partialsum;
	double ran, compare;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=250;
	
	number = ratioVec.size();
	clustSize = allClusts.size();//number / maxShift;
	if(number % (maxShift+1) !=0)
	{
		cout <<"ratio vector size problem."<<endl;
		exit(23);
	}
	sum = FACTOR_0;
	for(j=0;j<number;j++)
		sum = sum + ratioVec[j];	
//	ran = genunf(0,1);
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		shift = 0;//maxShift;
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		partialsum = partialsum + ratioVec[j];
		if(partialsum > compare)
		{
			decide = j/maxShift;
			shift = j%maxShift;
			return(decide + 1);
		}//end of if
	}//end of j
	cout <<"error in decide"<<endl;
	exit(50);
}//end of compare

int compareInvert(const vector <double> ratioVec, const int maxShift, int &shift, bool &outsign,const int invert)
{
	int j, decide;
	int na,nb,nc;
	int number, clustSize;
	double sum, partialsum;
	int leftover;
	double ran, compare;

//	seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=250;
	
	number = ratioVec.size();
	clustSize = allClusts.size();//number / maxShift;
	if(number % ((maxShift+1)*2) !=0)
	{
		cout <<"ratio vector size problem."<<endl;
		exit(23);
	}
	sum = FACTOR_0;
	for(j=0;j<number;j++)
		sum = sum + ratioVec[j];	
//	ran = genunf(0,1);
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		shift = 0;//maxShift;
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		partialsum = partialsum + ratioVec[j];
		if(partialsum > compare)
		{
			decide = j/((maxShift + 1)*2);//12/25/05
			leftover = j%((maxShift + 1)*2);//12/25/05
			shift = leftover % (maxShift + 1);
			if(leftover < (maxShift + 1))//> 12/25/05
				outsign = true;
			else
				outsign = false;
			return(decide + 1);
			if(!invert)
			{
				cout <<"error"<<endl;
				exit(55);
			}
		}//end of if
	}//end of j
	cout <<"error in decide"<<endl;
	exit(50);
}//end of compareInvert

double bayesratio(const int geneID, const int shift, const int length, const double *ysum, 
				  const double *ysum2, const int nsize, const double priorbeta0, const double priora, double priorb)// const int column,
{
	int j;
	double logsum, logratio;
	double *ysumnew, *ysum2new;
	double part1,part2,part3, part4, part5,part6;
	list <clusterClass>::iterator mem;
	
	ysumnew = new double[length];
	ysum2new = new double[length];
	/*
	if(self)//cluster contain the gene in operation
	{
		for(j=0;j<length;j++)
		{
			ysumnew[j] = ysum[j] - expre[j];
			ysum2new[j] = ysum2[j] - expre[j]*expre[j];
			if(ysum2new[j]<0)
			{	
				cout <<"error"<<endl;
				exit(2);
			}//end of if
		}//end of j
	}//end of if
	else
	{
		for(j=0;j<length;j++)
		{			
			ysumnew[j] = ysum[j];
			ysum2new[j] = ysum2[j];
			ysum[j] = ysum[j] + expre[j];
			ysum2[j] = ysum2[j] + expre[j]*expre[j];
		}//end of j
		nsize ++;
	}//end of else
	*/
	for(j=0;j<length;j++)
	{			
		if(expre[geneID][j + shift]> -900)// test if it is missing, added 12/11/05
		{
			ysumnew[j] = ysum[j] + expre[geneID][j + shift];
			ysum2new[j] = ysum2[j] + expre[geneID][j + shift] *expre[geneID][j + shift];
		}//end of if
		else
		{
			ysumnew[j] = ysum[j];
			ysum2new[j] = ysum2[j]; 
		}//end of else
	}//end of j
	logsum = 0;
	for(j=0;j<length;j++)
	{
		if(expre[geneID][j + shift]> -900)// test if it is missing, added 12/11/05
		{
			//part1 = -(double)nsize/2*log(2*pie)-0.5*log((double)nsize+1);
			part1 = ((double)nsize/2+priora)*log(priorb + 0.5*(ysum2new[j]+priorbeta0*priorbeta0-(ysumnew[j]+priorbeta0)*(ysumnew[j]+priorbeta0)/((double)nsize + 1)));
			//part3 = -(double)(nsize-1)/2*log(2*pie)-0.5*log((double)nsize);
			double d = ysum2[j]+priorbeta0*priorbeta0-(ysum[j]+priorbeta0)*(ysum[j]+priorbeta0)/((double)nsize);
			part2 = ((double)(nsize-1)/2+priora)*log(priorb + 0.5*(ysum2[j]+priorbeta0*priorbeta0-(ysum[j]+priorbeta0)*(ysum[j]+priorbeta0)/((double)nsize)));
			if (part2<-1000)
			{
				cout <<"part2 wrong"<<part2<<endl;
				exit(0);
			}
			part3 = (priora + 0.5)*log(priorb + 0.5*(expre[geneID][j + shift]-priorbeta0)*(expre[geneID][j + shift]-priorbeta0));//added +shift 11/30
			part4 = 0.5*(log(2.0) + log((double)nsize) - log((double)nsize + 1)); 
			part5 = priora*log(priorb) - gammaln(priora);
			part6 = gammaln(priora+0.5*nsize) -gammaln(priora +0.5) - gammaln(priora + 0.5*(nsize-1));
			logsum = logsum + part2 + part3 - part1 + part4 - part5 + part6; 
			/*if (logsum<-1000)
			{
				mem= member[geneID];
				list < int >::iterator geneIDIter;
				double sum=0;
				for(geneIDIter = (*mem).geneID.begin();geneIDIter != (*mem).geneID.end(); geneIDIter++)
					sum = sum + expre[(*geneIDIter)][shift];
				cout <<sum<<" "<<ysum<<" "<<expre[geneID][shift]<<endl;
				cout <<ysum[j]<<" " <<ysumnew[j]<<" "<<part1<<" "<<part2<<" "<<part3<<" "<<part4<<" "<<part5<<" "<<part6<<endl<<endl;
				cout <<"testing only. logsum <-1000."<<endl;
				exit(3);
			}*/
		}//end of if
	}//end of j
	logratio = logsum;
	delete [] ysumnew;
	delete [] ysum2new;
	return(logratio);
}//end of bayesratio

double bayesratioinvert(const int geneID, const int shift, const int length, const double *ysum, 
				  const double *ysum2, const int nsize, const double priorbeta0, const double priora, double priorb)// const int column,
{
	int j;
	double logsum, logratio;
	double *ysumnew, *ysum2new, *invexpre;
	double part1,part2,part3, part4, part5,part6;
	list <clusterClass>::iterator mem;
	
	ysumnew = new double[length];
	ysum2new = new double[length];
	invexpre = new double[length];
	for(j=0;j<length;j++)
	{			
		if(expre[geneID][j + shift]> -900)// test if it is missing, added 12/11/05
		{
			invexpre[j] = 2*average[geneID] - expre[geneID][j + shift];
			ysumnew[j] = ysum[j] + invexpre[j];
			ysum2new[j] = ysum2[j] + invexpre[j] * invexpre[j];
		}//end of if
		else
		{
			invexpre[j] = 0;
			ysumnew[j] = ysum[j];
			ysum2new[j] = ysum2[j]; 
		}//end of else
	}//end of j
	logsum = 0;
	for(j=0;j<length;j++)
	{
		if(expre[geneID][j + shift]> -900)// test if it is missing, added 12/11/05
		{
			//part1 = -(double)nsize/2*log(2*pie)-0.5*log((double)nsize+1);
			part1 = ((double)nsize/2+priora)*log(priorb + 0.5*(ysum2new[j]+priorbeta0*priorbeta0-(ysumnew[j]+priorbeta0)*(ysumnew[j]+priorbeta0)/((double)nsize + 1)));
			//part3 = -(double)(nsize-1)/2*log(2*pie)-0.5*log((double)nsize);
			//double d = ysum2[j]+priorbeta0*priorbeta0-(ysum[j]+priorbeta0)*(ysum[j]+priorbeta0)/((double)nsize);
			part2 = ((double)(nsize-1)/2+priora)*log(priorb + 0.5*(ysum2[j]+priorbeta0*priorbeta0-(ysum[j]+priorbeta0)*(ysum[j]+priorbeta0)/((double)nsize)));
			part3 = (priora + 0.5)*log(priorb + 0.5*(invexpre[j]-priorbeta0)*(invexpre[j]-priorbeta0));//added +shift 11/30
			part4 = 0.5*(log(2.0) + log((double)nsize) - log((double)nsize + 1)); 
			part5 = priora*log(priorb) - gammaln(priora);
			part6 = gammaln(priora+0.5*nsize) -gammaln(priora +0.5) - gammaln(priora + 0.5*(nsize-1));
			logsum = logsum + part2 + part3 - part1 + part4 - part5 + part6; 
		}//end of if
	}//end of j
	logratio = logsum;
	delete [] ysumnew;
	delete [] ysum2new;
	delete [] invexpre;
	return(logratio);
}//end of bayesratioinvert

void readData(const char dataFileName[], vector <double> &data, int *nCol)
{
	int count=0;
	int nRow;
	double temVal;
	istringstream iss;
    string lineString;
	string name;
//308	int clusterID;//*****308

	ifstream inFile(dataFileName);
    if (! inFile)
    {
        cout << "Error opening input file" <<dataFileName<< endl;
        exit(0);
		//return -1;
    }
// The first line is just a header, so I read in, and then ignore it.
// Note, I could also use the "ignore" function.
//	getline(inFile, lineString);
	count=0;
	while (inFile)
	{
//***** added to allow input operon, gene name and true class label 12/09/05
//308		inFile >> name >> clusterID;//for galactose.data//*****//308
		inFile >> name;// >> clusterID;//for yeast data//*****//412
		geneName.push_back(name);
//308		trueClass.push_back(clusterID);//*****//308
//***** added to allow input operon, gene name and true class label 12/09/05
		if(inFile)
		{
			getline(inFile,lineString);
			iss.clear();
			iss.str(lineString +" ");
			while(iss)
			{
				iss >> temVal;
				if(iss)
				{
					// *****for yeast data, take log,
/*					if(temVal>0)
						temVal = log(temVal);
					else
						temVal = -999;
					// *****for yeast data, take log
*/					data.push_back(temVal);
				}//end of if
			}//end of while
		}//end of if
		count ++;
	}//end of while
	nRow=count-1;
	*nCol=data.size()/nRow;
}//end of readData

void standard(vector <double> &data, const int nCol, const int nRow)
{
	int j,k,len;
	double sum, sum2;
	double mean, stddev;
	int count;
	vector <double> newdata;

	len = nRow *nCol;
	if(len != data.size())
	{
		cout <<"no match"<<endl;
		exit(0);
	}
	for(j=0;j<len;j++)
		newdata.push_back(data[j]);
	for(j=0;j<nCol;j++)
	{
		sum =0;
		sum2 =0;
		count = 0;
		for(k=0;k<nRow;k++)
		{
			if(data[j] < -900)
				count++;
			else
			{	
				sum = sum + data[j + k*nCol];
				sum2 = sum2 + data[j + k*nCol]*data[j + k*nCol];
			}
		}//end of k
		mean = sum/((double) nRow - count);
		stddev = sqrt((sum2 - sum*mean)/((double) nRow - count -1));
		for(k=0;k<nRow;k++)
		{
			if(data[j] < -900)
				newdata[j + k*nCol] = (data[j + k*nCol]-mean)/stddev;
		}
	}//end of j
	for(j=0;j<len;j++)
		data[j] = newdata[j];
}//end of standard

void simudata()
{
	int j,k,m;
	int group, row, column;
	int total;
	double sum, sum2;
	double **mean,**stdev;
	ofstream dataFile;	
	double **a;

	group = 2;
	row = 10;
	column = 5;

	mean = new double*[group];
	stdev = new double*[group];
	for(j=0;j<group;j++)
	{
		mean[j] = new double[column];
		stdev[j] = new double[column]; 
		for(k=0;k<column;k++)
		{
			mean[j][k] = (double) 1 + 5*j + 2*k;
			stdev[j][k] = 1.0;
		}	
	}
	dataFile.open("input.txt");
	if(!dataFile)
	{
		cout << "ERROR: Unable to open file: input.txt." << endl;
	    exit(3);
	}//end of if
	total = row*group;
	a = new double*[total];
	for(j=0;j<total;j++)
	{
		a[j] = new double[column];
	}//end of j
	for(j=0;j<group;j++)
	{
		for(k=0;k<column;k++)
		{
			sum = 0;
			sum2 =0;
			for(m=0;m<row;m++)
			{
				a[j*row + m][k] = gennor(mean[j][k],stdev[j][k]);
				sum = sum + a[j*row+m][k];
				sum2 = sum2 + a[j*row+m][k] * a[j*row+m][k];
			}//end of m
			cout <<sum/row<<"  "<<sqrt((sum2-sum*sum/row)/(row-1))<<endl;
		}//end of k
	}//end of j
	//output data
	for(j=0;j<total;j++)
	{
		for(k=0;k<column;k++)
			dataFile <<a[j][k]<<" ";
		dataFile <<endl;
	}
	for(j=0;j<group;j++)
	{
		delete [] mean[j];
		delete [] stdev[j];
	}
	delete [] mean;
	delete [] stdev;
	for(j=0;j<total;j++)
		delete [] a[j];
	delete [] a;
	dataFile.close();
}//end of simudata

double gammaln(double xx)
{
	double ser,stp,tmp,x,y,cof[6],gam;
	int j;
	cof[0]=76.18009172947146;
	cof[1]=-86.50532032941677;
	cof[2]=24.01409824083091;
	cof[3]=-1.231739572450155;
	cof[4]=0.1208650973866179*0.01;
	cof[5]=-0.5395239384953*0.00001;
	stp=2.5066282746310005;
	x=xx;
	y=x;
	tmp=x+5.5;
	tmp=(x+0.5)*log(tmp)-tmp;
	ser=1.000000000190015;
	for (j=0;j<6;j++) 
	{
		y=y+1.0;
		ser=ser+cof[j]/y;
	}
	gam=tmp+log(stp*ser/x);
	return gam;
}

double unran(int *na, int *nb, int *nc)
{
	double random;
	*na=(171*(*na))%30269;
	*nb=(172*(*nb))%30307;
	*nc=(170*(*nc))%30323;
	random=(double) *na/30269.0+(double) *nb/30307.0+(double) *nc/30323.0;
	random=random-floor(random);
	return random;
}

void cooccur(const int row, int **&occur)
{// determine co-occurance matrix using cluster result
	int j;
	list < clusterClass >::iterator clustPtr;
	list < int >::iterator allGeneIDIter;

// go through each gene, first determine which cluster it is in, 
// then add all its neighbors to the matrix. 
	for(j=0;j<row;j++)
	{
		clustPtr = member[j];
		for(allGeneIDIter = (*clustPtr).geneID.begin();allGeneIDIter != 
			(*clustPtr).geneID.end(); allGeneIDIter++)
		{
			if((*allGeneIDIter)>j)
				occur[j][(*allGeneIDIter)] ++;
			else if((*allGeneIDIter)<j)
				occur[(*allGeneIDIter)][j] ++;
		}//end of allGeneIDIter
	}//end of j
}//end of cooccur
