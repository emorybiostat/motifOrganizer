// clustering program BMC 2.0
// origin from bmc by Steve Qin, 01/15/02
// modified on 03/16/04 using STL list and map
// modified 02/18/06 to work for long list
// 2.3
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <ctype.h>
#define ANNEALING 0
#define INITIAL_CLUSTER 100
//#define IN_WIDTH 16//16//6//2
#define HEAPSORT class sortclustClass
#define KEY postprob
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

using namespace std;

double gammaln(double xx);
static double beta[4]={1.0,1.0,1.0,1.0};
static double betasum=beta[0]+beta[1]+beta[2]+beta[3];
static double sumbeta=gammaln(beta[0])+gammaln(beta[1])+gammaln(beta[2])+gammaln(beta[3]);
static double betasumhalf=0.5*betasum;
static double sumbetahalf=gammaln(0.5*beta[0])+gammaln(0.5*beta[1])+gammaln(0.5*beta[2])+
				gammaln(0.5*beta[3]);
static short kstar[4]={3,2,1,0};//{1,0,3,2};

class motifClass
{
private:
public:
	int id;
	string name;
	bool isEven;
	int width;
	int depth;
	int motifshift;
	int clustshift;
	double postprob;
	vector <string> basepair;
	vector <int> profile[4];
	void printSeq();
	void genProfile();
//	friend void readSeq(list <clustClass > &clustersList);
};

vector <motifClass> allMotifs;

class clustClass
{
private:
public:
	list <int> motifsID;
	vector <bool> power;
	int width;
	int left;
	int right;
	double remain;
	vector <int> profile[4];
	void printClust();
	void clustClass::iniProfile(short bin);
	void clustClass::addMotif(/*vector <motifClass> allMotifs,*/ int motifOrder);
	int clustClass::minusMotif(/*vector <motifClass> allMotifs,*/ int motifOrder);
	void clustClass::updateProfile(short sign, /*vector <motifClass> allMotifs,*/ int motifOrder);
	double clustClass::loglikelihood(const int stage);
};

class sortclustClass
{
private:
public:
	int order;
	list <int> motifsID;
	int size;
	int left;
	int right;
	int width;
	vector <bool> power;
	double remain;
	double postprob;
};

void motifClass::printSeq()
{
	int j;

	for(j=0;j<depth;j++)
		cout << basepair[j] <<endl;
}

void motifClass::genProfile()
{
	int j,k;
	
// initialize
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			profile[j].push_back(0);
	}
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			profile[j][k]=0;
	}
// get count
	for(j=0;j<depth;j++)
	{
		for(k=0;k<width;k++)
		{
			//cout << basepair[j]<<endl;
			switch (basepair[j][k])
			{
			case 'a': profile[0][k]++; break;
			case 'c': profile[1][k]++; break;
			case 'g': profile[2][k]++; break;
			case 't': profile[3][k]++; break;
			case 'A': profile[0][k]++; break;
			case 'C': profile[1][k]++; break;
			case 'G': profile[2][k]++; break;
			case 'T': profile[3][k]++; break;
			default: cout << "error in motif profile" << name <<endl;
			}//end of switch
		}//end of k
	}//end of j
/*	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << profile[j][k] << " ";
		cout << endl;
	}//end of j
	cout << endl;
*/
}//end of motifClass::genProfile

void clustClass::printClust()
{
	int j,k;
	int count=1;//j,k;
	list <int>::iterator motifsIDIter;

	count=1;
	for(motifsIDIter = motifsID.begin(); motifsIDIter != motifsID.end(); motifsIDIter++)
	{
//		outputFile << "ID " << count;
//		outputFile << "Name " << (*motifsIDIter).name << " " <<(*motifsIDIter).width << 
//			" " <<(*motifsIDIter).depth << " " <<(*motifsIDIter).isEven << endl;
		//(*motifsIDIter).printSeq();
//218		cout << "ID " << count;
//218		cout << "Name " << allMotifs[(*motifsIDIter)].name << " " <<allMotifs[(*motifsIDIter)].width << 
//218			" " << allMotifs[(*motifsIDIter)].depth << " " << allMotifs[(*motifsIDIter)].isEven << endl;
		count++;
	}//end of motifsIDIter
/* 218
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
		{
//			outputFile << profile[j][k] <<" ";
			cout << profile[j][k] << " ";
		}
		cout << endl;
//		outputFile << endl;
	}//end of j
218 */
//	outputFile << endl;
//	outputFile.close();
}//end of clustClass::printClust()

//----- initialize the profile matrix when the first motif join this cluster

void clustClass::iniProfile(short bin)
{
	int j,k;
	list <int>::iterator motifsIDIter;
//	int motifWidth;

	left=0;
	right=0;
	
	motifsIDIter=motifsID.begin();
	//int a = (*motifsIDIter).start;
	//int b = (*motifsIDIter).shift; 
//	start=(*motifsIDIter).start + (*motifsIDIter).shift;//added (*motifsIDIter).shift 
	width=allMotifs[(*motifsIDIter)].width;
// ******* moved from before motifsIDIter
	if(bin==2)//added 08/17/05
		remain = width;//added 08/17/05
	else//added 08/17/05
		remain=(double) width*0.5;
// ******* moved from before motifsIDIter
// first set profile matrix to 0
	if(width<0)
	{
		cout <<"cluster width likely not specified yet."<<endl;
		exit(20);
	}
	for(j=0;j<width;j++)
	{
		power.push_back(1);
		for(k=0;k<4;k++)
			profile[k].push_back(0);
	}
	for(j=0;j<width;j++)
	{
		//power[j]=1;
		for(k=0;k<4;k++)
		{
			int a = allMotifs[(*motifsIDIter)].profile[k][j];
			profile[k][j]=allMotifs[(*motifsIDIter)].profile[k][j];//+(*motifsIDIter).profile[kstar[k]][motifWidth-1-j-start];
		}
	}//end of j
}//end of clustClass::iniProfile

void clustClass::addMotif(/*vector <motifClass> allMotifs,*/ int motifOrder)
{
	motifsID.push_back(motifOrder);
	updateProfile(1, motifOrder);
//218	cout <<"shift= "<<allMotifs[0].motifshift<<endl;
}

int clustClass::minusMotif(/*vector <motifClass> allMotifs,*/ int motifOrder)
{
	if(motifsID.size()>1)// there are ore han one element in the cluster
	{
		motifsID.pop_front();
		updateProfile(-1, motifOrder);
		return(0);
	}
	else if (motifsID.size() == 1)//there is only one element in this cluster
		return(1);
	else 
	{
		cout << "error with element number"<<endl;
		return(-1);
	}
}

void clustClass::updateProfile(short sign, /*vector <motifClass> allMotifs,*/ int motifOrder)
{
	int j,k;
	list <int>::iterator motifsIDIter;
	int newwidth;
	int motifshift,clustshift;
	bool begin = true;

	motifshift = allMotifs[motifOrder].motifshift;
	clustshift = allMotifs[motifOrder].clustshift;
/* 02/18/06
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << profile[j][k] << " ";//removed +clustshift
		cout <<endl;
	}
	cout <<endl;
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << allMotifs[motifOrder].profile[j][k+motifshift] << " ";
		cout <<endl;
	}
	cout <<endl;
02/18/06 */
// update cluster width as needed
	if((sign == 1)&&(allMotifs[motifOrder].width < width))//add motif, cluster width may decrease
	{
		newwidth = allMotifs[motifOrder].width;
	}//end of if
	else if((sign == -1)&&(motifsID.size()==1))//remove motif, cluster width may increase 
	{
		motifsIDIter = motifsID.begin();
		newwidth = allMotifs[(*motifsIDIter)].width;
		allMotifs[(*motifsIDIter)].motifshift = 0;
		allMotifs[(*motifsIDIter)].clustshift = 0;
	}//end of else if
	else
		newwidth = width;
	if((sign == -1)&&(motifsID.size()==1))
	{
		motifsIDIter = motifsID.begin();
		for(j=0;j<4;j++)
		{
			profile[j].clear();
			for(k=0;k<newwidth;k++)
			{	
				profile[j].push_back(allMotifs[(*motifsIDIter)].profile[j][k]);
//218				cout << profile[j][k] << " ";	
			}//end of k
//218			cout <<endl;
		}//end of j

	}//end of if
	else if(width == newwidth)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<width;k++)
			{
				//int e= allMotifs[motifOrder].profile[j][k];
				//int f = profile[j][k];
				profile[j][k]=profile[j][k]+sign * allMotifs[motifOrder].profile[j][k+motifshift];
				if(profile[j][k] <0)
				{
					int e = profile[j][k];
					int f = allMotifs[motifOrder].profile[j][k+motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
 						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}
			}//end of k
		}//end of j
	}//end of if
	else if(width < newwidth)//need to add profile, increase only when there is just 1 motif in the cluster
	{
		for(j=0;j<4;j++)
		{
			for(k=width;k<newwidth;k++)
				profile[j].push_back(0);
		}//end of j
		for(j=0;j<4;j++)
		{
			for(k=0;k<width;k++)
			{
				int g = profile[j][k];
				int h = allMotifs[motifOrder].profile[j][k+motifshift];
				profile[j][k]=profile[j][k] - allMotifs[motifOrder].profile[j][k+motifshift];
				if(profile[j][k] <0)
				{
					int e = profile[j][k];
					int f = allMotifs[motifOrder].profile[j][k+motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}//end of if
			}//end of k
			motifsIDIter = motifsID.begin();
			for(k=width;k<newwidth;k++)
			{
				int g = allMotifs[(*motifsIDIter)].profile[j][k];
				profile[j][k] = allMotifs[(*motifsIDIter)].profile[j][k];
			}//end of k
		}//end of j
	}//end of else if
	else//need to subtract profile, since a shorter motif was added to the cluster
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				profile[j][k] = profile[j][k + clustshift]+allMotifs[motifOrder].profile[j][k];
			// need to update motif shift
		}//end of j
		for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
		{
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
			int e = (*motifsIDIter);
			int f = allMotifs[(*motifsIDIter)].motifshift;
			if((*motifsIDIter) != motifOrder)
				allMotifs[(*motifsIDIter)].motifshift = allMotifs[(*motifsIDIter)].motifshift + clustshift;
//				allMotifs[(*motifsIDIter).id].motifshift = allMotifs[(*motifsIDIter).id].motifshift + clustshift; 
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
		}//end of motifsIDIter
	}//end of else
	width = newwidth;
}//end of clustClass::updateProfile

double clustClass::loglikelihood(const int stage)
{
	int j,k;
	int motifsSize;
	int la, lsum = 0, sumin,start;
	double logsum,tem;
	list <int>::iterator motifsIDIter;

	tem=gammaln(betasum)-sumbeta;
	motifsSize = motifsID.size();
	lsum=0;
	logsum =0;
 	for(j=0;j<4;j++)
		lsum = lsum + profile[j][0];
	for(j=0;j<width;j++)
	{
		if((stage!=1)||(power[j]))
		{
			for(k=0;k<4;k++)
			{
				int a = profile[k][j];
				la = (double) profile[k][j];
				logsum=logsum+gammaln(la+beta[k]);
			}//end of k
			logsum=logsum-gammaln((double) lsum + betasum);
		}//end of if
		else if ((stage==1)&&(!power[j]))
		{
			for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
			{
				sumin=allMotifs[(*motifsIDIter)].depth;
				start=allMotifs[(*motifsIDIter)].motifshift;// add (*motifsIDIter).shift 
				for(k=0;k<4;k++)
				{
					la = (double) allMotifs[(*motifsIDIter)].profile[k][j+start];//+beta[k];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) sumin + betasum);
			}//end of 
			logsum=logsum+(double) (motifsSize-1)*tem;
		}//end of else if
	}//end of j
	logsum = logsum + (double) width *tem;
	return(logsum);
}//end of clustClass::loglikelihood

list <clustClass> allClusts;
vector <list <clustClass>::iterator> motifHome;

int readSeq(vector <motifClass> &fullMotifVector, const char outPutFileName[]);

void initialCluster(const short bin, const int motifCount, const vector <motifClass> fullMotifVector);

double bmc(const short bin, const int numChains, const int numItera, const int W, const int WMIN, 
		 const int WMAX, const int MAX_SHIFT, const double factor_0, list <clustClass> &finalClusters);

double loglikelihoodSingle(const short stage);

double chain(const int NUM_CYCLE, const int W, const int WMIN, const int WMAX, const int MAX_SHIFT, 
			 const short bin, const double FACTOR_0);

double operation(const int motifOrder, const short bin, const double FACTOR_0,
				 int &del, const short stage);

double clustInsertMotif(const int motifOrder, const int decision, const short stage);

double clustRemoveMotif(const int motifOrder, const list <clustClass>::iterator clustPoint, const short stage);

double loglikelihoodUpdateSingle(const int motifOrder, const list <clustClass>::iterator allClustsIter, const short stage);

//double diffWidthLoglikelihoodUpdateSingle (const int oldwidth, const int newwidth, const int motifOrder, 
//										   const list <clustClass>::iterator allClustsIter,const short stage);

double onoffMetro(const list <clustClass>::iterator allClustsIter, int W, int WMIN, short bin);

double metEven(const list <clustClass>::iterator allClustsIter, short indel, short position, short W);

double metOdd(const list <clustClass>::iterator allClustsIter, short indel, short position, short W);

double metSingle(const list <clustClass>::iterator allClustsIter, short indel, short position, short W);

double fragment(const short bin, const int WMIN, const int WMAX, const int W);

double bayesRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart, 
						const short stage);

int decide(vector <double> &ratioVec, const double FACTOR_0, double &postprob);

int shiftdecide(vector <vector <double> > &ratioAll, const double FACTOR_0,int &shift, double &postprob);

double unran(int *na, int *nb, int *nc);

double emptyLoglikelihoodupdateSingle(const int motifOrder, const list <clustClass>::iterator allClustsIter, 
			   const short stage);

void initWidth(const short bin, const int W);

void resultDisplay(const int motifCount, const double loglike, const short bin, 
				   const double factor_0, list <clustClass> clusters, const char outPutFileName[]);

void mqheap(HEAPSORT *p[], short total, short begin, short end);

static void sift(HEAPSORT *p[], short k, short j, short m);

double finalRatioEven(list <clustClass>::iterator allClustsIter);

double finalRatioOdd(list <clustClass>::iterator allClustsIter);

double finalRatioSingle(list <clustClass>::iterator allClustsIter);

int checkWidth(	list <clustClass>::iterator clustPoint);

int main()
{
	list <clustClass> finalClusters;
	vector <motifClass> fullMotifVector;
	int W,WMIN,WMAX;
	int  MAX_SHIFT;
	int numChains,numIter;
	double factor_0;
	int motifCount;
	short bin=0;
	double loglike;
	int IN_WIDTH;
	char inputFileName[]="tra.txt";

	IN_WIDTH = 6;//5;//20;
	W=IN_WIDTH;//6;
	WMIN=6;//5;//16;//16;//4;//less than WMAX if consider fragmentation
	WMAX=W;//16;//6;
	MAX_SHIFT=1;//5;
	numChains=1;
	numIter=50;
	factor_0=1.0;
	bin=2;//0;

//	srand((unsigned)time(NULL));
// read in data, initial assign clusters and calculate cluster profiles
	motifCount=readSeq(fullMotifVector,inputFileName);
	initialCluster(bin,motifCount, fullMotifVector);
	//exit(0);
// clustering
//-	cout << "cluster size= " << allClusts.size();
	loglike = bmc(bin,numChains,numIter, W, WMIN,WMAX,MAX_SHIFT,factor_0,finalClusters);
// output result
	resultDisplay(motifCount, loglike, bin, factor_0, finalClusters,"new01.txt");
	return(0);
}//end of main

void resultDisplay(const int motifCount, const double loglike, const short bin, 
				   const double factor_0, list <clustClass> clusters, const char outPutFileName[])
{
	int j,k,m;
	list <clustClass>::iterator clustIter;
	list <int>::iterator motifsIDIter;
	ofstream outPutFile;
	ofstream motifblockFile;// FASTA format motif sequence alignment file
	int finalClustNumber;
	int tot,count;
	double value =0;
	int motifsize;
	sortclustClass *clustarray, **sorted;
	int start;
	int minWidth, maxWidth;
	bool begin;

	outPutFile.open(outPutFileName);
	if(!outPutFile)
	{
		cout << "ERROR: Unable to open file: " << outPutFileName << endl;
	    exit(3);
	}//end of if	
	motifblockFile.open("motifblock.txt");
	if(!motifblockFile)
	{
		cout << "ERROR: Unable to open file: " << "motifblock.txt"<<endl;
	    exit(3);
	}//end of if
	finalClustNumber = clusters.size();
// find out minimum and maximum cluster width
	begin = true;
	for(clustIter = clusters.begin();clustIter != clusters.end();clustIter ++)
	{
		if(begin)
		{
			minWidth = (*clustIter).width;
			maxWidth = minWidth;
			begin = false;
		}
		else 
		{
			if ((*clustIter).width > maxWidth)
				maxWidth = (*clustIter).width;
			if ((*clustIter).width < minWidth)
				minWidth = (*clustIter).width;
		}//end of else
	}
	outPutFile << "**************************************** \n";
	outPutFile << "*                                      * \n";
	outPutFile << "*       BMC 2.1 clustering Result      * \n";
	outPutFile << "*                                      * \n";
	outPutFile << "**************************************** \n\n";
	
	outPutFile << " Total number of all motifs: " << motifCount << "\n";
	    cout << "\n Total number of all motifs: " << motifCount << "\n";
	outPutFile << "\nThere are total of " << finalClustNumber << " clusters\n";
    cout <<"\n There are total of of " << finalClustNumber << " clusters, see output file " << outPutFileName << " for details\n\n";
	outPutFile << "q=" << factor_0 << ", motif width range = " <<minWidth <<"--"<<maxWidth<<endl;
	outPutFile << "log likelihood = " << loglike <<endl<<endl;
    cout <<"output log likelihood = "<<loglike <<endl<<endl;//-

//    finalprob(finalClustNumber,beta,W,MAX_SHIFT,bin);        
//sort	
	count = 0;
	clustarray = new sortclustClass[finalClustNumber];
	sorted = new sortclustClass *[finalClustNumber];
	for(clustIter = clusters.begin();clustIter != clusters.end();clustIter ++)
	{
		clustarray[count].order = count;
		clustarray[count].left = (*clustIter).left;
		clustarray[count].right = (*clustIter).right;
		clustarray[count].remain = (*clustIter).remain;
		clustarray[count].motifsID = (*clustIter).motifsID;
		clustarray[count].power = (*clustIter).power;
		clustarray[count].width = (*clustIter).width;
		/*if(bin==0)
			value=finalRatioEven(clustIter);
		else if(bin==1)
			value=finalRatioOdd(clustIter);
		else
		*/	value=finalRatioSingle(clustIter);
		motifsize = (*clustIter).motifsID.size();
		clustarray[count].size = motifsize;
		clustarray[count].postprob = value /motifsize;
		sorted[count] = &clustarray[count];
		count++;
	}//end of clustIter
/*	for(j=0;j<finalClustNumber;j++)
    {
        clust[j].order=j;
        p[j]=&clust[j];
    }//end of j
*/
    mqheap(sorted,finalClustNumber,0,finalClustNumber);
    for (j=finalClustNumber;j>0;j--)
    {
		outPutFile << "cluster " << finalClustNumber-j+1 << ", size= "<< (*sorted[j-1]).size <<"  ";// <<"Bayes Ratio=";
		outPutFile << "width = "<<(*sorted[j-1]).width;//02/18/06<<", number of On columns = " << (*sorted[j-1]).remain << "   ";
//***** delete 02/18/06
/*		if (bin ==0)
		{
			for(k=0;k<(*sorted[j-1]).width/2;k++)
				outPutFile << (*sorted[j-1]).power[k];
		}//end of if
		else if (bin ==1)
		{
			for(k=0;k<(*sorted[j-1]).width/2+1;k++)
				outPutFile << (*sorted[j-1]).power[k];
		}//end of if
		else
		{
			for(k=0;k<(*sorted[j-1]).width;k++)
				outPutFile << (*sorted[j-1]).power[k];
		}//end of if
*/
//***** delete 02/18/06
		outPutFile << " Average Bayes ratio = " << (*sorted[j-1]).postprob <<endl;
		tot=0;
		for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
		{
			tot++;
			//if((tot%2==1)&&(tot!=1)) 
			//	outPutFile << endl;                        ";
            outPutFile << tot << " " << allMotifs[(*motifsIDIter)].name << " [" << allMotifs[(*motifsIDIter)].motifshift 
						<< "]  ("<< allMotifs[(*motifsIDIter)].postprob <<")	";
		}//end of motifsIDIter
		outPutFile <<endl<<endl;
		for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
		{
			start = allMotifs[(*motifsIDIter)].motifshift;
			for(k=0;k<allMotifs[(*motifsIDIter)].depth;k++)
			{
				motifblockFile << ">" << allMotifs[(*motifsIDIter)].name<<"  cluster "<<finalClustNumber-j+1<<endl;
				for(m=0;m<(*sorted[j-1]).width;m++)
					motifblockFile <<allMotifs[(*motifsIDIter)].basepair[k][m+start];
				motifblockFile <<endl;
			}//end of k
			//motifblockFile <<endl;
		}//end of motifsIDIter
		motifblockFile <<endl;
	}//end of j
	delete [] clustarray;
	delete [] sorted;
	outPutFile.close();
	motifblockFile.close();
	/*
	count =0;
	for(clustIter = clusters.begin();clustIter != clusters.end();clustIter ++)
	{
		motifsize = (*clustIter).motifs.size();
		outPutFile << "cluster " << count << ", size= "<< motifsize <<"  ";// <<"Bayes Ratio=";
		outPutFile << "number of On columns = " << (*clustIter).remain << "   ";
		if (bin ==0)
		{
			for(k=0;k<W/2;k++)
				outPutFile << (*clustIter).power[k];
			value=finalRatioEven(W,clustIter);
			outPutFile << " Average Bayes ratio = " <<value/motifsize <<endl;
		}//end of if
		count++;
		tot=0;
		for(motifsIDIter=(*clustIter).motifs.begin();motifsIDIter !=(*clustIter).motifs.end();motifsIDIter++)
		{
			tot++;
		//	if((tot%2==1)&&(tot!=1)) 
		//		outPutFile << endl;                        ";
            outPutFile << tot << " " << (*motifsIDIter).name << " [" << (*motifsIDIter).shift << "]  ("<< (*motifsIDIter).postprob <<")	";
		}//end of motifsIDIter
		outPutFile <<endl<<endl;
	}//end of clustIter
	*/
}//end of resultDisplay

void mqheap(HEAPSORT *p[], short total, short begin, short end)
{
	short j,k;
	HEAPSORT *w;
	//void sift();
	if(begin<0) begin=0;
	if(end>total-1) end=total-1;
	k=end-begin+1;
	for(j=k/2-1;j>=0;j--)
		sift(p,begin,j,k);
	for(j=k-1;j>=1;j--)
	{
		w=p[begin];
		p[begin]=p[j+begin];
		p[j+begin]=w;
		sift(p,begin,0,j);
	}//end of j
	return;
}//end of mqheap

static void sift(HEAPSORT *p[], short k, short j, short m)
{
	short a,b;
	double t; //short
	HEAPSORT *w;
	t=(*p[j+k]).KEY;
	w=p[j+k];
	b=2*(j+1)-1;
	a=j;
	while(b<=m-1)
	{
		if((b<m-1)&&((*p[b+k]).KEY<(*p[b+k+1]).KEY))
			b=b+1;
		if(t<(*p[b+k]).KEY)
		{
			p[a+k]=p[b+k];
			a=b;
			b=2*(a+1)-1;
		}//end of if
		else b=m;
	}//end of while
	p[a+k]=w;
	return;
}//end of shift

double finalRatioEven(list <clustClass>::iterator allClustsIter)
{
    int j,k, start;
	list <int>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, sumex;
	int suminclust,la;
	double difsum;
	double tem;
	int width;

	width = (*allClustsIter).width;
	tem=gammaln(betasum)-sumbeta;
	suminclust=0;
	for(j=0;j<4;j++)
		suminclust=suminclust+(*allClustsIter).profile[j][0];
	suminclust=2*suminclust;
	sumex=0;
	difsum=0;
	if((*allClustsIter).motifsID.size()==1)
	{
	  return 0;
	}
	else
	{
		for(j=0;j<width/2;j++)
		{
			if((*allClustsIter).power[j])//==1
			{
				sumex=0;
				for(k=0;k<4;k++)
				{	
					sumex=sumex+gammaln((*allClustsIter).profile[k][j]
						+(*allClustsIter).profile[kstar[k]][width-1-j]+beta[k]);
				}
				oldProfile=sumex-gammaln((double) suminclust+betasum);
				newProfile=0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = 2*allMotifs[(*motifsIDIter)].depth;
					start = allMotifs[(*motifsIDIter)].motifshift;
					sumex=0;
					for(k=0;k<4;k++)
					{
						la=allMotifs[(*motifsIDIter)].profile[k][j+start]
							+allMotifs[(*motifsIDIter)].profile[kstar[k]][width-j-1+start];
						sumex=sumex+gammaln((double) la+beta[k]);
					}//end of m
					newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
			    }//end of k
				difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
			}//end of if
		}//end of j
		difsum=difsum-(double)((*allClustsIter).remain)*(double)((*allClustsIter).motifsID.size()-1)*tem;
		return(difsum);
	}//end of else
}//end of finalRatioEven

double finalRatioOdd(list <clustClass>::iterator allClustsIter)
{
    int j,k, start;
	list <int>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, suminhalf;
	double sumex,sumexhalf;
	int suminclust,suminclusthalf;
	int la,lahalf;
	double difsum;
	double temhalf,tem;
	int motifsSize;
	int width;

	width = (*allClustsIter).width;

	motifsSize = (*allClustsIter).motifsID.size();
	temhalf = gammaln(betasumhalf)-sumbetahalf;
	tem=gammaln(betasum)-sumbeta;
	suminclusthalf=0;
	for(j=0;j<4;j++)
		suminclusthalf = suminclusthalf+(*allClustsIter).profile[j][0];
	suminclust=2*suminclusthalf;
	difsum=0;
	if(motifsSize==1)
		return 0;
	else
	{
		for(j=0;j<width/2;j++)
		{
			if((*allClustsIter).power[j])//==1
			{
				sumex=0;
				for(k=0;k<4;k++)
				{	
					sumex=sumex+gammaln((*allClustsIter).profile[k][j]
						+(*allClustsIter).profile[kstar[k]][width-1-j]+beta[k]);
				}
				oldProfile=sumex-gammaln((double) suminclust+betasum);
				newProfile=0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = 2*allMotifs[(*motifsIDIter)].depth;
					start = allMotifs[(*motifsIDIter)].motifshift;
					sumex=0;
					for(k=0;k<4;k++)
					{
						la=allMotifs[(*motifsIDIter)].profile[k][j+start]
							+allMotifs[(*motifsIDIter)].profile[kstar[k]][width-j-1+start];
						sumex=sumex+gammaln((double) la+beta[k]);
					}//end of m
					newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
			    }//end of k
				difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
			}//end of if
		}//end of j
//middle column		
		if((*allClustsIter).power[width/2])
		{
			sumexhalf=0;
			for(k=0;k<4;k++)
				sumexhalf = sumexhalf + gammaln((*allClustsIter).profile[k][width/2] + 0.5*beta[k]);
			oldProfile=sumexhalf - gammaln((double) suminclusthalf + betasumhalf);
			newProfile=0;
			for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
			{
				suminhalf = allMotifs[(*motifsIDIter)].depth;
				start = allMotifs[(*motifsIDIter)].motifshift;
				sumexhalf=0;
				for(k=0;k<4;k++)
				{
					lahalf = allMotifs[(*motifsIDIter)].profile[k][width/2+start];
					sumexhalf=sumexhalf+gammaln((double) lahalf + 0.5*beta[k]);
				}//end of m
				newProfile=newProfile+sumexhalf-gammaln((double) suminhalf+betasumhalf);
		    }//end of k
			difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
		}//end of if
		difsum=difsum-(double)((*allClustsIter).remain)*(double)(motifsSize-1)*tem;
		if((*allClustsIter).power[width/2])
			difsum=difsum-(double)(motifsSize-1)*temhalf;
		return(difsum);
	}//end of else
}//end of finalRatioOdd

double finalRatioSingle(list <clustClass>::iterator allClustsIter)
{
    int j,k, start;
	list <int>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, sumex;
	int suminclust,la;
	double difsum;
	double tem;
	int motifsSize;
	int width;

	width = (*allClustsIter).width;
	tem=gammaln(betasum)-sumbeta;
	suminclust=0;
	for(j=0;j<4;j++)
		suminclust=suminclust+(*allClustsIter).profile[j][0];
	sumex=0;
	difsum=0;
	motifsSize = (*allClustsIter).motifsID.size();
	if(motifsSize==1)
	  return 0;
	else
	{
		for(j=0;j<width;j++)
		{
			if((*allClustsIter).power[j])//==1
			{
				sumex=0;
				for(k=0;k<4;k++)
				{	
					sumex=sumex+gammaln((*allClustsIter).profile[k][j]+beta[k]);
				}
				oldProfile=sumex-gammaln((double) suminclust+betasum);
				newProfile=0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
					start = allMotifs[(*motifsIDIter)].motifshift;
					sumex=0;
					for(k=0;k<4;k++)
					{
						la = allMotifs[(*motifsIDIter)].profile[k][j+start];
						sumex=sumex+gammaln((double) la+beta[k]);
					}//end of m
					newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
			    }//end of k
				difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
			}//end of if
		}//end of j
		difsum=difsum-(double)((*allClustsIter).remain)*(double)(motifsSize-1)*tem;
		return(difsum);
	}//end of else
}//end of finalRatioSingle

double bmc(const short bin, const int numChains, const int numItera, const int W, const int WMIN, 
		 const int WMAX, const int MAX_SHIFT, const double factor_0, list <clustClass> &finalClusters)
{
	int j;
	double like,large;

	cout << "start iteration"<<endl;
	for(j=0;j<numChains;j++)
	{
		cout << "chain" << j+1 << endl;
		like=chain(numItera, W, WMIN,WMAX, MAX_SHIFT, bin, factor_0);
		cout << "log likelihood =" << like <<endl;
		if ((j==0) || (like > large))
		{
			large=like;
			finalClusters=allClusts;
		}//end of if
	}//end of j
	return(large);
}//end of bmc

double chain(const int NUM_CYCLE, const int W, const int WMIN, const int WMAX, const int MAX_SHIFT, 
			 const short bin, const double FACTOR_0)
{
	int j;
	int cycle;
	short stage = 0;
	double like,like1,difference;
	list <clustClass>::iterator allClustsIter;
	list <motifClass>::iterator motifsIDIter;
	int currentClusterOrder;
	int del;
	vector <int> delVector;
	double diff =0;

//*****
//  initialize
//*****

//*****
//  get initial likelihood
//*****
	like=loglikelihoodSingle(0);
//	double sum =0;
//	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//		sum = sum + (*allClustsIter).loglikelihood(0);
	//-cout << "original like=" << like <<endl;
//*****
//  start iteration
//*****
	stage = 0;///whether in annealing step
	for(cycle = 0; cycle < NUM_CYCLE+ANNEALING;cycle ++ )
	{
		like=loglikelihoodSingle(0);
		double sum =0;
		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
			sum = sum + (*allClustsIter).loglikelihood(0);

		cout << "iteration=" << cycle+1 <<"\r";//"cluster=" << allClusts.size() << endl;
//***** special check on clusters
 		//-for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		//-	(*allClustsIter).printClust(IN_WIDTH);		
		//***** removed 08/17/05
		//-if(cycle == (NUM_CYCLE/3*2))
		//-{
		//-	stage=1;
		//-	initWidth(bin,W);
		//-}
		//***** moved to the end of the function
		currentClusterOrder=0;
		//-cout <<"cluster total = " << allClusts.size() <<endl;
		for(j=0;j<allMotifs.size();j++)// go through every motif
		{
			cout <<"motif= "<<j<<" cluster size = "<<allClusts.size()<<endl;
//218			cout << "motif= " << j <<" " << "name=" <<allMotifs[j].name /*<< "clust=" << allMotifs[j].cluster*/ <<endl;
//			cout <<"shift= "<<allMotifs[23].motifshift<<endl;
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//			{
//				cout <<"width= "<<(*allClustsIter).width<<endl;
//				cout <<"remain= "<<(*allClustsIter).remain<<endl;
//218			}
			difference=operation(j, bin, FACTOR_0,del,stage);
			like1=loglikelihoodSingle(stage);
			like=like+difference;
//218			cout << "like=" <<like << "like1=" <<like1<<endl;
			if (fabs(like-like1)>0.00001)
			{
				cout << "problem here: motif = "<< j<<" dif=" << difference << "new like=" << like << " like1=" <<like1<<endl;
				exit(18);
			}
//***** special check on clusters
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//				(*allClustsIter).printClust(16);
//***** end special check
		}
//
// fragmentation step
//			
		//-like1=loglikelihoodEven(clusters,W,stage);
//		cout << "like=" <<like << "like1=" <<like1;
		if(cycle == (NUM_CYCLE/3*2))//moved 08/17/05
		{
			stage = 1;//added 08/17/05
			initWidth(bin,W);//added 08/17/05
		}
		if (stage ==1)
		{
		//	diff=fragment(bin,WMIN, WMAX,W);
			//printf("diff=%10.5f\n",diff);
			like=like+diff;
		}//end of if
		//-like1=loglikelihoodEven(clusters,W, stage);
		//-cout << "like=" <<like << "like1=" <<like1;
	}//end of cycle
	//like1=loglikelihoodEven(clusters,WMIN);
	//cout << "like=" << like+difference << " like1=" <<like1<<endl;
	return(like+difference);	
}//end of chain

double operation(const int motifOrder, const short bin, const double FACTOR_0,  
				 int &del, const short stage)
{
	int j;
	int clustOrder = -1;
	vector <double> ratioVec;
	vector < vector <double> > ratioAll;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustPoint;
	int clusterNumber = 0;
	double ratio;
	int decision = -1, currentCluster = -1,shift = 0;
	//int start;
	//motifClass temMotif;
	//list <motifClass> tempMotifs;
	clustClass tempClust;
	double likeminus,likeadd;
	bool empty = false;
	double postprob;
	int motifStart, clustStart;
	int motifWidth, clustWidth;
	double max;
	int newwidth = 0, oldiwdth = 0;
	vector <bool> clustlong;
	int oldmotifshift, oldclustshift;
	double likeold, likenew;

//	the pointer pointed to the cluster that the motif belongs to
	clustPoint = motifHome[motifOrder];
// ********** temp add 02/18/06	
/*	int k;
	for(j=0;j<4;j++)
	{
		int wid = allMotifs[motifOrder].width;
		int shf = allMotifs[motifOrder].motifshift;
		for(k=0;k<allMotifs[motifOrder].width;k++)
			cout << allMotifs[motifOrder].profile[j][k+allMotifs[motifOrder].motifshift] << " ";
		cout <<endl;
	}
	cout <<endl;
*/
//*********** temp add 02/18/06
//  translate clustPoint to clustOrder 
	currentCluster = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if (allClustsIter == clustPoint)
			break;
		else
			currentCluster ++;
	}//end of allClustsIter
	del=-1;
//218	cout << "start of operation, motif = " <<motifOrder<<"size = " <<(*clustPoint).motifsID.size() <<endl;
	// remove motif from cluster 
	int a = (*clustPoint).motifsID.size();
	if((*clustPoint).motifsID.size() == 1) //for cluster with only 1 member, remove one exisitng cluster
	{
		empty = true;
		//likeminus=loglikelihoodUpdateSingle(motifOrder, clustPoint,stage);	
		likeminus = - (*clustPoint).loglikelihood(0);
		allClusts.erase(clustPoint);
//218		cout <<"empty likeminus = "<<likeminus<<endl;
	}//end of if
	else
		likeminus=clustRemoveMotif(motifOrder, clustPoint, stage);
//start trying to fit all clusters
	int tmp = allClusts.size();
	int count =0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		count ++;
// decide shift
		clustWidth = (*allClustsIter).width;
		motifWidth = allMotifs[motifOrder].width;

		if(clustWidth >=motifWidth)//cluster is longer
		{
			clustlong.push_back(true);
			shift = clustWidth - motifWidth +1;
			max = 0;
			for(j=0;j<shift;j++)
			{
				motifStart = 0;
				clustStart = j;
				ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart,stage));
				//if((j==0)||(ratio >max))//removed 09/24/05
				//	max = ratio;			//removed 09/24/05
				//ratioVec.push_back(max); //removed 09/24/05
				ratioVec.push_back(ratio);
			}//end of j
		}//end of if
		else//motif is longer
		{
			clustlong.push_back(false);
			shift = motifWidth - clustWidth +1;
			for(j=0;j<shift;j++)
			{
				motifStart = j;
				clustStart = 0;
				ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart,stage));
				ratioVec.push_back(ratio);
			}//end of j
		}//end of else
		ratioAll.push_back(ratioVec);
		ratioVec.clear();
	}//end of allClustiter
// decide which cluster fit the motif best
	decision=shiftdecide(ratioAll, FACTOR_0,shift,postprob);
	allMotifs[motifOrder].postprob = postprob;
//	cout << "decision=" << decision <<endl;
// if the current cluster is the best, then do nothing
	if((!empty)&&(decision==(currentCluster + 1)))//added !empty 08/07/05
	{// same cluster
		//(*clustPoint).motifs.pop_front();
		oldmotifshift = allMotifs[motifOrder].motifshift;
		oldclustshift = allMotifs[motifOrder].clustshift;
		if(clustlong[currentCluster])
		{
			allMotifs[motifOrder].motifshift = 0;
			allMotifs[motifOrder].clustshift = shift;
		}
		else
		{
			allMotifs[motifOrder].motifshift = shift;
			allMotifs[motifOrder].clustshift = 0;
		}
		//(*clustPoint).updateProfile(1,allMotifs, motifOrder);
		//(*clustPoint).motifs.push_back(allMotifs[motifOrder]);
		clustlong.clear();
		if((oldmotifshift == allMotifs[motifOrder].motifshift)&&(oldclustshift == allMotifs[motifOrder].clustshift))
		{
			(*clustPoint).addMotif( motifOrder);
			return(0);
		}
		else 
		{
			likeold = (*clustPoint).loglikelihood(0);
			(*clustPoint).addMotif( motifOrder);
			likenew = (*clustPoint).loglikelihood(0);
            likeadd = likeold - likenew;
//218       cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
            return(likeminus - likeadd) ;
		}
	}//end of if
// if a different cluster is desired
	else 
	{
// remove it from original cluster
		//(*clustPoint).motifs.pop_front();
		if(decision > 0)
		{
// insert it into the desired cluster
			//likeadd=clustInsertMotif(motifOrder,decision-1,clustWidth, stage);
			//if(MAX_SHIFT>1)
			//allMotifs[motifOrder].motifshift = shift;//update shift, added 08/03/05 //remove -MAX_SHIFT/2 08/06/05
			if(clustlong[decision-1])
			{
				allMotifs[motifOrder].motifshift = 0;
				allMotifs[motifOrder].clustshift = shift;
			}
			else
			{
				allMotifs[motifOrder].motifshift = shift;
				allMotifs[motifOrder].clustshift = 0;
			}
			likeadd=clustInsertMotif(motifOrder,decision-1, stage);// moved down from before if 08/08/05
			// check if need to update cluster width 09/13/05
			//-if(allMotifs[motifOrder].width <oldiwdth)
			//-	newiwdth = allMotifs[motifOrder].width;
		//-	cout << " likeadd1=" <<likeadd<<endl;
		}//end of if
// if a new cluster is desired
		else if (decision == 0)
		{
// add a new cluster to list cluster, with this motif as its only member.
			allMotifs[motifOrder].motifshift = 0;
			allMotifs[motifOrder].clustshift = 0;
//			tempMotifs.push_back(allMotifs[motifOrder]);//(*motifsIDIter);
			tempClust.motifsID.push_back(motifOrder);
			//tempClust.width = allMotifs[motifOrder].width;
			tempClust.iniProfile(bin);
			allClusts.push_back(tempClust);
			clustPoint = allClusts.end();
			clustPoint --;
			motifHome[motifOrder] =clustPoint;
//			tempMotifs.clear();
			//likeadd=loglikelihoodNewEven(clustWidth,motifsIDIter,stage);//***)
			if (empty)
				likeadd=likeminus;
			else
			{
				likeadd=emptyLoglikelihoodupdateSingle(motifOrder,clustPoint,stage);
			}//end of else
		//-	cout << " likeadd2=" <<likeadd<<endl;
		}//end of else if		
// remove it from its current cluster
		del=motifOrder;//.first
		clustlong.clear();
//218		cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
		return(likeminus-likeadd);
	}//end of else
}//end of operation

int checkWidth(	list <clustClass>::iterator clustPoint)
{
	list <int>::iterator allmotifsIDIter;
	int newiwdth;
	bool begin = true;

	begin = true;
	for(allmotifsIDIter = (*clustPoint).motifsID.begin();allmotifsIDIter != (*clustPoint).motifsID.end(); allmotifsIDIter++)
	{
		if(begin)
		{
			newiwdth = allMotifs[(*allmotifsIDIter)].width;
			begin = false;
		}
		else if(allMotifs[(*allmotifsIDIter)].width < newiwdth)
			newiwdth = allMotifs[(*allmotifsIDIter)].width;	
	}//end of allmotifsIDIter
	return(newiwdth);
}

void initWidth(const short bin, const int W)
{
	int k;
	list <clustClass>::iterator clustIter;

	for(clustIter = (allClusts).begin();clustIter != (allClusts).end();clustIter ++)
	{
		(*clustIter).left=0;
		(*clustIter).right=0;
		if(bin==0)//even 
		{
			(*clustIter).remain=(double) W/2;
			for(k=0;k<W/2;k++)
			{
				(*clustIter).power[k]= true;//1;
			}//end of k
		}
		else if(bin==1)//odd
		{
			(*clustIter).remain=(double) W*0.5;
			for(k=0;k<W/2+1;k++)
			{
				(*clustIter).power[k]=true;//1;
			}//end of k
		}
		else// single column
		{
			(*clustIter).remain=(double) W;
			for(k=0;k<W;k++)
			{
				(*clustIter).power[k]=true;//1;
			}//end of k
		}
//		for(motifsIDIter=(*clustIter).motifs.begin();motifsIDIter !=(*clustIter).motifs.end();motifsIDIter++)
//			clust[j].member[k]=-1;
	}//end of j
/*
	for(j=0;j<numOfCluster;j++)
	{
		count=0;
		for(k=0;k<number_of_motif.total;k++)
		{
			if((motif[k].cluster==j)&&((bin==3)||(motif[k].oddeven==bin)))
			{
				clust[j].member[count]=k;
				count++;
			}//end of if
		}//end of k
	}//end of j
*/
}//end of initWidth

double metEven(const list <clustClass>::iterator clustIter, short indel, short position, short W)
{
	int m;
	double oldProfile, newProfile;
	double sumin, sumex, proba;
	double ran;//,decide,partialsum;
	int suminclust,la;
	int na,nb,nc;
	double dif;
	double tem;
	list <int>::iterator motifsIDIter;
	int shift;

	/* seed random number generator */
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=19;

	tem=gammaln(betasum)-sumbeta;
	suminclust=0;
	for(m=0;m<4;m++)
		suminclust=suminclust+(*clustIter).profile[m][0];
	suminclust=2*suminclust;
	sumex=0;
	for(m=0;m<4;m++)
	{
		sumex=sumex+gammaln((*clustIter).profile[m][position]
						+(*clustIter).profile[kstar[m]][W-1-position]+beta[m]);
	}
	oldProfile=sumex-gammaln((double) suminclust+betasum);
	newProfile=0;
	for(motifsIDIter=(*clustIter).motifsID.begin();motifsIDIter !=(*clustIter).motifsID.end();motifsIDIter++)
	{
		sumin=2*allMotifs[(*motifsIDIter)].depth;
		sumex=0;
		shift=allMotifs[(*motifsIDIter)].motifshift;
		for(m=0;m<4;m++)
		{
			la=allMotifs[(*motifsIDIter)].profile[m][position+shift]
			  +allMotifs[(*motifsIDIter)].profile[kstar[m]][W-position-1+shift];
			sumex=sumex+gammaln((double) la+beta[m]);
		}
		newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
	}//end of k
	dif=oldProfile-newProfile-((*clustIter).motifsID.size()-1)*tem;
	//proba=exp(dif);// changed 08/19/05
	proba=1-1/(1+exp(dif));// changed 08/19/05
	ran=unran(&na, &nb, &nc);
	if((indel==1)&&(proba>ran))
	{
		return dif;
	}
	else if((indel==-1)&&((1.0-proba)>ran))//1.0/proba
	{
		return -dif;
	}//end of else if
	return 0;
}//end of metEven

double metOdd(const list <clustClass>::iterator clustIter, short indel, short position, short W)
{
	int m;
	double oldProfile, newProfile;
	double sumin, suminhalf;
	double sumex, sumexhalf, proba;
	double ran;//,decide,partialsum;
	int suminclust,suminclusthalf;
	int la, lahalf;
	int na,nb,nc;
	double dif;
	double tem,temhalf;
	list <int>::iterator motifsIDIter;
	int shift,motifsSize;

	/* seed random number generator */
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=19;

	motifsSize=(*clustIter).motifsID.size();
	tem=gammaln(betasum)-sumbeta;
	temhalf=gammaln(betasumhalf)-sumbetahalf;
	suminclusthalf = 0;
	for(m=0;m<4;m++)
		suminclusthalf=suminclusthalf+(*clustIter).profile[m][0];
	suminclust=2*suminclusthalf;
	if(position==(W/2))
	{
		sumex=0;
		for(m=0;m<4;m++)
		{
			sumex=sumex+gammaln((*clustIter).profile[m][position]+0.5*beta[m]);
		}
		oldProfile=sumex-gammaln((double) suminclusthalf+betasumhalf);
		newProfile=0;
		for(motifsIDIter=(*clustIter).motifsID.begin();motifsIDIter !=(*clustIter).motifsID.end();motifsIDIter++)
		{
			suminhalf=allMotifs[(*motifsIDIter)].depth;
			sumexhalf=0;
			shift=allMotifs[(*motifsIDIter)].motifshift;
			for(m=0;m<4;m++)
			{
				lahalf=allMotifs[(*motifsIDIter)].profile[m][position + shift];
				sumexhalf=sumexhalf+gammaln((double) lahalf+0.5*beta[m]);
			}
			newProfile=newProfile+sumexhalf-gammaln((double) suminhalf+betasumhalf);
		}//end of motifsIDIter
		dif=oldProfile-newProfile-(motifsSize-1)*temhalf;	
	}//end of if
	else
	{
		sumex=0;
		for(m=0;m<4;m++)
		{
			sumex=sumex+gammaln((*clustIter).profile[m][position]
						+(*clustIter).profile[kstar[m]][W-1-position]+beta[m]);
		}
		oldProfile=sumex-gammaln((double) suminclust+betasum);
		newProfile=0;
		for(motifsIDIter=(*clustIter).motifsID.begin();motifsIDIter !=(*clustIter).motifsID.end();motifsIDIter++)
		{
			sumin=2*allMotifs[(*motifsIDIter)].depth;
			sumex=0;
			shift=allMotifs[(*motifsIDIter)].motifshift;
			for(m=0;m<4;m++)
			{
				la=allMotifs[(*motifsIDIter)].profile[m][position+shift]
				  +allMotifs[(*motifsIDIter)].profile[kstar[m]][W-position-1+shift];
				sumex=sumex+gammaln((double) la+beta[m]);
			}
			newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
		}//end of motifsIDIter
		dif=oldProfile-newProfile-(motifsSize-1)*tem;
	}//end of else
	//proba=exp(dif);// changed 08/19/05
	proba=1-1/(1+exp(dif));// changed 08/19/05
	ran=unran(&na, &nb, &nc);
	if((indel==1)&&(proba>ran))
	{
		return dif;
	}
	else if((indel==-1)&&((1.0-proba)>ran))//1.0/proba
	{
		return -dif;
	}//end of else if
	return 0;
}//end of metroOdd

double metSingle(const list <clustClass>::iterator clustIter, short indel, short position, short W)
{
	int m;
	double oldProfile, newProfile;
	double sumin, sumex, proba;
	double ran;//,decide,partialsum;
	int suminclust,la;
	int na,nb,nc;
	double dif;
	double tem;
	list <int>::iterator motifsIDIter;
	int shift;

	/* seed random number generator */
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=19;

	tem=gammaln(betasum)-sumbeta;
	suminclust=0;
	for(m=0;m<4;m++)
		suminclust=suminclust+(*clustIter).profile[m][0];
	sumex=0;
	for(m=0;m<4;m++)
	{
		sumex=sumex+gammaln((*clustIter).profile[m][position]+beta[m]);
	}
	oldProfile=sumex-gammaln((double) suminclust+betasum);
	newProfile=0;
	for(motifsIDIter=(*clustIter).motifsID.begin();motifsIDIter !=(*clustIter).motifsID.end();motifsIDIter++)
	{
		sumin=allMotifs[(*motifsIDIter)].depth;//removed 2*
		sumex=0;
		shift=allMotifs[(*motifsIDIter)].motifshift;
		for(m=0;m<4;m++)
		{
			la = allMotifs[(*motifsIDIter)].profile[m][position+shift];
			  //+(*motifsIDIter).profile[kstar[m]][W-position-1+shift+start];
			sumex=sumex+gammaln((double) la+beta[m]);
		}
		newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
	}//end of k
	dif=oldProfile-newProfile-((*clustIter).motifsID.size()-1)*tem;
	//proba=exp(dif);// changed 08/19/05
	proba=1-1/(1+exp(dif));//changed 08/19/05
	ran=unran(&na, &nb, &nc);
	if((indel==1)&&(proba>ran))
	{
		return dif;
	}
	else if((indel==-1)&&((1.0-proba)>ran))//1.0/proba
	{
		return -dif;
	}//end of else if
	return 0;
}//end of metSingle

double onoffMetro(const list <clustClass>::iterator clustIter, int W, int WMIN, short bin)
{
	int na,nb,nc;
	double dif,ran,lower;
	int current, leftMargin,rghtMargin;

	/* seed random number generator */
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=50;

	dif=0;
	ran=unran(&na, &nb, &nc);
	if(bin==0)
	{
		lower=(double) WMIN/2;
		current=W/2;
	}
	else if(bin==1)
	{
		lower=(double)0.5*WMIN;
		current=W/2+1;
	}
	else
	{
		lower=(double) WMIN;
		current=W;
	}
	leftMargin=(*clustIter).left-1;
	rghtMargin=(*clustIter).right;
	if(ran<0.5)
	{
		if((ran<0.25)&&(leftMargin>1))// whether add a nt to the left
		{
			if(bin==0)
				dif=metEven(clustIter,1,leftMargin,W);
			else if(bin==1)
				dif=metOdd(clustIter,1,leftMargin,W);
			else if(bin==2)
				dif=metSingle(clustIter,1,leftMargin,W);
			if(fabs(dif)>0.00001)
			{
				(*clustIter).power[leftMargin]=true;//1;
				(*clustIter).remain++;
				(*clustIter).left=leftMargin;
			}
		}
		else if((ran>0.25)&&((*clustIter).remain>lower))// whether remove a nt in the left
		{
			if(bin==0)
				dif=metEven(clustIter,-1,leftMargin+1,W);
			else if(bin==1)
				dif=metOdd(clustIter,-1,leftMargin+1,W);
			else 
				dif=metSingle(clustIter,-1,leftMargin+1,W);
			if(fabs(dif)>0.00001)
			{
				(*clustIter).power[leftMargin+1]=false;//0;
				(*clustIter).remain--;
				(*clustIter).left=leftMargin+2;
			}
		}
	}//end of if
	else
	{
		if((ran<0.75)&&((*clustIter).right>0))//whether add a nt to the right
		{
			if(bin==0)
				dif=metEven(clustIter,1,current-rghtMargin,W);
			else if(bin==1)
				dif=metOdd(clustIter,1,current-rghtMargin,W);
			else 
				dif=metSingle(clustIter,1,current-rghtMargin,W);
			if(fabs(dif)>0.00001)
			{
				(*clustIter).power[current-rghtMargin]=true;//1;
				if((bin==1)&&(rghtMargin==1))
					(*clustIter).remain=(*clustIter).remain+0.5;
				else
					(*clustIter).remain++;
				(*clustIter).right=rghtMargin-1;
			}
		}
		else if((ran>0.75)&&((*clustIter).remain>lower))//whether remove a nt to the right
		{
			if(bin==0)
				dif=metEven(clustIter,-1,current-rghtMargin-1,W);
			else if(bin==1)
				dif=metOdd(clustIter,-1,current-rghtMargin-1,W);
			else 
				dif=metSingle(clustIter,-1,current-rghtMargin-1,W);
			if(fabs(dif)>0.00001)
			{
				(*clustIter).power[current-rghtMargin-1]=false;//0;
				if((bin==1)&&(rghtMargin==0))
					(*clustIter).remain=(*clustIter).remain-0.5;
				else
					(*clustIter).remain--;
				(*clustIter).right=rghtMargin+1;
			}
		}
	}//end of else
	return(dif);
}//end of onoffMetro

double fragment(const short bin, const int WMIN, const int WMAX, const int W)
{
	double dif;
	list <clustClass>::iterator clustIter;

	dif=0;
	for(clustIter = (allClusts).begin();clustIter != (allClusts).end();clustIter ++)
	{			
		dif=dif+onoffMetro(clustIter,W,WMIN,bin);
		//int k;
		//for(k=0;k<W;k++)
		//	cout << (*clustIter).power[k];
		//cout <<endl;
	}
	return dif;
}//end of fragmentconst

double clustRemoveMotif(const int motifOrder, const list <clustClass>::iterator clustPoint, const short stage)
{
//	int oldwidth = 0;
//	int newwidth = 0;
	double likenew,likeold;

//	oldwidth = (*clustPoint).width;
	likeold = (*clustPoint).loglikelihood(0);
	(*clustPoint).minusMotif(motifOrder);
	likenew = (*clustPoint).loglikelihood(0);
	return(likenew-likeold);
	cout <<"remove: old= "<<likeold<<"new= "<<likenew<<endl;
}// end of clustRemoveMotif

double clustInsertMotif(const int motifOrder, const int decision, const short stage)
{
	int count;
	list <clustClass>::iterator clustIter;
	double likenew,likeold;
	int clustWidth = 0;
	int newwidth = 0;
			
	count=0;
	for(clustIter = allClusts.begin();clustIter != allClusts.end();clustIter ++)
	{
		if(count==decision)
		{
			likeold = (*clustIter).loglikelihood(0);
			(*clustIter).addMotif( motifOrder);
			likenew = (*clustIter).loglikelihood(0);
// assign this cluster to motifHome
			motifHome[motifOrder] = clustIter;
//218			cout <<"insert: old= "<<likeold<<"new= "<<likenew<<endl;
			return(likeold-likenew);		
		}//end of if
		else
			count++;
	}//end of clustIter
	return(-1);//wierd result
}// end of clustInsertMotif

int shiftdecide(vector <vector <double> > &ratioAll,const double FACTOR_0, int &shift, double &postprob)
{
	int j,k;
	double sum, insum;
	int number;
	int na,nb,nc;
	double ran;
	int decide = 0;
	double compare, partialsum,accumu;

	//seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=50;

	sum=FACTOR_0;
	number=ratioAll.size();
	for(j=0;j<number;j++)
	{
		for(k=0;k<(ratioAll[j]).size();k++)				
			sum=sum+ratioAll[j][k];
	}
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		shift = 0;
		postprob = FACTOR_0/sum;
		return(0);
	}
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		insum = 0;
		for(k=0;k<ratioAll[j].size();k++)
			insum = insum + ratioAll[j][k];
		partialsum = partialsum + insum;
		if(partialsum > compare)
		{
			decide = j;
			postprob = insum / sum;
			break;
		}//end of if
	}//end of j
	//to decide on shift 08/02/05
	if(ratioAll[decide].size() ==1)
	{
		shift = 0;
		return (decide + 1);
	}//end of if
	else
	{
		ran=unran(&na, &nb, &nc);
		compare = insum *ran;
		accumu = 0;
		for(k=0;k<ratioAll[decide].size();k++)
		{
			accumu = accumu + ratioAll[decide][k];
			if(accumu > compare)
			{	
				shift = k;
				return (decide + 1);				
			}//end of if
		}//end of k
	}//end of else
	cout << "errors"<<endl;
	exit(0);
}//end of shiftdecide

/*
int decide(vector <double> &ratioVec, const double FACTOR_0, double &postprob)
{
	int j;
	double sum;
	int number;
	int na,nb,nc;
	double ran;
	int decide = 0;
	double compare, partialsum;

// seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=50;

	sum=FACTOR_0;
	number=ratioVec.size();
	for(j=0;j<number;j++)
		sum=sum+ratioVec[j];
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	if (FACTOR_0 > compare)
	{
		postprob = FACTOR_0/sum;
		return(0);
	}//end of if
	partialsum = FACTOR_0;
	for(j=0;j<number;j++)
	{
		partialsum = partialsum + ratioVec[j];
		if(partialsum > compare)
		{
			decide = j + 1;
			postprob = ratioVec[j]/sum;
			return (decide);
		}//end of if
	}//end of j
	cout << "errors"<<endl;
	exit(0);
}//end of decide
*/

double emptyLoglikelihoodupdateSingle(const int motifOrder, const list <clustClass>::iterator allClustsIter,
						    const short stage)
{
	short ** newProfile;
	short j,k;
	double logsum;//,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
	short lsummotif,lb;
	double difference,tem;
	int start;
	int width;

	width = (*allClustsIter).width;
	newProfile = new short*[4];
	for(j=0;j<4;j++)
		newProfile[j] = new short[width];
	start= allMotifs[motifOrder].motifshift;
	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	//lsum=0.0;
	motifsum=0;
	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
	for (j=0;j<width;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+start];
				newProfile[k][j] = allMotifs[motifOrder].profile[k][j+start];//***shift to start//remove (*allClustsIter).profile[k][j]
			}//end of k
		}//end of if
	}//end of j
	for(j=0;j<width;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				logsum=logsum+gammaln(beta[k]);
				lanew=newProfile[k][j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		}//end of if
		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	for(j=0;j<4;j++)
		delete [] newProfile[j];
	delete [] newProfile;
	return (-difference);
}//end of emptyLoglikelihoodupdateSingle

//-----
// calculate likelihood for clusters of single motifs
//-----
double loglikelihoodSingle(const short stage)
{
	int j,k;
	double tem;
	int lsum,sumin;
	list <clustClass>::iterator allClustsIter;
	list <int>::iterator motifsIDIter;
	double logsum,la;
	int start;
	int motifsSize;
	double logsumold =0;

	tem=gammaln(betasum)-sumbeta;
	logsum=0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end();allClustsIter ++)
	{
/*		for(j=0;j<4;j++)
		{
			for(k=0;k<(*allClustsIter).width;k++)
				cout << (*allClustsIter).profile[j][k];
			cout <<"+++++"<<endl;
		}	
*/
		motifsSize = (*allClustsIter).motifsID.size();
		lsum=0;
 		for(j=0;j<4;j++)
			lsum=lsum+(*allClustsIter).profile[j][0];
		//lsum=2*lsumhalf;//removed 08/12/05
		int a = (*allClustsIter).width;
		for(j=0;j<(*allClustsIter).width;j++)
		{
			if((stage!=1)||((*allClustsIter).power[j]))
			{
				for(k=0;k<4;k++)
				{
					int a=(*allClustsIter).profile[k][j];
					la = (double) (*allClustsIter).profile[k][j];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) lsum + betasum);
			}//end of if
			else if ((stage==1)&&(!(*allClustsIter).power[j]))
			{
				for(motifsIDIter = (*allClustsIter).motifsID.begin();motifsIDIter != (*allClustsIter).motifsID.end();motifsIDIter ++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
					start = allMotifs[(*motifsIDIter)].motifshift;// add (*motifsIDIter).shift 
					for(k=0;k<4;k++)
					{
						la = (double) allMotifs[(*motifsIDIter)].profile[k][j+start];//+beta[k];
						logsum=logsum+gammaln(la+beta[k]);
					}//end of k
					logsum=logsum-gammaln((double) sumin + betasum);
				}//end of 
				logsum=logsum+(double) (motifsSize-1)*tem;
			}//end of else if
		}//end of j
		logsum = logsum + (double) (*allClustsIter).width *tem;
//218		cout <<"cluster like="<<logsum-logsumold<<endl;
        logsumold = logsum;
	}//end of allClustsIter
//218	cout <<endl;
	//logsum=logsum+allClusts.size()*(double) width *tem;
	return(logsum);
}//end of loglikelihoodSingle

double bayesRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart, 
						const short stage)
{
	int j,k;
	int Lold,Lz,Lnew,sum;
	double sumbayes,tem;
	double newProfile;
	double bayesRatio;
	int a,b;
	int width;
	
	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);//changed 9/13/05
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum+(*allClustsIter).profile[j][0];
	}
	Lold=sum;
	Lz=allMotifs[motifOrder].depth;
	Lnew=Lold+Lz;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
		if((stage !=1)||((*allClustsIter).power[j]))//removed ==1 08/10/05
		{
			for(k=0;k<4;k++)
			{
				a=(*allClustsIter).profile[k][j+clustStart];
				b=allMotifs[motifOrder].profile[k][j+motifStart];
				newProfile=(*allClustsIter).profile[k][j+clustStart] + allMotifs[motifOrder].profile[k][j+motifStart];
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) (*allClustsIter).profile[k][j+clustStart]+beta[k])
					-gammaln((double) allMotifs[motifOrder].profile[k][j+motifStart]+beta[k]);
			}//end of k
		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum)+gammaln((double) Lz+betasum)
		-gammaln(betasum)+sumbeta;
	if(stage == 1)
		bayesRatio = sumbayes+tem*(double) (*allClustsIter).remain;
	else 
		bayesRatio = sumbayes+tem*(double) width;
	return(bayesRatio);
}//end of bayesRatioSingle

double loglikelihoodUpdateSingle (const int motifOrder, const list <clustClass>::iterator allClustsIter,
						    const short stage)
{
	short ** newProfile;
	short j,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
	short lsummotif,lb;
	double difference,tem;
	int width;
	int start;
	
	width = (*allClustsIter).width;
	// need to find out if the cluster width changed during the add/removal of the motif
	newProfile = new short*[4];
	for(j=0;j<4;j++)
		newProfile[j] = new short[width];
	start = allMotifs[motifOrder].motifshift;
	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
	motifsum=0;
	for (j=0;j<4;j++)
	{
		lsum=lsum+(*allClustsIter).profile[j][0];		   		
	}
	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=lsum+lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);

	for (j=0;j<width;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+start];
				newProfile[k][j]=(*allClustsIter).profile[k][j]
							 +allMotifs[motifOrder].profile[k][j+start];//***shift to start
			}//end of k
		}//end of if
	}//end of j
	for(j=0;j<width;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				la=(*allClustsIter).profile[k][j];
				logsum=logsum+gammaln((double)la+beta[k]);
				lanew=newProfile[k][j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		}//end of if
		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	for(j=0;j<4;j++)
		delete [] newProfile[j];
	delete [] newProfile;
	return difference;
}//end of loglikelihoodupdateSingle

/*
double diffWidthLoglikelihoodUpdateSingle (const int oldwidth, const int newwidth, const int motifOrder, 
										   const list <clustClass>::iterator allClustsIter,const short stage)
{
	short ** newProfile;
	short j,k;
	double logsum,la,lsum;
	double logsumold,lanew, lsumold;
	double logsumnew, lsumnew;
	double motifsum,sumin;
	short lsummotif,lb;
	double difference,tem;
	int clustwidth;
	int start;
	int motifwidth = 0;
	int motifsSize = 0;
	int clustshift, motifshift;
	list <motifClass>::iterator motifsIDIter;
	
	clustshift = allMotifs[motifOrder].clustshift;
	motifshift = allMotifs[motifOrder].motifshift;
	newProfile = new short*[4];
	for(j=0;j<4;j++)
		newProfile[j] = new short[newwidth];
	start = allMotifs[motifOrder].motifshift;
	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
	motifsum=0;
	for (j=0;j<4;j++)
	{
		lsum=lsum+(*allClustsIter).profile[j][0];		   		
	}
	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=lsum+lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);

	if(oldwidth < newwidth)//the only possibility is that just one motif is left, and a long one.
	{
// likelihood for new cluster
		lsum=0;
		logsum = 0;
 		for(j=0;j<4;j++)
			lsum=lsum+(*allClustsIter).profile[j][0];
		for(j=0;j<newwidth;j++)
		{
			for(k=0;k<4;k++)
			{
				//int a=(*allClustsIter).profile[k][j];
				la = (double) (*allClustsIter).profile[k][j];
				logsum=logsum+gammaln(la+beta[k]);
			}//end of k
			logsum=logsum-gammaln((double) lsum + betasum);
		}//end of j
		logsum = logsum + (double) newwidth *tem;
// likelihood for old cluster
		logsumold = 0;
		lsumold = lsum + allMotifs[motifOrder].depth;
		for(j=0;j<oldwidth;j++)
		{
			if((stage!=1)||((*allClustsIter).power[j]))
			{
				for(k=0;k<4;k++)
				{
					la=(*allClustsIter).profile[k][j+start] + allMotifs[motifOrder].profile[k][j];
					logsumold = logsumold + gammaln((double)la+beta[k]);
//					printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
				}//end of k
				logsumold = logsumold - gammaln(lsumold+betasum);
//				printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
			}//end of if
			else if((stage==1)&&(!(*allClustsIter).power[j]))
			{
				for(k=0;k<4;k++)
				{
					lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
					motifsum=motifsum+gammaln((double) lb+beta[k]);
				}//end of k
				motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
			}//end of else if
		}//end of j	

		difference=logsumnew-logsum+motifsum;
	}//end of if
	else if(oldwidth > newwidth)//a shofter motif is added to the cluster
	{

	}
	else
	{
		cout <<"error, old and new widths are supposed to be different."<<endl;
	}
	for (j=0;j<clustshift;j++)// the part only old profile exist
	{
		for (k=0;k<4;k++)
		{
			//a=allMotifs[motifOrder].profile[k][j+start];
			newProfile[k][j]=0;
		}//end of k
	}//end of j
	for (j=clustshift;j<clustshift+newwidth;j++)//the part both old and new profile exist
	{
		if((stage!=1)||((*allClustsIter).power[j+clustshift]))
		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+start];
				newProfile[k][j]=(*allClustsIter).profile[k][j]
						 +allMotifs[motifOrder].profile[k][j];
			}//end of k
		}//end of if
	}//end of j
	for (j=newwidth+clustshift;j<clustwidth;j++)// the part only old profile exist
	{
		for (k=0;k<4;k++)
		{
			//a=allMotifs[motifOrder].profile[k][j+start];
			newProfile[k][j]=0;
		}//end of k
	}//end of j

// likelihood for this cluster
	motifsSize = (*allClustsIter).motifs.size();
	lsum=0;
	logsum = 0;
 	for(j=0;j<4;j++)
		lsum=lsum+(*allClustsIter).profile[j][0];
	int a = (*allClustsIter).width;
	for(j=0;j<clustwidth;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				//int a=(*allClustsIter).profile[k][j];
				la = (double) (*allClustsIter).profile[k][j];
				logsum=logsum+gammaln(la+beta[k]);
			}//end of k
			logsum=logsum-gammaln((double) lsum + betasum);
		}//end of if
		else if ((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(motifsIDIter = (*allClustsIter).motifs.begin();motifsIDIter != (*allClustsIter).motifs.end();motifsIDIter ++)
			{
				sumin=(*motifsIDIter).depth;
				start=(*motifsIDIter).motifshift;// add (*motifsIDIter).shift 
				for(k=0;k<4;k++)
				{
					la = (double) (*motifsIDIter).profile[k][j+start];//+beta[k];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) sumin + betasum);
			}//end of motifsIDIter
			logsum=logsum+(double) (motifsSize-1)*tem;
		}//end of else if
	}//end of j
	logsum = logsum + (double) (*allClustsIter).width *tem;
// likelihood for updated cluster
	logsumnew = 0;
	for(j=0;j<newwidth;j++)
	{
		if((stage!=1)||((*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				la=(*allClustsIter).profile[k][j];
				logsum=logsum+gammaln((double)la+beta[k]);
				lanew=newProfile[k][j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		}//end of if
		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	for(j=0;j<4;j++)
		delete [] newProfile[j];
	delete [] newProfile;
	return difference;
}//end of diffWidthloglikelihoodupdateSingle
*/

//----- read in motif block data and put into map and list -----//
int readSeq(vector <motifClass> &fullMotifVector, const char inputFileName[])
{
	int j;
	list <motifClass> clust;
	list < clustClass >::iterator allClustsIter;
	list < motifClass >::iterator motifsIDIter;
	motifClass temMotif;
	clustClass tempClust;
	vector <string> seq;
//	const string inputFilename = inputFileName;//"fileIn1.txt";//fileIn1
//	const string outputFilename = "newout.txt";
	ifstream inFile;
	string lineString;
	istringstream iss;
	string motif;
	char c;
	int rowCount;
	int motifCount = 0;
	string tempName;
	ofstream outputFile;
	int signal=0;//signal for if it is the first line of the sfirst motif
	int minwidth = 0;
	int maxwidth = 0;

	cout << "start reading in data" << endl;
	inFile.open(inputFileName);//.c_str());
	if(!inFile)
	{
		cout << "ERROR: Unable to open file: " << inputFileName << endl;
	    exit(2);
	}//end of if
// start reading in data
	motifCount=0;
	//temMotif.name=" ";
	signal=0;
	while (inFile)
	{
// get name
		getline(inFile,lineString);
// get DNA sequence
		iss.clear();
		iss.str(lineString + " ");
		iss >> c;
		if(!iss) // at the end of a motif just read in, summerize information
		{
			//cout << "blank line" <<endl;
//			temMotif.cluster = motifCount % INITIAL_CLUSTER;
			temMotif.motifshift = 0;
			temMotif.width=motif.length();
			temMotif.depth=rowCount;
			if(temMotif.width %2 ==0)
				temMotif.isEven = true;
			else
				temMotif.isEven =  false;
			temMotif.basepair=seq;
			seq.clear();
			temMotif.genProfile();
			fullMotifVector.push_back(temMotif);
			motifCount++;
		}
		else if(c=='>') // first line of FASTA format, take name
		{
			iss >> tempName;
			if((signal==0)||(tempName!=temMotif.name))
			{
				signal=1;
				temMotif.id=motifCount;
				temMotif.name=tempName;
				rowCount=0;
				iss.ignore();
			}
		}
		else // subsequent lines of FASTA files, take sequences
		{
			iss.clear();
			iss.str(lineString + " ");
			iss >> motif;
			seq.push_back(motif);
			rowCount++;
		}
	}//end of while
	inFile.close();
// print information
/*
	outputFile.open(outputFilename.c_str());
	if(!outputFile)
	{
		cout << "ERROR: Unable to open file: " << outputFilename << endl;
	    exit(2);
	}//end of if
	j = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		cout << "cluster " << j << endl;
		j++;
		(*allClustsIter).printClust(clustWidth, outputFile);
		for(motifsIDIter = (*allClustsIter).motifs.begin();motifsIDIter != (*allClustsIter).motifs.end(); motifsIDIter++)
		{
			//(*motifsIDIter).clustIter = allClustsIter;
			//motifHome[(*motifsIDIter).id]=allClustsIter;
			cout << "name= " << (*motifsIDIter).name << endl;
		}// end of motifsIDIter
	}// end of allClustsIter
	outputFile.close();
*/
	minwidth = fullMotifVector[0].width;
	maxwidth = minwidth;
	for(j=0;j<fullMotifVector.size();j++)
	{
		if((j>0)&&(fullMotifVector[j].width > maxwidth))
		{
			maxwidth = fullMotifVector[j].width;
		}
		else if((j>0)&&(fullMotifVector[j].width < minwidth))
		{
			cout << minwidth<<endl;
			minwidth = fullMotifVector[j].width;
		}
		//cout << "motif " << j << endl; 		
		//cout <<"width = "<< fullMotifVector[j].width<<endl;
		//fullMotifVector[j].printSeq();
	}
	cout << "total sites" << fullMotifVector.size()<<endl;
	cout <<"min width = "<< minwidth<<endl;
	cout <<"max width = "<< maxwidth<<endl;
	//exit(0);
	return(motifCount);
}//end of readSeq

void initialCluster(const short bin, const int motifCount, const vector <motifClass> fullMotifVector)
{
	int j;
	list <motifClass> clust;
	list < clustClass >::iterator allClustsIter;
	list < motifClass >::iterator motifsIDIter;
	clustClass tempClust;
	int motifsSize;
//	int *width, group;

	for(j=0;j<motifCount;j++)
	{
		if((bin ==2)||((bin ==0)&&(fullMotifVector[j].isEven))||((bin ==1)&&(!fullMotifVector[j].isEven)))
			allMotifs.push_back(fullMotifVector[j]);
	}//end of j
	motifsSize = allMotifs.size();
//	width = new int[INITIAL_CLUSTER];
//	for(j=0;j<INITIAL_CLUSTER;j++)// initialize clusters
//		width[j] = allMotifs[j].width;
//	for(j=INITIAL_CLUSTER;j<motifsSize;j++)
//	{
//		group = j%INITIAL_CLUSTER;
//		if(allMotifs[j].width < width[group])
//			width[group] = allMotifs[j].width;
//	}//end of j
	for(j=0;j<INITIAL_CLUSTER;j++)// initialize clusters
	{
		allMotifs[j].clustshift =0;
		tempClust.motifsID.push_back(j);
//		tempClust.width = width[j];
		tempClust.width = allMotifs[j].width;
		tempClust.iniProfile(bin);
		int b = tempClust.profile[0][0];
		allClusts.push_back(tempClust);
		tempClust.motifsID.clear();
	}//end of j
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		motifHome.push_back(allClustsIter);
	//allClustsIter = allClusts.begin();
	//(*allClustsIter).addMotif(allMotifs, INITIAL_CLUSTER, clustWidth);
	//motifHome.push_back(allClustsIter);
	for(j=INITIAL_CLUSTER;j<motifsSize;j++)
	{
		allMotifs[j].clustshift =0;
		//group = j % INITIAL_CLUSTER;
		if ((j % INITIAL_CLUSTER) == 0)// set allClusts pointer to the beginning
		{
			allClustsIter = allClusts.begin();
			int a = (*allClustsIter).profile[0][0];
			(*allClustsIter).addMotif( j);
			motifHome.push_back(allClustsIter);
		}
		else // add motifs to existing clusters
		{
			allClustsIter ++;
			allMotifs[j].clustshift =0;		
			(*allClustsIter).addMotif(j);
			motifHome.push_back(allClustsIter);
		}
	}//end of if
//	delete []width;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		(*allClustsIter).printClust();
}//end of initialCluster

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
