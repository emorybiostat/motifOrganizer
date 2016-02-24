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
//#define INITIAL_CLUSTER 200
//#define IN_WIDTH 16//16//6//2
#define MAX_LENGTH 100
#define HEAPSORT class sortclustClass
#define KEY postprob
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

using namespace std;

double gammaln(double xx);
static double beta[4]={1.0,1.0,1.0,1.0};
static double betasum=beta[0]+beta[1]+beta[2]+beta[3];
static double sumbeta=gammaln(beta[0])+gammaln(beta[1])+gammaln(beta[2])+gammaln(beta[3]);
static double betasumdif = gammaln(betasum)-sumbeta;
//static double betasumhalf=0.5*betasum;
//static double sumbetahalf=gammaln(0.5*beta[0])+gammaln(0.5*beta[1])+gammaln(0.5*beta[2])+gammaln(0.5*beta[3]);
//static short kstar[4]={3,2,1,0};//{1,0,3,2};

class motifClass
{
private:
public:
	string name;
	short width;
	int start;
	short depth;
	short motifshift;
	short clustshift;
	double ratio;
	double postprob;
//-01/19/07	void printSeq();
};

vector <motifClass> allMotifs;
int **allMotifProfiles;
double **allRatio;

class clustClass
{
private:
public:
	list <short> motifsID;
//-	vector <bool> power;
	short width;
///	short left;
///	short right;
///	double remain;
	vector <short> profile[4];
	void printClust();
	void iniProfile();//short bin);
	void addMotif(/*vector <motifClass> allMotifs,*/ short motifOrder);
	short minusMotif(/*vector <motifClass> allMotifs,*/ short motifOrder);
	void updateProfile(short sign, /*vector <motifClass> allMotifs,*/ short motifOrder);
	double loglikelihood();//const short stage);
	~clustClass()
	{
		motifsID.clear();
//-		power.clear();
		profile[0].clear();
		profile[1].clear();
		profile[2].clear();
		profile[3].clear();
	}
};

class sortclustClass
{
private:
public:
	short order;
	list <short> motifsID;
	short size;
///	short left;
///	short right;
	short width;
//-	vector <bool> power;
///	double remain;
	double postprob;
	~sortclustClass()
	{
		motifsID.clear();
//-		power.clear();
	}
};

//***** removed 01/19/07
/*
void motifClass::printSeq()
{
	short j;

	for(j=0;j<depth;j++)
		cout << basepair[j] <<endl;
}
*/
//***** removed 01/19/07 

void clustClass::printClust()
{
	short j,k;
	short count=1;//j,k;
	list <short>::iterator motifsIDIter;

	count=1;
	for(motifsIDIter = motifsID.begin(); motifsIDIter != motifsID.end(); motifsIDIter++)
	{
//		outputFile << "ID " << count;
//		outputFile << "Name " << (*motifsIDIter).name << " " <<(*motifsIDIter).width << 
//			" " <<(*motifsIDIter).depth << " " <<(*motifsIDIter).isEven << endl;
		//(*motifsIDIter).printSeq();
		cout << "ID " << count;
		cout << "Name " << allMotifs[(*motifsIDIter)].name << " " <<allMotifs[(*motifsIDIter)].width << 
			" " << allMotifs[(*motifsIDIter)].depth << " " << /*allMotifs[(*motifsIDIter)].isEven <<*/ endl;
		count++;
	}//end of motifsIDIter

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
//	outputFile << endl;
//	outputFile.close();
}//end of clustClass::printClust()

//----- initialize the profile matrix when the first motif join this cluster

void clustClass::iniProfile()//short bin)
{
	short j,k;
	list <short>::iterator motifsIDIter;
//	short motifWidth;

///	left=0;
///	right=0;
	
	motifsIDIter=motifsID.begin();
	//short a = (*motifsIDIter).start;
	//short b = (*motifsIDIter).shift; 
//	start=(*motifsIDIter).start + (*motifsIDIter).shift;//added (*motifsIDIter).shift 
	width=allMotifs[(*motifsIDIter)].width;
// ******* moved from before motifsIDIter
/*
	if(bin==2)//added 08/17/05
		remain = width;//added 08/17/05
	else//added 08/17/05
		remain=(double) width*0.5;
*/
// ******* moved from before motifsIDIter
// first set profile matrix to 0
	if(width<0)
	{
		cout <<"cluster width likely not specified yet."<<endl;
		exit(20);
	}
	for(j=0;j<width;j++)
	{
//-		power.push_back(1);
		for(k=0;k<4;k++)
			profile[k].push_back(0);
	}
	for(j=0;j<width;j++)
	{
		//power[j]=1;
		for(k=0;k<4;k++)
		{
	//		short b = (*motifsIDIter);
	//		short a = allMotifs[(*motifsIDIter)].profile[k][j];
			profile[k][j]=allMotifProfiles[k][allMotifs[(*motifsIDIter)].start + j];//+(*motifsIDIter).profile[kstar[k]][motifWidth-1-j-start];
		}
	}//end of j
}//end of clustClass::iniProfile

void clustClass::addMotif(/*vector <motifClass> allMotifs,*/ short motifOrder)
{
	motifsID.push_back(motifOrder);
	updateProfile(1, motifOrder);
//218	cout <<"shift= "<<allMotifs[0].motifshift<<endl;
}

short clustClass::minusMotif(/*vector <motifClass> allMotifs,*/ short motifOrder)
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

void clustClass::updateProfile(short sign, /*vector <motifClass> allMotifs,*/ short motifOrder)
{
	short j,k;
	list <short>::iterator motifsIDIter;
	short motifWidth, newwidth;
	short motifshift,clustshift;
	bool begin = true;

	motifshift = allMotifs[motifOrder].motifshift;
	clustshift = allMotifs[motifOrder].clustshift;
	motifWidth = allMotifs[motifOrder].width;
/* 02/18/06 ---*/
/* ********	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
		{
//-			cout << profile[j][k] << " ";//removed +clustshift
			if(profile[j][k]<0)
			{
				cout <<"in clust, motif="<<motifOrder<<endl;
				exit(0);
			}
		}
//-		cout <<endl;
	}
//-	cout <<endl;
	for(j=0;j<4;j++)
	{
		for(k=0;k<motifWidth;k++)
		{
//-			cout << allMotifs[motifOrder].profile[j][k]<<" ";//+motifshift] << " ";
			if(allMotifProfiles[j][allMotifs[motifOrder].start + k]<0)//+motifshift]<0)
			{
				cout <<"motif="<<motifOrder<<endl;
				exit(0);
			}
		}
//-		cout <<endl;
	}
//-	cout <<endl;
******** */
/*--- 02/18/06 */
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
				profile[j].push_back(allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k]);
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
				short e= allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
				short f = profile[j][k];
				profile[j][k]=profile[j][k]+sign * allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
/* *****
				if(profile[j][k] <0)
				{
					for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
					{
						cout <<(*motifsIDIter)<<" ";	
					}
					cout <<endl;
					//short e = profile[j][k];
					//short f = allMotifs[motifOrder].profile[j][k+motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
 						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}
***** */
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
				//short g = profile[j][k];
				//short h = allMotifs[motifOrder].profile[j][k+motifshift];
				profile[j][k]=profile[j][k] - allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
/* *****
				if(profile[j][k] <0)
				{
					short e = profile[j][k];
					short f = allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}//end of if
***** */
			}//end of k
			motifsIDIter = motifsID.begin();
			for(k=width;k<newwidth;k++)
			{
				//short g = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
				profile[j][k] = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
			}//end of k
		}//end of j
		cout <<"*****"<<endl;
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				cout <<profile[j][k]<<" ";
			cout <<endl;
		}
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				cout << allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k]<<" ";
			cout <<endl;
		}
		cout <<"*****"<<endl;
	}//end of else if
	else//need to subtract profile, since a shorter motif was added to the cluster
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				profile[j][k] = profile[j][k + clustshift] + allMotifProfiles[j][allMotifs[motifOrder].start + k];
			// need to update motif shift
		}//end of j
		for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
		{
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
			//short e = (*motifsIDIter);
			//short f = allMotifs[(*motifsIDIter)].motifshift;
			if((*motifsIDIter) != motifOrder)
				allMotifs[(*motifsIDIter)].motifshift = allMotifs[(*motifsIDIter)].motifshift + clustshift;
//				allMotifs[(*motifsIDIter).id].motifshift = allMotifs[(*motifsIDIter).id].motifshift + clustshift; 
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
		}//end of motifsIDIter
	}//end of else
	width = newwidth;
}//end of clustClass::updateProfile

double clustClass::loglikelihood()//const short stage)
{
	short j,k;
	short motifsSize;
	short la, lsum = 0;//, sumin,start;
	double logsum;//,tem;
	list <short>::iterator motifsIDIter;

//+	tem=gammaln(betasum)-sumbeta;
	motifsSize = motifsID.size();
	lsum=0;
	logsum =0;
 	for(j=0;j<4;j++)
		lsum = lsum + profile[j][0];
	for(j=0;j<width;j++)
	{
//-		if((stage!=1)||(power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				//+short a = profile[k][j];
				la = profile[k][j];
				logsum=logsum+gammaln((double) la+beta[k]);
			}//end of k
			logsum=logsum-gammaln((double) lsum + betasum);
//-		}//end of if
/*		else if ((stage==1)&&(!power[j]))
		{
			for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
			{
				sumin=allMotifs[(*motifsIDIter)].depth;
				start=allMotifs[(*motifsIDIter)].motifshift;// add (*motifsIDIter).shift 
				for(k=0;k<4;k++)
				{
					la = allMotifs[(*motifsIDIter)].profile[k][j+start];//+beta[k];
					logsum=logsum+gammaln((double) la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) sumin + betasum);
			}//end of 
			logsum=logsum+(double) (motifsSize-1)*tem;
		}//end of else if
*/
	}//end of j
	logsum = logsum + (double) width *betasumdif;//tem;
	return(logsum);
}//end of clustClass::loglikelihood

list <clustClass> allClusts;
vector <list <clustClass>::iterator> motifHome;

short readSeq(const char outPutFileName[], const short MAX_SHIFT, short &widthmin);

void initialCluster(const short INITIAL_SIZE);//const short motifCount, const short bin,

//void initialCluster(const short bin, const short motifCount, const vector <motifClass> fullMotifVector,const short INITIAL_CLUSTER);

double bmc(const int initialclusternum, const short numChains, const short numItera, //const short bin,
//+		   const short W, const short WMIN, const short WMAX, 
		   const short MAX_SHIFT, const double factor_0, list <clustClass> &finalClusters, vector <motifClass> &finalallMotifs);

double loglikelihoodSingle();//const short stage);

double chain(const short NUM_CYCLE, 
//+			 const short W, const short WMIN, const short WMAX, 
			 const short MAX_SHIFT, const double FACTOR_0, list <clustClass> &savedClusters, vector <motifClass> &savedallMotifs);
//const short bin,

double operation(const short motifOrder, const double FACTOR_0, const short MAX_SHIFT);//const short bin,
				 //short &del, const short stage);

//+double clustInsertMotif(const short motifOrder, const short decision);//, const short stage);

double clustInsertMotif(const short motifOrder, const list <clustClass>::iterator clustpoint);//, const short stage)

double clustRemoveMotif(const short motifOrder, const list <clustClass>::iterator clustpoint);//, const short stage);

double loglikelihoodUpdateSingle(const short motifOrder, const list <clustClass>::iterator allClustsIter);//, const short stage);

//double diffWidthLoglikelihoodUpdateSingle (const short oldwidth, const short newwidth, const short motifOrder, 
//										   const list <clustClass>::iterator allClustsIter,const short stage);

//double onoffMetro(const list <clustClass>::iterator allClustsIter, short W, short WMIN, short bin);

//+double metSingle(const list <clustClass>::iterator allClustsIter, short indel, short position, short W);

//double fragment(const short bin, const short WMIN, const short WMAX, const short W);

double bayesRatioSingle(short motifOrder, list <clustClass>::iterator allClustsIter,short motifStart, short clustStart);//, const short stage);

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,short motifStart, short clustStart);

double bayesSpecialRatioSingle(int motifOrder, list <short>::iterator motifsIDIter, short motifStart, short clustStart);

short decide(vector <double> &ratioVec, const double FACTOR_0, double &postprob);

short shiftdecide(vector <vector <double> > &ratioAll, const double FACTOR_0,short &shift, double &postprob);

double unran(int *na, int *nb, int *nc);

double emptyLoglikelihoodupdateSingle(const short motifOrder, const list <clustClass>::iterator allClustsIter);//, const short stage);

//+void initWidth(const short W);//const short bin,

//void resultDisplay(const short motifCount, const double loglike, const short bin, 
//		   const double factor_0, list <clustClass> clusters, const char outPutFileName[], 
//		   const char memberFileName[], const char shiftFileName[]);

void resultDisplay(const short motifCount, const double loglike, const double factor_0, list <clustClass> clusters, vector <motifClass> finalallmotifs,
                   const char inputFileName[], const char outPutFileName[], const char memberFileName[], const char shiftFileName[], const char widthFileName[], 
		   const short numChains, const short numIter, const short initialclusternum, const int MAX_SHIFT, const double diftime);
//                 const short numChains, const short numIter, const short initialclusternum, const double diftime);

void mqheap(HEAPSORT *p[], short total, short begin, short end);

static void sift(HEAPSORT *p[], short k, short j, short m);

//double finalRatioEven(list <clustClass>::iterator allClustsIter);

//double finalRatioOdd(list <clustClass>::iterator allClustsIter);

double finalRatioSingle(list <clustClass>::iterator allClustsIter);

//short checkWidth(	list <clustClass>::iterator clustpoint);

void genProfile(vector <string> seq, const int start, const int motifOrder, const short MAX_SHIFT);//const int width, const int depth);

int main(short argc, char **argv)
{
	short j;
	list <clustClass> finalClusters;
	vector <motifClass> finalallMotifs;
//	vector <motifClass> fullMotifVector;
//+	short W,WMIN,WMAX;
	short  MAX_SHIFT = 2;
	short numChains,numIter;
	double factor_0;
	short motifCount;
//+	short bin=0;
	double loglike;
//	short IN_WIDTH;
	char inputFileName[MAX_LENGTH]="a2.txt";//"1015.txt";//profile2.txt//90.txt
        char outputFileName[MAX_LENGTH]="out.txt";
//-        char motifFileName[MAX_LENGTH]="motifblock.txt";
	char memberFileName[MAX_LENGTH]="member.txt";
	char shiftFileName[MAX_LENGTH]="shift.txt";
        char widthFileName[MAX_LENGTH]="width.txt";
	short initialclusternum =10;
	short widthmin;
        time_t start,end;
        double dif;

//	IN_WIDTH = 6;//5;//20;
//	W=IN_WIDTH;//6;
//	WMIN=6;//5;//16;//16;//4;//less than WMAX if consider fragmentation
//	WMAX=W;//16;//6;
//	MAX_SHIFT=0;//5;
//	numChains=1;
//	numIter=5;
	factor_0=1.0;
//+	bin=2;//0;

	srand((unsigned)time(NULL));

        if(argc<10)
        {
          printf("9 options need to be specified:\n\tinput filename,\n\toutput filename,\n\tcluster member filename,\n\tmotif shift filename,\n\tcluster width filename,\n\tnumber of chains,\n\tnumber of iterations,\n\tinitial cluster size,\n\tmaximum allowed shift.\n");
            exit(0);
        }
        for(j=0;j<MAX_LENGTH;j++)
        {
            inputFileName[j]=argv[1][j];
            outputFileName[j]=argv[2][j];
//-            motifFileName[j]=argv[3][j];
	    memberFileName[j]=argv[3][j];
	    shiftFileName[j]=argv[4][j];
 	    widthFileName[j]=argv[5][j];
        }

	numChains = atoi(argv[6]);//1;
	numIter = atoi(argv[7]);//5;
	initialclusternum = atoi(argv[8]);//10;
	MAX_SHIFT = atoi(argv[9]);//2;
// read in data, initial assign clusters and calculate cluster profiles
	motifCount=readSeq(inputFileName,MAX_SHIFT, widthmin);
//+	W = widthmin;
//+	WMIN = W;
//+	WMAX = W;
//+++	initialCluster(initialclusternum);//motifCount, bin,
//	fullMotifVector.clear();//added 01/19/07
	//exit(0);
// clustering
//-	cout << "cluster size= " << allClusts.size();
        time (&start);	
	loglike = bmc(initialclusternum, numChains,numIter, /*W, WMIN,WMAX,*/ MAX_SHIFT,factor_0,finalClusters,finalallMotifs);//bin,
        time(&end);              
	dif = difftime (end,start) /60;
        allClusts.clear();//added 02/11/07
// output result
//	resultDisplay(motifCount, loglike, bin, factor_0, finalClusters,outputFileName, memberFileName,shiftFileName);
        resultDisplay(motifCount, loglike, factor_0, finalClusters, finalallMotifs, inputFileName, outputFileName, memberFileName,
                      shiftFileName, widthFileName, numChains, numIter, initialclusternum,MAX_SHIFT, dif);
	return(0);
}//end of main

//void resultDisplay(const short motifCount, const double loglike, const short bin, 
//		   const double factor_0, list <clustClass> clusters, const char outPutFileName[], 
//		   const char memberFileName[], const char shiftFileName[])
void resultDisplay(const short motifCount, const double loglike, const double factor_0, list <clustClass> clusters, vector <motifClass> finalallmotifs,
                   const char inputFileName[], const char outPutFileName[], const char memberFileName[], const char shiftFileName[], const char widthFileName[], 
		   const short numChains, const short numIter, const short initialclusternum, const int MAX_SHIFT, const double diftime)
{
	short j;//,k,m;
	list <clustClass>::iterator clustIter;
	list <short>::iterator motifsIDIter;
	ofstream outPutFile;
//	ofstream motifblockFile;// FASTA format motif sequence alignment file.
	ofstream memberFile;// cluster member motif information.
	ofstream shiftFile;// cluster motif profile alignment information.
        ofstream widthFile;// cluster width information.
	short finalClustNumber;
	short tot,count;
	double value =0;
	short motifsize;
	sortclustClass *clustarray, **sorted;
//	short start;
	short minWidth, maxWidth;
	bool begin;
        time_t rawtime;
        struct tm * timeinfo;

        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

	outPutFile.open(outPutFileName);
	if(!outPutFile)
	{
		cout << "ERROR: Unable to open file: " << outPutFileName << endl;
	    exit(3);
	}//end of if	
	memberFile.open(memberFileName);
	if(!memberFile)
	{
		cout << "ERROR: Unable to open file: " << memberFileName<<endl;
	    exit(4);
	}//end of if
        shiftFile.open(shiftFileName);   
        if(!shiftFile)    
        {
                cout << "ERROR: Unable to open file: " << shiftFileName<<endl;  
            exit(5);
        }//end of if
        widthFile.open(widthFileName);
        if(!widthFile)
        {
                cout << "ERROR: Unable to open file: " << widthFileName<<endl;
            exit(6);
        }//end of if

//***** deleted 01/19/07
	finalClustNumber = clusters.size();
// find out minimum and maximum cluster width
	begin = true;
	for(clustIter = clusters.begin();clustIter != clusters.end();clustIter ++)
	{
//+		(*clustIter).printClust();//added 03/20/07
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
	outPutFile << "*      BMC - VAR clustering Result     * \n";
	outPutFile << "*                                      * \n";
	outPutFile << "**************************************** \n\n";

        outPutFile <<"***************************************************"<<endl;
        outPutFile <<"The local time now is "<< asctime(timeinfo) <<endl;
        outPutFile <<"The following parameters are in effect:"<<endl;
        outPutFile <<"Input file name:           " <<inputFileName<<endl;
        outPutFile <<"Output file name:          " <<outPutFileName<<endl;
        outPutFile <<"Member file name:          " <<memberFileName<<endl;
        outPutFile <<"Shift position file name:  " <<shiftFileName<<endl;
        outPutFile <<"Motif width file name:     " <<widthFileName<<endl;
        outPutFile <<"Number of chains:          " <<numChains<<endl;
        outPutFile <<"Number of cycles:          " <<numIter<<endl;
        outPutFile <<"Number of initial clusters:" <<initialclusternum<<endl;
        outPutFile <<"Maximum allowed shift:	 " <<MAX_SHIFT<<endl;
        outPutFile <<"Values of alpha:           " <<factor_0<<endl;
        outPutFile <<"(Tunning parameter in DP model)."<<endl<<endl;
        outPutFile <<"It took BMC-VAR "<<diftime<<" minutes to finish this job."<<endl;
        outPutFile <<"***************************************************"<<endl<<endl;
	
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
///		clustarray[count].left = (*clustIter).left;
///		clustarray[count].right = (*clustIter).right;
///		clustarray[count].remain = (*clustIter).remain;
		clustarray[count].motifsID = (*clustIter).motifsID;
//		clustarray[count].power = (*clustIter).power;
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
		widthFile << (*sorted[j-1]).width<<endl;
		for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
		{
			tot++;
			//if((tot%2==1)&&(tot!=1)) 
			//	outPutFile << endl;                        ";
	        //+    outPutFile << tot << " " << allMotifs[(*motifsIDIter)].name << " [" << allMotifs[(*motifsIDIter)].motifshift 
		//+				<< "]  ("<< allMotifs[(*motifsIDIter)].postprob <<")"<<endl;
			outPutFile << tot << " " << finalallmotifs[(*motifsIDIter)].name << " [" << finalallmotifs[(*motifsIDIter)].motifshift
                                             << "]  ("<< finalallmotifs[(*motifsIDIter)].postprob <<")"<<endl;

			memberFile << (*motifsIDIter) <<"  ";
//			shiftFile <<allMotifs[(*motifsIDIter)].motifshift<< " ";
                        shiftFile <<finalallmotifs[(*motifsIDIter)].motifshift<< " ";
		}//end of motifsIDIter
		outPutFile <<endl<<endl;
		memberFile <<endl;
		shiftFile <<endl;
// ***** deleted 01/19/07
/*
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
*/
// ***** deleted 01/19/07
	}//end of j
	delete [] clustarray;
	delete [] sorted;
	outPutFile.close();
	memberFile.close();
	shiftFile.close();
	widthFile.close();
//-01/19/07	motifblockFile.close();
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

/*
double finalRatioEven(list <clustClass>::iterator allClustsIter)
{
    short j,k, start;
	list <short>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, sumex;
	short suminclust,la;
	double difsum;
	double tem;
	short width;

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
//-			if((*allClustsIter).power[j])//==1
//-			{
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
//-			}//end of if
		}//end of j
		difsum=difsum-(double)((*allClustsIter).remain)*(double)((*allClustsIter).motifsID.size()-1)*tem;
		return(difsum);
	}//end of else
}//end of finalRatioEven

double finalRatioOdd(list <clustClass>::iterator allClustsIter)
{
    short j,k, start;
	list <short>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, suminhalf;
	double sumex,sumexhalf;
	short suminclust,suminclusthalf;
	short la,lahalf;
	double difsum;
	double temhalf,tem;
	short motifsSize;
	short width;

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
//-			if((*allClustsIter).power[j])//==1
//-			{
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
//-			}//end of if
		}//end of j
//middle column		
//-		if((*allClustsIter).power[width/2])
//-		{
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
//-		}//end of if
		difsum=difsum-(double)((*allClustsIter).remain)*(double)(motifsSize-1)*tem;
//-		if((*allClustsIter).power[width/2])
//-			difsum=difsum-(double)(motifsSize-1)*temhalf;
		return(difsum);
	}//end of else
}//end of finalRatioOdd
*/

double finalRatioSingle(list <clustClass>::iterator allClustsIter)
{
    short j,k, begin;
	list <short>::iterator motifsIDIter;
//+  	double oldProfile, newProfile;
	double sumin, oldsumex,newsumex;
	short suminclust,la;
	double difsum;
//+	double tem;
	short motifsSize;
	short width;

	width = (*allClustsIter).width;
//+	tem=gammaln(betasum)-sumbeta;
	suminclust=0;
	for(j=0;j<4;j++)
		suminclust=suminclust+(*allClustsIter).profile[j][0];
//+	sumex=0;
	difsum=0;
	motifsSize = (*allClustsIter).motifsID.size();
	if(motifsSize==1)
	  return 0;
	else
	{
		for(j=0;j<width;j++)
		{
//-			if((*allClustsIter).power[j])//==1
//-			{
				oldsumex=0;
				for(k=0;k<4;k++)
				{	
					oldsumex = oldsumex+gammaln((*allClustsIter).profile[k][j]+beta[k]);
				}
				//+ oldProfile=sumex-gammaln((double) suminclust+betasum);
				//+ newProfile=0;
				newsumex =0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
					begin = allMotifs[(*motifsIDIter)].motifshift;
					//+ sumex=0;
					for(k=0;k<4;k++)
					{
						la = allMotifProfiles[k][allMotifs[(*motifsIDIter)].start + j +begin];
						newsumex = newsumex+gammaln((double) la+beta[k]);
					}//end of m
					//+ newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
					newsumex = newsumex-gammaln((double) sumin+betasum);
				}//end of k
				//+ difsum=difsum+oldProfile-newProfile;
				difsum=difsum+oldsumex-newsumex;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
//-			}//end of if
		}//end of j
		//+ difsum=difsum-(double)((*allClustsIter).width/*remain*/)*(double)(motifsSize-1)*betasumdif;//tem;
		difsum=difsum-(double)width*(gammaln((double) suminclust+betasum) + (double)(motifsSize-1)*betasumdif);
		return(difsum);
	}//end of else
}//end of finalRatioSingle

double bmc(const int initialclusternum, const short numChains, const short numItera,//const short bin, 
//+		   const short W, const short WMIN, const short WMAX, 
		   const short MAX_SHIFT, const double factor_0, list <clustClass> &finalClusters, vector <motifClass> &finalallMotifs)
{
	short j,k;
	double like,large;
	list <clustClass> savedClusters;
	vector <motifClass> savedallMotifs;

//	initialCluster(initialclusternum);//+++++
	cout << "start iteration"<<endl;
	for(j=0;j<numChains;j++)
	{
		allClusts.clear();
		motifHome.clear(); 
//added 04/17/07
		for(k=0;k<allMotifs.size();k++)
			allMotifs[k].motifshift=0;
//added 04/17/07
		initialCluster(initialclusternum /*- numChains/2*/ +j);
		cout << "chain" << j+1 << endl;
		like=chain(numItera, /*W, WMIN,WMAX, */ MAX_SHIFT, factor_0,savedClusters,savedallMotifs);//bin,
		cout << "log likelihood =" << like <<", cluster= "<<savedClusters.size()<< allClusts.size()<<endl;
		if ((j==0) || (like > large))
		{
			large=like;
			finalClusters=savedClusters;//allClusts;
                        finalallMotifs=savedallMotifs;			
		}//end of if
	}//end of j
	return(large);
}//end of bmc

double chain(const short NUM_CYCLE, 
//			 const short W, const short WMIN, const short WMAX, 
			const short MAX_SHIFT, const double FACTOR_0, list <clustClass> &savedClusters, vector <motifClass> &savedallMotifs)
//const short bin,
{
	short j;
	short cycle;
//	short stage = 0;
	double like,difference;//,like1;
	list <clustClass>::iterator allClustsIter;
	list <motifClass>::iterator motifsIDIter;
	short currentClusterOrder;
//+	short del;
//+	vector <short> delVector;
	double diff =0;
	double large;

//*****
//  initialize
//*****
	cout << allClusts.size()<<endl;
//*****
//  get initial likelihood
//*****
	like=loglikelihoodSingle();//(0)
//	double sum =0;
//	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//		sum = sum + (*allClustsIter).loglikelihood(0);
	//-cout << "original like=" << like <<endl;
//*****
//  start iteration
//*****
//	stage = 0;///whether in annealing step
	for(cycle = 0; cycle < NUM_CYCLE+ANNEALING;cycle ++ )
	{
		cout <<"cycle="<<cycle<<endl;
//save		like=loglikelihoodSingle();//0);
//save	double sum =0;
//save	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//save	sum = sum + (*allClustsIter).loglikelihood();//0);

//219		cout << "iteration=" << cycle+1 <<"\r";//"cluster=" << allClusts.size() << endl;
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
//save			cout <<"motif= "<<j<<" cluster size = "<<allClusts.size()<<endl;
//218			cout << "motif= " << j <<" " << "name=" <<allMotifs[j].name /*<< "clust=" << allMotifs[j].cluster*/ <<endl;
//			cout <<"shift= "<<allMotifs[23].motifshift<<endl;
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//			{
//				cout <<"width= "<<(*allClustsIter).width<<endl;
//				cout <<"remain= "<<(*allClustsIter).remain<<endl;
//218			}
			difference=operation(j, FACTOR_0,MAX_SHIFT);//del);//,stage);
//+			cout <<difference<<endl;
//save ***
/**/
//save			like1=loglikelihoodSingle();//stage);
			like=like+difference;
//+			cout <<"like="<<like<<endl;
//save			cout << "like=" <<like << "like1=" <<like1<<endl;
/*
			if (fabs(like-like1)>0.00001)
			{
				cout << "problem here: motif = "<< j<<" dif=" << difference << "new like=" << like << " like1=" <<like1<<endl;
				exit(18);
			}
*/
/**/
//save ***			
//***** special check on clusters
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//				(*allClustsIter).printClust(16);
//***** end special check
			if(cycle >= NUM_CYCLE/2)
			{
				if(((cycle == NUM_CYCLE/2)&&(j==0))||(like>large))
				{
					large = like;
					savedClusters =allClusts;
					savedallMotifs = allMotifs;					
//+					cout <<"large="<<large<<" cluster size="<<savedClusters.size()<<endl;
				}
			}
		}
//
// fragmentation step
//			
		//-like1=loglikelihoodEven(clusters,W,stage);
//		cout << "like=" <<like << "like1=" <<like1;
/* 02/19/06
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
02/19/06 */
		//-like1=loglikelihoodEven(clusters,W, stage);
		//-cout << "like=" <<like << "like1=" <<like1;
	}//end of cycle
	//like1=loglikelihoodEven(clusters,WMIN);
//save	like1=loglikelihoodSingle();//0);
//save	cout << "like=" << like+difference << " like1=" <<like1<<endl;
	return(large);
	//return(like+difference);	
}//end of chain

double operation(const short motifOrder, const double FACTOR_0, const short MAX_SHIFT)//const short bin,
				 //short &del, const short stage)
{
	short j;//,k;
	short clustOrder = -1;
	vector <double> ratioVec;
	vector < vector <double> > ratioAll;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustpoint;
	short clusterNumber = 0;
	double ratio;
	short decision = -1, currentCluster = -1,shift = 0;
	clustClass tempClust;
	double likeminus,likeadd;
	bool empty = false;
	double postprob;
	short motifStart, clustStart;
	short motifWidth, clustWidth,newWidth;
//+	double max;
	short newwidth = 0, oldiwdth = 0;
	vector <bool> clustlong;
	short oldmotifshift, oldclustshift;
	short newmotifshift, newclustshift;
//+	short oldclustwidth,newclustwidth;
//	double likeold, likenew;
	vector <list <clustClass>::iterator> clustIterVec; 
	int count;
	bool unsolved = false;

//	the pointer pointed to the cluster that the motif belongs to
	clustpoint = motifHome[motifOrder];
//+	cout <<"shift = "<<allMotifs[144].motifshift<<endl;
//+	int oldclustwidth = (*clustpoint).width;
/*
	short a = (*clustpoint).width;
	short b = (*clustpoint).motifsID.size();
	short k;
	for(j=0;j<4;j++)
	{
		for(k=0;k<(*clustpoint).width;k++)
			cout << (*clustpoint).profile[j][k]<<" ";
		cout <<endl;
	}
	a = allClusts.size();
*/

//  translate clustpoint to clustOrder 

/*
	currentCluster = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if (allClustsIter == clustpoint)
			break;
		else
			currentCluster ++;
	}//end of allClustsIter
*/

// remove motif from cluster 
//	short a = (*clustpoint).motifsID.size();
//	short b = *((*clustpoint).motifsID.begin());//1/12/07
//	short aa = (*clustpoint).width;//1/12/07
//	short bb = allMotifs[motifOrder].width;//1/12/07
//+	if(motifOrder ==136)
//+		cout <<"here.";
//+	cout <<allClusts.size()<<" ";
//+	list <short>::iterator motifsIDIter;
//+	cout <<(*clustpoint).motifsID.size()<<endl;
	if((*clustpoint).motifsID.size() == 1) //for cluster with only 1 member, remove one exisitng cluster
	{
		empty = true;	
//+		likeminus = - (*clustpoint).loglikelihood();//0);
//+		allClusts.erase(clustpoint);
	}//end of if
	else
	{
//+		likeminus=clustRemoveMotif(motifOrder, clustpoint);//, stage);
		empty = false;
	}
//start trying to fit all clusters
	motifWidth = allMotifs[motifOrder].width;//move outside of loop

	count =0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
// decide shift
		clustWidth = (*allClustsIter).width;
/*		
		if((clustWidth -motifWidth)>MAX_SHIFT)//2)//cluster is longer
		{
			ratioVec.push_back(0);
			clustlong.push_back(true);
		}
		else if((motifWidth -clustWidth)>MAX_SHIFT)//2)//motif is longer
		{
			ratioVec.push_back(0);
			clustlong.push_back(false);
		}
		else
*/		
		if(abs(motifWidth -clustWidth) <= MAX_SHIFT)
		{
			list <short>::iterator motifsIDIter;
//+			for(motifsIDIter = (*allClustsIter).motifsID.begin();motifsIDIter != (*allClustsIter).motifsID.end();motifsIDIter ++)
//+				cout <<" "<<(*motifsIDIter) <<" "<<allMotifs[(*motifsIDIter)].motifshift;
//+			cout <<endl;
			if(allClustsIter == clustpoint)
				currentCluster = count;
			count ++;
			
			if(clustWidth >=motifWidth)//cluster is longer
			{
				clustlong.push_back(true);
				if((allClustsIter == clustpoint)&&((*clustpoint).motifsID.size()==2))
				{
					list <short>::iterator motifsIDIter;
					motifsIDIter = (*clustpoint).motifsID.begin();
					motifsIDIter ++;
					newWidth =  allMotifs[(*motifsIDIter)].width;
					shift = newWidth - motifWidth +1;
					for(j=0;j<shift;j++)
					{
						ratio=exp(bayesSpecialRatioSingle(motifOrder, motifsIDIter,0,j));//,motifStart,clustStart,stage));
						ratioVec.push_back(ratio);
					}//end of j
				}
				else
				{
					shift = clustWidth - motifWidth +1;
//+					max = 0;
					for(j=0;j<shift;j++)
					{
						motifStart = 0;
						clustStart = j;
						if(allClustsIter == clustpoint)
						{
							if(!empty)		
							{
								ratio=exp(bayesMinusRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,motifStart,clustStart,stage));
								ratioVec.push_back(ratio);	
							}
						}
						else
						{
							ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,stage));
							ratioVec.push_back(ratio);
						}
					}//end of j
				}//end of else
			}//end of if
			else//motif is longer
			{
				clustlong.push_back(false);
				if((allClustsIter == clustpoint)&&((*clustpoint).motifsID.size()==2))
				{
					list <short>::iterator motifsIDIter;
					motifsIDIter = (*clustpoint).motifsID.begin();
					motifsIDIter ++;
					newWidth =  allMotifs[(*motifsIDIter)].width;
//					shift = motifWidth - newWidth +1;
					if(motifWidth < newWidth)
					{
//						cout <<"error in shift. operation -> motif is longer ->..."<<endl;
//						exit(99);
						shift = newWidth - motifWidth +1;
						for(j=0;j<shift;j++)
						{
							ratio=exp(bayesSpecialRatioSingle(motifOrder, motifsIDIter,0,j));//,motifStart,clustStart,stage));
							ratioVec.push_back(ratio);
						}//end of j
					}
					else
					{
						shift = motifWidth - newWidth +1;
						for(j=0;j<shift;j++)
						{
							ratio=exp(bayesSpecialRatioSingle(motifOrder, motifsIDIter,j,0));//,motifStart,clustStart,stage));
							ratioVec.push_back(ratio);
						}//end of j
					}//end of else
				}//end of if
				else
				{
					shift = motifWidth - clustWidth +1;
					for(j=0;j<shift;j++)
					{
						motifStart = j;
						clustStart = 0;
						if(allClustsIter == clustpoint)
						{
							if(!empty)		
							{
								ratio=exp(bayesMinusRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,stage));
//+								cout <<ratio<<endl;
								ratioVec.push_back(ratio);
							}
						}
						else
						{
							ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,stage));
//+							cout <<ratio<<endl;
							ratioVec.push_back(ratio);
						}	
					}//end of j
				}//end of else
			}//end of else
//+		}//end of if
//+			if(ratioVec.size()>0)
//+			{
			ratioAll.push_back(ratioVec);
/*				
			for(j=0;j<ratioVec.size();j++)
			{
				cout <<ratioVec[j]<<" ";
			}
			cout <<endl;
*/	
			clustIterVec.push_back(allClustsIter);
			ratioVec.clear();
//+			}
		}//end of if
	}//end of allClustIter
// decide which cluster fit the motif best
	decision=shiftdecide(ratioAll, FACTOR_0,shift,postprob);
//+	cout <<"motif= "<<motifOrder<<" decision="<<decision<<" posterior proba = "<<postprob<<endl;
	allMotifs[motifOrder].postprob = postprob;
//+	int a = clustlong.size();
//+	int b = clustIterVec.size();
//+	if(a!=b)
//+		cout <<"error.";
// if the current cluster is the best, then do nothing
	if(empty)
	{
		if(decision ==0)
		{
			clustIterVec.clear();
			return(0);
		}
		else
		{
			if(clustlong[decision-1])//currentCluster
			{
				allMotifs[motifOrder].motifshift = 0;
				allMotifs[motifOrder].clustshift = shift;
			}
			else
			{
				allMotifs[motifOrder].motifshift = shift;
				allMotifs[motifOrder].clustshift = 0;
			}
			clustlong.clear();
			likeminus = -(*clustpoint).loglikelihood();
			likeadd=clustInsertMotif(motifOrder,clustIterVec[decision-1]);//[decision-1]
			allClusts.erase(clustpoint);
//			cout <<"cluster size= "<<allClusts.size()<<endl;
			clustIterVec.clear();
			return(likeminus-likeadd);
		}
	}//end of if
	else if((currentCluster>0)&&(decision==(currentCluster + 1)))//added !empty 08/07/05
	{// same cluster
		oldmotifshift = allMotifs[motifOrder].motifshift;
		oldclustshift = allMotifs[motifOrder].clustshift;
		unsolved = false;
		if((*clustpoint).motifsID.size()>2)
		{
			if(clustlong[currentCluster])
			{
//+				allMotifs[motifOrder].motifshift = 0;
//+				allMotifs[motifOrder].clustshift = shift;
				newmotifshift = 0;
				newclustshift = shift;
			}
			else
			{
//+				allMotifs[motifOrder].motifshift = shift;
//+				allMotifs[motifOrder].clustshift = 0;
				newmotifshift = shift;
				newclustshift = 0;
			}
			clustlong.clear();
			if(/*(oldclustwidth==newclustwidth)&&*/(oldmotifshift == newmotifshift)&&(oldclustshift == newclustshift))
			{
				(*clustpoint).motifsID.pop_front();
				(*clustpoint).motifsID.push_back(motifOrder);
				return(0);
			}//end of if
			else
			{
				unsolved = true;
			}
		}//end of if
		else//(*clustpoint).motifsID.size()==2
		{
			if(newWidth<0)
			{
				cout << "error in newWidth."<<endl;
				exit(199);
			}
//			int newclustwidth = (*clustpoint).width;// remove old and new clustwidth, 02/19/07
			if(newWidth >= motifWidth)//clustlong[currentCluster])
			{
//+				allMotifs[motifOrder].motifshift = 0;
//+				allMotifs[motifOrder].clustshift = shift;
				newmotifshift = 0;
				newclustshift = shift;
			}
			else
			{
//+				allMotifs[motifOrder].motifshift = shift;
//+				allMotifs[motifOrder].clustshift = 0;
				newmotifshift = shift;
				newclustshift = 0;
			}
			int oldclustwidth = (*clustpoint).width;
			int newclustwidth = MIN(motifWidth,newWidth);
			if((oldclustwidth==newclustwidth)&&(oldmotifshift == newmotifshift)&&(oldclustshift == newclustshift))
			{
				(*clustpoint).motifsID.pop_front();
				(*clustpoint).motifsID.push_back(motifOrder);
				return(0);
			}//end of if
			else
			{
				unsolved = true;
			}
		}//end of else
		if(unsolved)
		{
			likeminus=clustRemoveMotif(motifOrder, clustpoint);
			allMotifs[motifOrder].motifshift = newmotifshift;
			allMotifs[motifOrder].clustshift = newclustshift;
			likeadd=clustInsertMotif(motifOrder,clustIterVec[decision-1]);//[decision-1]
//218       cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
			clustIterVec.clear();
            return(likeminus - likeadd);
		}//end of else
	}//end of else if
// if a different cluster is desired
	else 
	{
		likeminus=clustRemoveMotif(motifOrder, clustpoint);
/*//removed 02/19/07		
		if(clustlong[decision-1])//currentCluster
		{
			allMotifs[motifOrder].motifshift = 0;
			allMotifs[motifOrder].clustshift = shift;
		}
		else
		{
			allMotifs[motifOrder].motifshift = shift;
			allMotifs[motifOrder].clustshift = 0;
		}
		clustlong.clear();
*/// removed 02/19/07
		if(decision > 0)
		{
			if(clustlong[decision-1])//currentCluster
			{
				allMotifs[motifOrder].motifshift = 0;
				allMotifs[motifOrder].clustshift = shift;
			}
			else
			{
				allMotifs[motifOrder].motifshift = shift;
				allMotifs[motifOrder].clustshift = 0;
			}
			clustlong.clear();
//+			likeadd=clustInsertMotif(motifOrder,decision-1);//, stage);// moved down from before if 08/08/05
			likeadd=clustInsertMotif(motifOrder,clustIterVec[decision-1]);//[decision-1]
		}//end of if
// if a new cluster is desired
		else if (decision == 0)
		{
			allMotifs[motifOrder].motifshift = 0;
			allMotifs[motifOrder].clustshift = 0;
			tempClust.motifsID.push_back(motifOrder);
			//tempClust.width = allMotifs[motifOrder].width;
			tempClust.iniProfile();//bin);
			allClusts.push_back(tempClust);
			tempClust.motifsID.clear();
			clustpoint = allClusts.end();
			clustpoint --;
			motifHome[motifOrder] =clustpoint;
//+			likeadd = allMotifs[motifOrder].ratio;
			likeadd = allRatio[0][motifOrder];
		}//end of else if		
//		cout <<"cluster size= "<<allClusts.size()<<endl;
		clustIterVec.clear();
		return(likeminus-likeadd);
	}//end of else
}//end of operation

/*
double operation(const short motifOrder, const short bin, const double FACTOR_0, const short MAX_SHIFT)
				 //short &del, const short stage)
{
	short j;//,k;
	short clustOrder = -1;
	vector <double> ratioVec;
	vector < vector <double> > ratioAll;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustpoint;
	short clusterNumber = 0;
	double ratio;
	short decision = -1, currentCluster = -1,shift = 0;
	//short start;
	//motifClass temMotif;
	//list <motifClass> tempMotifs;
	clustClass tempClust;
	double likeminus,likeadd;
	bool empty = false;
	double postprob;
	short motifStart, clustStart;
	short motifWidth, clustWidth;
	double max;
	short newwidth = 0, oldiwdth = 0;
	vector <bool> clustlong;
	short oldmotifshift, oldclustshift;
	short oldclustwidth,newclustwidth;
	double likeold, likenew;

//	the pointer pointed to the cluster that the motif belongs to
	clustpoint = motifHome[motifOrder];
	oldclustwidth = (*clustpoint).width;
// *****
//	short a = (*clustpoint).width;
//	short b = (*clustpoint).motifsID.size();
//	for(j=0;j<4;j++)
//	{
//		for(k=0;k<(*clustpoint).width;k++)
//			cout << (*clustpoint).profile[j][k]<<" ";
//		cout <<endl;
//	}	
// *****
// ********** temp add 02/18/06	
//	short k;
//	for(j=0;j<4;j++)
//	{
//		short wid = allMotifs[motifOrder].width;
//		short shf = allMotifs[motifOrder].motifshift;
//		for(k=0;k<allMotifs[motifOrder].width;k++)
//			cout << allMotifs[motifOrder].profile[j][k+allMotifs[motifOrder].motifshift] << " ";
//		cout <<endl;
//	}
//	cout <<endl;
//
//*********** temp add 02/18/06
//  translate clustpoint to clustOrder 
	currentCluster = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if (allClustsIter == clustpoint)
			break;
		else
			currentCluster ++;
	}//end of allClustsIter
//+	del=-1;
//218	cout << "start of operation, motif = " <<motifOrder<<"size = " <<(*clustpoint).motifsID.size() <<endl;
	// remove motif from cluster 
//	short a = (*clustpoint).motifsID.size();
//	short b = *((*clustpoint).motifsID.begin());//1/12/07
//	short aa = (*clustpoint).width;//1/12/07
//	short bb = allMotifs[motifOrder].width;//1/12/07
	if((*clustpoint).motifsID.size() == 1) //for cluster with only 1 member, remove one exisitng cluster
	{
		empty = true;
		//likeminus=loglikelihoodUpdateSingle(motifOrder, clustpoint,stage);	
		likeminus = - (*clustpoint).loglikelihood();//0);
		allClusts.erase(clustpoint);
//218		cout <<"empty likeminus = "<<likeminus<<endl;
	}//end of if
	else
		likeminus=clustRemoveMotif(motifOrder, clustpoint);//, stage);
//	a = (*clustpoint).width;
//start trying to fit all clusters
	short tmp = allClusts.size();
	short count =0;
	motifWidth = allMotifs[motifOrder].width;//move outside of loop
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		count ++;
//		if(count ==42)
//			cout <<"there"<<endl;
// decide shift
		clustWidth = (*allClustsIter).width;
//moved outside of loop		motifWidth = allMotifs[motifOrder].width;

		if((clustWidth -motifWidth)>MAX_SHIFT)//2)//cluster is longer
		{
			ratioVec.push_back(0);
			clustlong.push_back(true);
		}
		else if((motifWidth -clustWidth)>MAX_SHIFT)//2)//motif is longer
		{
			ratioVec.push_back(0);
			clustlong.push_back(false);
		}
		else if(clustWidth >=motifWidth)//cluster is longer
		{
			clustlong.push_back(true);
			shift = clustWidth - motifWidth +1;
			max = 0;
			for(j=0;j<shift;j++)
			{
				motifStart = 0;
				clustStart = j;
				ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,stage));
//+				cout <<ratio<<endl;
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
				ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,motifStart,clustStart));//,stage));
//+				cout <<ratio<<endl;
				ratioVec.push_back(ratio);
			}//end of j
		}//end of else
		ratioAll.push_back(ratioVec);
		ratioVec.clear();
	}//end of allClustiter
// decide which cluster fit the motif best
	decision=shiftdecide(ratioAll, FACTOR_0,shift,postprob);
//+	cout <<"motif= "<<motifOrder<<" decision="<<decision<<" posterior proba = "<<postprob<<endl;
	allMotifs[motifOrder].postprob = postprob;
//	cout << "decision=" << decision <<endl;
// if the current cluster is the best, then do nothing
	if((!empty)&&(decision==(currentCluster + 1)))//added !empty 08/07/05
	{// same cluster
		//(*clustpoint).motifs.pop_front();
		oldmotifshift = allMotifs[motifOrder].motifshift;
		oldclustshift = allMotifs[motifOrder].clustshift;
		newclustwidth = (*clustpoint).width;
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
		//(*clustpoint).updateProfile(1,allMotifs, motifOrder);
		//(*clustpoint).motifs.push_back(allMotifs[motifOrder]);
		clustlong.clear();
		if((oldclustwidth==newclustwidth)&&(oldmotifshift == allMotifs[motifOrder].motifshift)&&(oldclustshift == allMotifs[motifOrder].clustshift))
		{//used to be (clustWidth ==motifWidth)
			(*clustpoint).addMotif( motifOrder);
//			for(j=0;j<4;j++)
//			{
//				short k;
//				for(k=0;k<oldclustwidth;k++)
//					cout << (*clustpoint).profile[j][k]<<" ";
//				cout <<endl;
//			}	
//
			return(0);
		}
		else 
		{
			likeold = (*clustpoint).loglikelihood();//0);
			(*clustpoint).addMotif( motifOrder);
			likenew = (*clustpoint).loglikelihood();//0);
            likeadd = likeold - likenew;
//218       cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
            return(likeminus - likeadd) ;
		}
	}//end of if
// if a different cluster is desired
	else 
	{
// remove it from original cluster
		//(*clustpoint).motifs.pop_front();
		if(decision > 0)
		{
// insert it shorto the desired cluster
			//likeadd=clustInsertMotif(motifOrder,decision-1,clustWidth);//, stage);
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
			likeadd=clustInsertMotif(motifOrder,decision-1);//, stage);// moved down from before if 08/08/05
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
			clustpoint = allClusts.end();
			clustpoint --;
			motifHome[motifOrder] =clustpoint;
//			tempMotifs.clear();
			//likeadd=loglikelihoodNewEven(clustWidth,motifsIDIter,stage);//***)
			if (empty)
				likeadd=likeminus;
			else
			{
				likeadd=emptyLoglikelihoodupdateSingle(motifOrder,clustpoint);//,stage);
			}//end of else
		//-	cout << " likeadd2=" <<likeadd<<endl;
		}//end of else if		
// remove it from its current cluster
//+		del=motifOrder;//.first
		clustlong.clear();
//218		cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
		return(likeminus-likeadd);
	}//end of else
}//end of operation
*/

/*
short checkWidth(	list <clustClass>::iterator clustpoint)
{
	list <short>::iterator allmotifsIDIter;
	short newiwdth;
	bool begin = true;

	begin = true;
	for(allmotifsIDIter = (*clustpoint).motifsID.begin();allmotifsIDIter != (*clustpoint).motifsID.end(); allmotifsIDIter++)
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
*/

/*
void initWidth(const short bin, const short W)
{
	short k;
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
//	for(j=0;j<numOfCluster;j++)
//	{
//		count=0;
//		for(k=0;k<number_of_motif.total;k++)
//		{
//			if((motif[k].cluster==j)&&((bin==3)||(motif[k].oddeven==bin)))
//			{
//				clust[j].member[count]=k;
//				count++;
//			}//end of if
//		}//end of k
//	}//end of j
}//end of initWidth
*/

/*
double metSingle(const list <clustClass>::iterator clustIter, short indel, short position, short W)
{
	short m;
	double oldProfile, newProfile;
	double sumin, sumex, proba;
	double ran;//,decide,partialsum;
	short suminclust,la;
	int na,nb,nc;
	double dif;
	double tem;
	list <short>::iterator motifsIDIter;
	short shift;

	// * seed random number generator * 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
	
	na=nb=nc=19;
//
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
			la = allMotifProfiles[m][allMotifs[(*motifsIDIter)].start + position + shift];
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
*/

/*
double onoffMetro(const list <clustClass>::iterator clustIter, short W, short WMIN, short bin)
{
	int na,nb,nc;
	double dif,ran,lower;
	short current, leftMargin,rghtMargin;

	// seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
	na=nb=nc=50;
//
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
*/

/*
double fragment(const short bin, const short WMIN, const short WMAX, const short W)
{
	double dif;
	list <clustClass>::iterator clustIter;

	dif=0;
	for(clustIter = (allClusts).begin();clustIter != (allClusts).end();clustIter ++)
	{			
		dif=dif+onoffMetro(clustIter,W,WMIN,bin);
		//short k;
		//for(k=0;k<W;k++)
		//	cout << (*clustIter).power[k];
		//cout <<endl;
	}
	return dif;
}//end of fragmentconst
*/

double clustRemoveMotif(const short motifOrder, const list <clustClass>::iterator clustpoint)//, const short stage)
{
//	short oldwidth = 0;
//	short newwidth = 0;
	double likenew,likeold;

//	oldwidth = (*clustpoint).width;
	likeold = (*clustpoint).loglikelihood();//0);
	(*clustpoint).minusMotif(motifOrder);
	likenew = (*clustpoint).loglikelihood();//0);
	return(likenew-likeold);
	cout <<"remove: old= "<<likeold<<"new= "<<likenew<<endl;
}// end of clustRemoveMotif

double clustInsertMotif(const short motifOrder, const list <clustClass>::iterator clustpoint)//, const short stage)
{
//	short count;
//+	list <clustClass>::iterator clustIter;
	double likenew,likeold;
	short clustWidth = 0;
	short newwidth = 0;
			
	likeold = (*clustpoint).loglikelihood();//0);
	(*clustpoint).addMotif( motifOrder);
	likenew = (*clustpoint).loglikelihood();//0);
// assign this cluster to motifHome
	motifHome[motifOrder] = clustpoint;
//218			cout <<"insert: old= "<<likeold<<"new= "<<likenew<<endl;
	return(likeold-likenew);		
}// end of clustInsertMotif

/*
double clustInsertMotif(const short motifOrder, const short decision)//, const short stage)
{
	short count;
	list <clustClass>::iterator clustIter;
	double likenew,likeold;
	short clustWidth = 0;
	short newwidth = 0;
			
	count=0;
	for(clustIter = allClusts.begin();clustIter != allClusts.end();clustIter ++)
	{
		if(count==decision)
		{
			likeold = (*clustIter).loglikelihood();//0);
			(*clustIter).addMotif( motifOrder);
			likenew = (*clustIter).loglikelihood();//0);
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
*/

short shiftdecide(vector <vector <double> > &ratioAll,const double FACTOR_0, short &shift, double &postprob)
{
	short j,k;
	double sum, insum;
	short number;
	int na,nb,nc;
	double ran;
	short decide = 0;
	double compare, partialsum,accumu;

	//seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
//	na=nb=nc=193;
//
	sum=FACTOR_0;
	number=ratioAll.size();
	for(j=0;j<number;j++)
	{
		for(k=0;k<(ratioAll[j]).size();k++)				
		{
			sum=sum+ratioAll[j][k];
//+			cout <<"sum="<<sum<<endl;
		}//02/18/07
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
short decide(vector <double> &ratioVec, const double FACTOR_0, double &postprob)
{
	short j;
	double sum;
	short number;
	int na,nb,nc;
	double ran;
	short decide = 0;
	double compare, partialsum;

// seed random number generator 
	srand((unsigned)time(NULL));
	na=rand() +1;
	nb=rand() -1;
	nc=rand() ;
	na=nb=nc=50;
//
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

double emptyLoglikelihoodupdateSingle(const short motifOrder, const list <clustClass>::iterator allClustsIter)//, const short stage)
{
//+	short ** newProfile;
	short j,k;
	double logsum;//,la,lsum;
	double logsumnew,lanew, lsumnew;
//+	double motifsum;
	short lsummotif;//,lb;
	double difference;//,tem;
	short begin;
	short width;

	width = (*allClustsIter).width;
//+	newProfile = new short*[4];
//+	for(j=0;j<4;j++)
//+		newProfile[j] = new short[width];
	begin = allMotifs[motifOrder].motifshift;
//+	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	//lsum=0.0;
//+	motifsum=0;
	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
/* +
	for (j=0;j<width;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+begin];
				newProfile[k][j] = allMotifProfiles[k][allMotifs[motifOrder].start + j +begin];//***shift to start//remove (*allClustsIter).profile[k][j]
			}//end of k
//-		}//end of if
	}//end of j
+ */
	for(j=0;j<width;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				logsum=logsum+gammaln(beta[k]);
				lanew = allMotifProfiles[k][allMotifs[motifOrder].start + j +begin];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
//-		}//end of if
/*		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
*/
	}//end of j	
	difference=logsumnew-logsum;//+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
//+	for(j=0;j<4;j++)
//+		delete [] newProfile[j];
//+	delete [] newProfile;
	return (-difference);
}//end of emptyLoglikelihoodupdateSingle

//-----
// calculate likelihood for clusters of single motifs
//-----
double loglikelihoodSingle()//const short stage)
{
	short j,k;
//+	double tem;
	short lsum;//,sumin;
	list <clustClass>::iterator allClustsIter;
//+	list <short>::iterator motifsIDIter;
	double logsum,la;
//	short start;
	short motifsSize;
//+	double logsumold =0;

//	tem=gammaln(betasum)-sumbeta;
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
		//+short a = (*allClustsIter).width;
		for(j=0;j<(*allClustsIter).width;j++)
		{
//-			if((stage!=1)||((*allClustsIter).power[j]))
//-			{
				for(k=0;k<4;k++)
				{
					//+short a=(*allClustsIter).profile[k][j];
					la = (double) (*allClustsIter).profile[k][j];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
//+				logsum=logsum-gammaln((double) lsum + betasum);
//-			}//end of if
/*			else if ((stage==1)&&(!(*allClustsIter).power[j]))
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
*/
		}//end of j
		logsum = logsum + (double) (*allClustsIter).width*(betasumdif -gammaln((double) lsum + betasum));//*tem;
//218		cout <<"cluster like="<<logsum-logsumold<<endl;
//+        logsumold = logsum;
	}//end of allClustsIter
//218	cout <<endl;
	//logsum=logsum+allClusts.size()*(double) width *tem;
	return(logsum);
}//end of loglikelihoodSingle

double bayesRatioSingle(short motifOrder, list <clustClass>::iterator allClustsIter,short motifStart, short clustStart)//, const short stage)
{
	short j,k;
	short Lold,Lnew,sum;//,Lz;
	double sumbayes,tem;
	double newProfile;
	double bayesRatio;
	short clustProfile,motifProfile;
	short width,place, widthdif;
	
	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);//changed 9/13/05
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum+(*allClustsIter).profile[j][0];
	}
	Lold=sum;
//	Lz=allMotifs[motifOrder].depth;
	Lnew=Lold + allMotifs[motifOrder].depth;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
//-		if((stage !=1)||((*allClustsIter).power[j]))//removed ==1 08/10/05
//-		{
			for(k=0;k<4;k++)
			{
				clustProfile = (*allClustsIter).profile[k][j+clustStart];
				motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
				if(clustProfile<0)
					cout << "clust ="<<k<<" "<<j<<" "<<clustStart<<endl;
				if(motifProfile<0)
					cout << "clust = "<<k<<" "<<j<<" "<<motifStart<<endl;
				newProfile = clustProfile + motifProfile;
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) clustProfile + beta[k]);
//+					-gammaln((double) motifProfile + beta[k]);
		//		if((motifOrder ==136)&&(motifStart==0))
		//			cout <<clustProfile<<" "<<motifProfile<< " ";
			}//end of k
		//	if((motifOrder ==136)&&(motifStart==0))
		//		cout <<endl;
//-		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum);
//+		+gammaln((double) Lz+betasum) -gammaln(betasum)+sumbeta;
	if(allMotifs[motifOrder].width <= width)
		place = 0;
	else
	{
		widthdif = allMotifs[motifOrder].width - width;
		place = widthdif*(widthdif+1)/2 + motifStart;
	}
///	if(stage == 1)
///		bayesRatio = sumbayes+tem*(double) (*allClustsIter).remain;
///	else 
/*
	if(place >0)
	{
		for(j=0;j<6;j++)
			cout << allRatio[j][motifOrder]<<"  ";
	}
*/
	bayesRatio = sumbayes+tem*(double) width + allRatio[place][motifOrder];
	return(bayesRatio);
}//end of bayesRatioSingle

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,short motifStart, short clustStart)
{
	int j,k;
	int Lold,Lnew,sum;
	double sumbayes,tem;
	int newProfile;
	int difProfile,shift;
	double bayesRatio;
	short clustProfile,motifProfile;
	short width,place, widthdif;
	
	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum + (*allClustsIter).profile[j][0];
	}
	Lnew=sum;
//	Lz=allMotifs[motifOrder].depth;
	Lold=Lnew - allMotifs[motifOrder].depth;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
		for(k=0;k<4;k++)
		{
//+			clustProfile = (*allClustsIter).profile[k][j+clustStart];
//+			motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
//+			oldProfile = clustProfile - motifProfile;
			//oldProfile=(*allClustsIter).profile[k][j] - allMotifProfiles[k][motifOrder*fixedwidth + j];//[(j) + clustStart]
				//allMotifs[motifOrder].profile[k][j+motifStart];
			motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
			shift = allMotifs[motifOrder].motifshift;
			difProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + shift];
//+			clustProfile = allClustProfiles[k][(*allClustsIter).start + j+clustStart] - difProfile;
			clustProfile = (*allClustsIter).profile[k][j+clustStart] - difProfile;
			newProfile = clustProfile + motifProfile;
			
			sumbayes=sumbayes+gammaln((double) newProfile + beta[k])
				-gammaln((double) clustProfile + beta[k]);//[(j) + clustStart]
			//	-gammaln((double) allMotifProfiles[k][motifOrder*fixedwidth + j] + beta[k]);
		}//end of k
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum);	
	if(allMotifs[motifOrder].width <= width)
		place = 0;
	else
	{
		widthdif = allMotifs[motifOrder].width - width;
		place = widthdif*(widthdif+1)/2 + motifStart;
	}
	//double x = allMotifs[motifOrder].ratio;
	bayesRatio = sumbayes+tem*(double) width + allRatio[place][motifOrder];
	return(bayesRatio);
}//end of bayesMinusRatioSingle

double bayesSpecialRatioSingle(int motifOrder, list <short>::iterator motifsIDIter, short motifStart, short clustStart)
{//when clustsize =2, motfwidth <remaining single motif width. 
	int j,k;
	int Lold,Lnew;
	double sumbayes,tem;
	int newProfile;
	double bayesRatio;
	short clustProfile,motifProfile;
	short width,place, widthdif;
	short motifWidth, otherWidth;
	
	motifWidth = allMotifs[motifOrder].width;
	otherWidth = allMotifs[(*motifsIDIter)].width;
	width = MIN(motifWidth, otherWidth);
	Lold = allMotifs[(*motifsIDIter)].depth;
	Lnew = allMotifs[motifOrder].depth + allMotifs[(*motifsIDIter)].depth;
//	Lz=allMotifs[motifOrder].depth;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
		for(k=0;k<4;k++)
		{
			motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
			clustProfile = allMotifProfiles[k][allMotifs[(*motifsIDIter)].start + j + clustStart];
			newProfile = clustProfile + motifProfile;
			
			sumbayes=sumbayes+gammaln((double) newProfile + beta[k])
				-gammaln((double) clustProfile + beta[k]);//[(j) + clustStart]
			//	-gammaln((double) allMotifProfiles[k][motifOrder*fixedwidth + j] + beta[k]);
		}//end of k
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum);	
	if(otherWidth >= motifWidth)
		place = 0;
	else
	{
		widthdif = motifWidth - otherWidth;
		place = widthdif*(widthdif+1)/2 + motifStart;
	}
	bayesRatio = sumbayes+tem*(double) width + allRatio[place][motifOrder];
	return(bayesRatio);
}//end of bayesSpecialRatioSingle

double loglikelihoodUpdateSingle (const short motifOrder, const list <clustClass>::iterator allClustsIter)//, const short stage)
{
//+	short ** newProfile;
	short j,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
//+	double motifsum;
	short lsummotif;;//,lb;
	double difference;//,tem;
	short width;
	short begin;
	
	width = (*allClustsIter).width;
	// need to find out if the cluster width changed during the add/removal of the motif
//+	newProfile = new short*[4];
//+	for(j=0;j<4;j++)
//+		newProfile[j] = new short[width];
	begin = allMotifs[motifOrder].motifshift;
//+	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
//+	motifsum=0;
	for (j=0;j<4;j++)
	{
		lsum=lsum+(*allClustsIter).profile[j][0];		   		
	}
	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=lsum+lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
/* +
	for (j=0;j<width;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+begin];
				newProfile[k][j]=(*allClustsIter).profile[k][j] + allMotifProfiles[k][allMotifs[motifOrder].start + j + begin];
//							 +allMotifs[motifOrder].profile[k][j+begin];//***shift to start
			}//end of k
//-		}//end of if
	}//end of j
+ */
	for(j=0;j<width;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				la=(*allClustsIter).profile[k][j];
				logsum=logsum+gammaln((double)la+beta[k]);
				//+ lanew=newProfile[k][j];
				lanew = la + allMotifProfiles[k][allMotifs[motifOrder].start + j + begin];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
//-		}//end of if
/*		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
*/
	}//end of j	
	difference=logsumnew-logsum;//+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
//+	for(j=0;j<4;j++)
//+		delete [] newProfile[j];
//+	delete [] newProfile;
	return difference;
}//end of loglikelihoodupdateSingle

//----- read in motif block data and put shorto map and list -----//
short readSeq(const char inputFileName[], const short MAX_SHIFT, short &widthmin)
{
	int j,k;
	list <motifClass> clust;
	motifClass temMotif;
	vector <string> seq;
	ifstream inFile,inFile2;
	string lineString;
	istringstream iss;
	string motif;
	char c;
	int rowCount;
	int motifCount = 0;
	int count;
	string tempName;
	ofstream outputFile;
	int signal=0;//signal for if it is the first line of the sfirst motif
	int minwidth = 0;
	int maxwidth = 0;
	int totalcol;
	int start, width;//, order;
	int shifttotal;

// get names, etc
	cout << "start reading in data" << endl;
	inFile.open(inputFileName);//.c_str());
	if(!inFile)
	{
		cout << "ERROR: Unable to open file: " << inputFileName << endl;
	    exit(2);
	}//end of if
// start reading in data
	motifCount=0;
	signal=0;
	totalcol =0;
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
			temMotif.motifshift = 0;
			temMotif.width=motif.length();
			temMotif.depth=rowCount;
			allMotifs.push_back(temMotif);
			motifCount++;
			int a = temMotif.width;
			totalcol = totalcol + temMotif.width;
		}
		else if(c=='>') // first line of FASTA format, take name
		{
			iss >> tempName;
			if((signal==0)||(tempName!=temMotif.name))
			{
				signal=1;
//				temMotif.id=motifCount;
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
//			seq.push_back(motif);
			rowCount++;
		}
	}//end of while
	inFile.close();

	allMotifProfiles = new int*[4];
	for(j=0;j<4;j++)
	{
		allMotifProfiles[j] = new int[totalcol];
		for(k=0;k<totalcol;k++)
			allMotifProfiles[j][k]=0;
	}

	shifttotal = (MAX_SHIFT+1)*(MAX_SHIFT+2)/2;
	allRatio = new double*[shifttotal];
	for(j=0;j<shifttotal;j++)
	{
		allRatio[j] = new double[motifCount];
		for(k=0;k<motifCount;k++)
			allRatio[j][k]=0;
	}

// get profiles:
	inFile2.open(inputFileName);//.c_str());
	if(!inFile2)
	{
		cout << "ERROR 2: Unable to open file: " << inputFileName << endl;
	    exit(2);
	}//end of if
// start reading in data
	signal=0;
	start=0;
	count =0;
	while (inFile2)
	{
// get name
		getline(inFile2,lineString);
// get DNA sequence
		iss.clear();
		iss.str(lineString + " ");
		iss >> c;
		if(!iss) // at the end of a motif just read in, summerize information
		{
//			temMotif.motifshift = 0;
//			temMotif.width=motif.length();
//			temMotif.depth=rowCount;
			width = motif.length();
			genProfile(seq,start,count, MAX_SHIFT);
			seq.clear();
			allMotifs[count].start = start;
			start = start + width;
//			allMotifs.push_back(temMotif);
			count++;
		}
		else if(c=='>') // first line of FASTA format, take name
		{
			iss >> tempName;
			if((signal==0)||(tempName!=temMotif.name))
			{
				signal=1;
//				temMotif.name=tempName;
//				rowCount=0;
				iss.ignore();
			}
		}
		else // subsequent lines of FASTA files, take sequences
		{
			iss.clear();
			iss.str(lineString + " ");
			iss >> motif;
			seq.push_back(motif);
//			rowCount++;
		}
	}//end of while
	inFile2.close();
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
	minwidth = allMotifs[0].width;
	maxwidth = minwidth;
	for(j=0;j<allMotifs.size();j++)
	{
		if((j>0)&&(allMotifs[j].width > maxwidth))
		{
			maxwidth = allMotifs[j].width;
		}
		else if((j>0)&&(allMotifs[j].width < minwidth))
		{
//-			cout << minwidth<<endl;
			minwidth = allMotifs[j].width;
		}

	}
	cout << "total sites" << allMotifs.size()<<endl;
	cout <<"min width = "<< minwidth<<endl;
	cout <<"max width = "<< maxwidth<<endl;
	//exit(0);
	widthmin = minwidth;
	return(motifCount);
}//end of readSeq

void genProfile(vector <string> seq, const int start, const int motifOrder, const short MAX_SHIFT)
{
	int j,k,m,n;
	int depth, width;
	double sum, tem, ratio;
	short count;

	depth = allMotifs[motifOrder].depth;
	width = allMotifs[motifOrder].width;
// get count
	for(j=0;j<depth;j++)
	{
		for(k=0;k<width;k++)
		{
			switch (seq[j][k])
			{
			case 'a': allMotifProfiles[0][start + k]++; break;
			case 'c': allMotifProfiles[1][start + k]++; break;
			case 'g': allMotifProfiles[2][start + k]++; break;
			case 't': allMotifProfiles[3][start + k]++; break;
			case 'A': allMotifProfiles[0][start + k]++; break;
			case 'C': allMotifProfiles[1][start + k]++; break;
			case 'G': allMotifProfiles[2][start + k]++; break;
			case 'T': allMotifProfiles[3][start + k]++; break;
			default: cout << "error in motif profile" <<endl;
			}//end of switch
		}//end of k
	}//end of j
/*	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << allMotifProfiles[j][k] << " ";
		cout << endl;
	}//end of j
	cout << endl;
*/
	tem = gammaln((double) depth +betasum) - betasumdif;
	count =0;
	for(j=0;j<(MAX_SHIFT+1);j++)
	{
		for(k=0;k<(j+1);k++)
		{
			sum=0;
			for(m=0;m<(width-j);m++)
			{
				for(n=0;n<4;n++)					
					sum = sum + gammaln((double) allMotifProfiles[n][start + m + k] + beta[n]);
			}//end of m
			ratio = -sum + tem*(double)(width-j);
			allRatio[count][motifOrder] = ratio;
			count ++;
		}//end of k
	}//end of j
/*
	sum = 0;
	for(j=0;j<width;j++)
	{
		for(k=0;k<4;k++)
			sum = sum + gammaln((double) allMotifProfiles[k][start + j] + beta[k]);
	}
	tem = gammaln((double) depth +betasum) - betasumdif;//-gammaln(betasum)+sumbeta;
	ratio = -sum + tem*(double)width;
	allMotifs[motifOrder].ratio = ratio;
*/
}//end of genProfile

void initialCluster(const short INITIAL_SIZE)//const short motifCount, const short bin,
{
	short j;
	list <motifClass> clust;
	list < clustClass >::iterator allClustsIter;
	list <short>::iterator motifsIDIter;
//	list < motifClass >::iterator motifsIDIter;
	clustClass tempClust;
	short motifsSize;
	short oldwidth, order, count;

/*
	for(j=0;j<motifCount;j++)
	{
		if((bin ==2)||((bin ==0)&&(fullMotifVector[j].isEven))||((bin ==1)&&(!fullMotifVector[j].isEven)))
			allMotifs.push_back(fullMotifVector[j]);
	}//end of j
*/
	motifsSize = allMotifs.size();
//	width = new short[INITIAL_CLUSTER];
//	for(j=0;j<INITIAL_CLUSTER;j++)// initialize clusters
//		width[j] = allMotifs[j].width;
//	for(j=INITIAL_CLUSTER;j<motifsSize;j++)
//	{
//		group = j%INITIAL_CLUSTER;
//		if(allMotifs[j].width < width[group])
//			width[group] = allMotifs[j].width;
//	}//end of j
// ***** added 02/22/06	
	count =0;
	oldwidth = 0;
	for(j=0;j<motifsSize;j++)
	{
//-		cout <<allMotifs[j].name<<endl;
		allMotifs[j].clustshift =0;
		if((oldwidth != allMotifs[j].width)||(count == (INITIAL_SIZE -1)))
		{// need to start a new cluster
			if(j>0)//wrap up old cluster
			{
				allClusts.push_back(tempClust);
				tempClust.motifsID.clear();
				count =0;
			}
			tempClust.motifsID.push_back(j);
			oldwidth = allMotifs[j].width;
			tempClust.width = oldwidth;
			tempClust.iniProfile();//bin);
		}
		else// continue to add motifs to the current clusters
		{
			count ++;
			tempClust.addMotif(j);
		}
	}
	allClusts.push_back(tempClust);
	tempClust.motifsID.clear();
// ***** set motifHome
	for(j=0;j<motifsSize;j++)
		motifHome.push_back(allClusts.begin());
//		motifHome.push_back(NULL);
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		for(motifsIDIter = (*allClustsIter).motifsID.begin();motifsIDIter != (*allClustsIter).motifsID.end(); motifsIDIter++)
		{
			order = *motifsIDIter;
			motifHome[order] = allClustsIter;
		}
	}
// ***** added 02/22/06
// ***** deleted 02/22/06
/*
	for(j=0;j<INITIAL_CLUSTER;j++)// initialize clusters
	{
		allMotifs[j].clustshift =0;
		tempClust.motifsID.push_back(j);
//		tempClust.width = width[j];
		tempClust.width = allMotifs[j].width;
		tempClust.iniProfile(bin);
		short b = tempClust.profile[0][0];
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
			short a = (*allClustsIter).profile[0][0];
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
*/
// ***** deleted 02/22/06
//	delete []width;
//222	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//222		(*allClustsIter).printClust();
	//exit(0);
}//end of initialCluster

/*
void initialCluster(const short bin, const short motifCount, const vector <motifClass> fullMotifVector,const short INITIAL_CLUSTER)
{
	short j;
	list <motifClass> clust;
	list < clustClass >::iterator allClustsIter;
	list < motifClass >::iterator motifsIDIter;
	clustClass tempClust;
	short motifsSize;
//	short *width, group;

	for(j=0;j<motifCount;j++)
	{
		if((bin ==2)||((bin ==0)&&(fullMotifVector[j].isEven))||((bin ==1)&&(!fullMotifVector[j].isEven)))
			allMotifs.push_back(fullMotifVector[j]);
	}//end of j
	motifsSize = allMotifs.size();
//	width = new short[INITIAL_CLUSTER];
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
		short b = tempClust.profile[0][0];
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
			short a = (*allClustsIter).profile[0][0];
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
//219	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//219		(*allClustsIter).printClust();
}//end of initialCluster
*/

double gammaln(double xx)
{
	double ser,stp,tmp,x,y,cof[6],gam;
	short j;
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
