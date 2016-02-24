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
#define ANNEALING 0
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
static double betasumhalf=0.5*betasum;
static double sumbetahalf=gammaln(0.5*beta[0])+gammaln(0.5*beta[1])+gammaln(0.5*beta[2])+
				gammaln(0.5*beta[3]);
static int kstar[4]={3,2,1,0};//{1,0,3,2};

class motifClass
{
private:
public:
	int start;
	string name;
	int width;
	int depth;
	int clust;
	int motifshift;
	int clustshift;
	double postprob;
};

vector <motifClass> allMotifs;

int **allMotifProfiles;

int **allClustProfiles;

class clustClass
{
private:
public:
	list <int> motifsID;
	int start;
	int width;
//	vector <int> profile[4];
	void printClust();
//	void clustClass::iniProfile();
//+	void clustClass::addMotif(/*vector <motifClass> allMotifs,*/ int motifOrder);
//+	int clustClass::minusMotif(/*vector <motifClass> allMotifs,*/ int motifOrder);
//+	void clustClass::updateProfile(int sign, /*vector <motifClass> allMotifs,*/ int motifOrder);
//	void clustClass::getProfile();
	double loglikelihood();//const int stage);
	~clustClass()
	{
		motifsID.clear();
//		profile[0].clear();
//		profile[1].clear();
//		profile[2].clear();
//		profile[3].clear();
	}
};

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
		cout << "ID " << count;
		cout << "Name " << allMotifs[(*motifsIDIter)].name << " " <<allMotifs[(*motifsIDIter)].width << 
			" " << allMotifs[(*motifsIDIter)].depth << " " << /*allMotifs[(*motifsIDIter)].isEven <<*/ endl;
		count++;
	}//end of motifsIDIter

	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
		{
//			outputFile << allClustProfiles[j][k + start] <<" ";
			cout << allClustProfiles[j][k + start] << " ";
		}
		cout << endl;
//		outputFile << endl;
	}//end of j
//	outputFile << endl;
//	outputFile.close();
}//end of clustClass::printClust()

//----- initialize the profile matrix with the entry of the first motif.
/*
void clustClass::iniProfile()
{
	int j,k;
	list <int>::iterator motifsIDIter;
//	int motifWidth;

///	left=0;
///	right=0;
	
	motifsIDIter=motifsID.begin();
	//int a = (*motifsIDIter).start;
	//int b = (*motifsIDIter).shift; 
//	start=(*motifsIDIter).start + (*motifsIDIter).shift;//added (*motifsIDIter).shift 
	width=allMotifs[(*motifsIDIter)].width;
//- ******* moved from before motifsIDIter
//-
//-	if(bin==2)//added 08/17/05
//-		remain = width;//added 08/17/05
//-	else//added 08/17/05
//-		remain=(double) width*0.5;
//-
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
	//		int b = (*motifsIDIter);
	//		int a = allMotifs[(*motifsIDIter)].profile[k][j];
			profile[k][j]=allMotifProfiles[k][allMotifs[(*motifsIDIter)].start + j];//+(*motifsIDIter).profile[kstar[k]][motifWidth-1-j-start];
		}
	}//end of j
}//end of clustClass::iniProfile
*/

/*
void clustClass::addMotif(/*vector <motifClass> allMotifs,*//* int motifOrder)
{
	motifsID.push_back(motifOrder);
	updateProfile(1, motifOrder);
//218	cout <<"shift= "<<allMotifs[0].motifshift<<endl;
}

int clustClass::minusMotif(/*vector <motifClass> allMotifs,*//* int motifOrder)
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

void clustClass::updateProfile(int sign, /*vector <motifClass> allMotifs,*//* int motifOrder)
{
	int j,k;
	list <int>::iterator motifsIDIter;
	int motifWidth, newwidth;
	int motifshift,clustshift;
	bool begin = true;

	motifshift = allMotifs[motifOrder].motifshift;
	clustshift = allMotifs[motifOrder].clustshift;
	motifWidth = allMotifs[motifOrder].width;
//***** 02/18/06 ---
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
******** *//*
//***** 02/18/06 
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
//			profile[j].clear();
			for(k=0;k<newwidth;k++)
			{	
				allClustProfiles[j][start +k] = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
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
				int e= allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
				int f = allClustProfiles[j][start +k];
				allClustProfiles[j][start +k] = allClustProfiles[j][start +k] + 
					sign * allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
/* *****
				if(profile[j][k] <0)
				{
					for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
					{
						cout <<(*motifsIDIter)<<" ";	
					}
					cout <<endl;
					//int e = profile[j][k];
					//int f = allMotifs[motifOrder].profile[j][k+motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
 						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}
***** *//*
			}//end of k
		}//end of j
	}//end of if
	else if(width < newwidth)//need to add profile, increase only when there is just 1 motif in the cluster
	{
//-		for(j=0;j<4;j++)
//-		{
//-			for(k=width;k<newwidth;k++)
//-				profile[j].push_back(0);
//-		}//end of j
		for(j=0;j<4;j++)
		{
			for(k=0;k<width;k++)
			{
				//int g = profile[j][k];
				//int h = allMotifs[motifOrder].profile[j][k+motifshift];
				//-profile[j][k]=profile[j][k] - allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
				allClustProfiles[j][start +k] = allClustProfiles[j][start +k] - allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
/* *****
				if(profile[j][k] <0)
				{
					int e = profile[j][k];
					int f = allMotifProfiles[j][allMotifs[motifOrder].start + k + motifshift];
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
						"motifshift = "<<motifshift <<"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}//end of if
***** *//*
			}//end of k
			motifsIDIter = motifsID.begin();
			for(k=width;k<newwidth;k++)
			{
				//int g = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
				//profile[j][k] = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
				allClustProfiles[j][start +k] = allMotifProfiles[j][allMotifs[(*motifsIDIter)].start + k];
			}//end of k
		}//end of j
		cout <<"*****"<<endl;
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				cout << allClustProfiles[j][start +k] <<" ";
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
	else//need to subtract profile, since a inter motif was added to the cluster
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<newwidth;k++)
				allClustProfiles[j][start +k] = allClustProfiles[j][start +k + clustshift] 
					+ allMotifProfiles[j][allMotifs[motifOrder].start + k];
				//profile[j][k] = profile[j][k + clustshift] + allMotifProfiles[j][allMotifs[motifOrder].start + k];
			// need to update motif shift
		}//end of j
		for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
		{
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
			//int e = (*motifsIDIter);
			//int f = allMotifs[(*motifsIDIter)].motifshift;
			if((*motifsIDIter) != motifOrder)
				allMotifs[(*motifsIDIter)].motifshift = allMotifs[(*motifsIDIter)].motifshift + clustshift;
//				allMotifs[(*motifsIDIter).id].motifshift = allMotifs[(*motifsIDIter).id].motifshift + clustshift; 
//218			cout <<"shift= "<<allMotifs[0].motifshift<<endl;
		}//end of motifsIDIter
	}//end of else
	width = newwidth;
}//end of clustClass::updateProfile
*/

/*
void getProfile()
{
	int j,k;
	list <int>::iterator motifsIDIter;
	
	for(j=0;j<width;j++)
	{
		for(k=0;k<4;k++)
			profile[k].push_back(0);
	}
	for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
	{
	
	}
}
*/

double clustClass::loglikelihood()//const int stage)
{
	int j,k;
	int motifsSize;
	int la, lsum = 0;//, sumin,start;
	double logsum,tem;
	list <int>::iterator motifsIDIter;

	tem=gammaln(betasum)-sumbeta;
	motifsSize = motifsID.size();
	lsum=0;
	logsum =0;
 	for(j=0;j<4;j++)
		lsum = lsum + allClustProfiles[j][start];
	for(j=0;j<width;j++)
	{
//-		if((stage!=1)||(power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				//+int a = profile[k][j];
//				la = profile[k][j];
				la = allClustProfiles[k][start + j];
				logsum=logsum+gammaln( (double) la+beta[k]);
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
					logsum=logsum+gammaln( (double)la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) sumin + betasum);
			}//end of 
			logsum=logsum+(double) (motifsSize-1)*tem;
		}//end of else if
*/
	}//end of j
	logsum = logsum + (double) width *tem;
	return(logsum);
}//end of clustClass::loglikelihood

list <clustClass> allClusts;

class sortclustClass
{
private:
public:
	int order;
	list <int> motifsID;
	int size;
///	int left;
///	int right;
	int width;
//-	vector <bool> power;
///	double remain;
	double postprob;
	~sortclustClass()
	{
		motifsID.clear();
//-		power.clear();
	}
};

double clustRemoveMotif(const int motifOrder, const list <clustClass>::iterator clustpoint);

int readResult(const char resultFileName[], const char shiftFileName[],
			   const char widthFileName[], const double factor_0, const int MAX_SHIFT);

int readSeq(const char inputFileName[], int &widthmin);

void genMotifProfile(vector <string> seq, const int start, const int width, const int order);

void genClustProfile(const int start, const int width, list <int> motifsID);

void resultDisplay(const int motifCount, const double loglike,// const int bin, 
		   const double factor_0, list <clustClass> clusters, const char inputFileName[], 
		   const char outPutFileName[], const char cleanOutPutFileName[], const char cleanMemberFileName[],
		   const char memberFileName[], const char shiftFileName[],const char widthFileName[],
		   const double posprobacutoff, const double bayesratiocutoff, const double diftime);
			
double finalRatioSingle(list <clustClass>::iterator allClustsIter);

double bayesRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart);//, const int stage);

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart);//, const int stage);

void mqheap(HEAPSORT *p[], int total, int begin, int end);

static void sift(HEAPSORT *p[], int k, int j, int m);

int main(int argc, char **argv)
{
	int j;
	int widthmin,motifCount;
	int totalclust;
	list <clustClass>::iterator allClustsIter;
	double factor_0 =1.0;
	double loglike;
	char inputFileName[MAX_LENGTH]="12";//"cisred12";//"fasta.23";//wrong";//"test.txt";//"fasta.wrong";
    char memberFileName[MAX_LENGTH]= "member.12";//.all";//"mot.23";//wrong";//"test.mot";//.wrong";//"mot.wrong";
    char shiftFileName[MAX_LENGTH]="shift.12";//.all";//"shift.23";//wrong";//"test.shift";//"shift.wrong";
	char widthFileName[MAX_LENGTH]="width.12";
	char outputFileName[MAX_LENGTH]="out.txt";
    	char cleanoutFileName[MAX_LENGTH]="cleanout.txt";
        char cleanMemberFileName[MAX_LENGTH]="out.clean.txt";
	double posprobacutoff,bayesratiocutoff;
        time_t start,end;
	double dif;
	int MAX_SHIFT =2;

        if(argc<11)
        {
          printf("10 options need to be specified:\n\tinput sequence filename,\n\tinput member filename,\n\tinput shift filename,\n\tinput width filename,\n\toutput filename,\n\toutput cleanup filename,\n\tmember cleanup filename,\n\tposterior probability cutoff,\n\tBayes ratio cutoff,\nMaximum allowed shifts.\n");
            exit(0);
        }
        for(j=0;j<MAX_LENGTH;j++)
        {
            inputFileName[j]=argv[1][j];
            memberFileName[j]=argv[2][j];
            shiftFileName[j]=argv[3][j];
            widthFileName[j]=argv[4][j];
            outputFileName[j]=argv[5][j];
	    cleanoutFileName[j]=argv[6][j];
	    cleanMemberFileName[j]=argv[7][j];
        }
	posprobacutoff = atof(argv[8]);
	bayesratiocutoff = atof(argv[9]);
        MAX_SHIFT = atoi(argv[10]);
	time (&start);  
	motifCount=readSeq(inputFileName,widthmin);
	cout <<"motifCount ="<<motifCount<<endl; 
	
//	read in result.
	totalclust = readResult(memberFileName,shiftFileName,widthFileName,factor_0,MAX_SHIFT);//readResult(memberFileName,shiftFileName,factor_0,MAX_SHIFT);
// calculate loglikelihood
	loglike = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end();allClustsIter ++)
	{
		loglike = loglike + (*allClustsIter).loglikelihood();
	}
//+	cout <<loglike<<endl;
        time(&end); 
        dif = difftime (end,start) /60;
// display results
	resultDisplay(motifCount, loglike, factor_0, allClusts,inputFileName, outputFileName,cleanoutFileName,cleanMemberFileName,memberFileName, shiftFileName, widthFileName, posprobacutoff,bayesratiocutoff,dif);
	delete [] allMotifProfiles;
	return(0);
}//end of main

/*
double clustRemoveMotif(const int motifOrder, const list <clustClass>::iterator clustpoint)//, const int stage)
{
//	int oldwidth = 0;
//	int newwidth = 0;
	double likenew,likeold;

//	oldwidth = (*clustpoint).width;
	likeold = (*clustpoint).loglikelihood();//0);
	(*clustpoint).minusMotif(motifOrder);
	likenew = (*clustpoint).loglikelihood();//0);
	return(likenew-likeold);
	cout <<"remove: old= "<<likeold<<"new= "<<likenew<<endl;
}// end of clustRemoveMotif
*/

int readResult(const char resultFileName[], const char shiftFileName[],
			   const char widthFileName[], const double factor_0,const int MAX_SHIFT)
{
	int j,k;
	int count; 
	int motif,shift;
	istringstream iss, issshift,isswidth;
	ifstream inFile, inshiftFile,inwidthFile;
	string id,lineString;
	string shiftlineString;
	string widthlineString;
	clustClass temClust;
	int totalclust;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustplace;
	list <int>::iterator motifsIDIter;
	double ratio,numerator;
	double sum;//,likeminus;
	bool single;
	int motifOrder;
	int motifWidth, clustWidth;
//+	int minwidth;
	int clustwidth;
	int clustorder =0;
	int start,width,totalcol;

	cout << "start reading in result" << endl;
	inFile.open(resultFileName);//.c_str());
	if(!inFile)
	{
		cout << "ERROR: Unable to open file: " << resultFileName << endl;
	    exit(1);
	}//end of if
        inshiftFile.open(shiftFileName);//.c_str());
        if(!inshiftFile)
        {
                cout << "ERROR: Unable to open file: " << shiftFileName << endl;
            exit(2);
        }//end of if
		inwidthFile.open(widthFileName);//.c_str());
	    if(!inwidthFile)
        {
                cout << "ERROR: Unable to open file: " << widthFileName << endl;
            exit(3);
        }//end of if
	count =0;
//	for(j=0;j<allMotifs.size();j++)
//		cout <<j<<" "<<allMotifs[j].width<<endl;
//	exit(0);
	start =0;
	while (inFile)
	{
// get name
		getline(inFile,lineString);
		getline(inshiftFile,shiftlineString); 
		getline(inwidthFile,widthlineString); 
		iss.clear();
		issshift.clear();
		isswidth.clear();
		iss.str(lineString);// + " ");
		issshift.str(shiftlineString);// + " ");
		isswidth.str(widthlineString);
		isswidth >> clustwidth;
		temClust.width = clustwidth;
//+		cout << temClust.width<<endl;
		while (iss)
		{
			iss >> motif;
			if (iss)
			{
				allMotifs[motif].clust = count;
				issshift >> shift;
				allMotifs[motif].motifshift = shift;
				temClust.motifsID.push_back(motif);
//				cout <<id<<" "<<motif<<endl;
/*			
				if(temClust.motifsID.size() ==0)
				{
					temClust.motifsID.push_back(motif);
					minwidth = allMotifs[motif].width;
					//temClust.iniProfile();
				}
				else
				{
					if(allMotifs[motif].width < minwidth)
						minwidth = allMotifs[motif].width;  
					//temClust.addMotif(motif);
					temClust.motifsID.push_back(motif);
				}
*/
			}//end of if
		}//end of while
		if(temClust.motifsID.size()>0)
		{	
//+			temClust.width = minwidth;
			temClust.start = start;
			allClusts.push_back(temClust);
//			genClustProfile(start,minwidth,temClust.motifsID);
			start = start + temClust.width;//minwidth;
			temClust.motifsID.clear();
			count ++;
		}
	}//end of while
	inFile.close();
	totalclust = allClusts.size();
	totalcol = start;
//+	cout <<"total column ="<< totalcol<<endl;

// get cluster profiles.
	allClustProfiles = new int*[4];
	for(j=0;j<4;j++)
	{
		allClustProfiles[j] = new int[totalcol];
		for(k=0;k<totalcol;k++)
			allClustProfiles[j][k]=0;
	}
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end();allClustsIter ++)
	{
		start = (*allClustsIter).start;
		width = (*allClustsIter).width;
		genClustProfile(start,width,(*allClustsIter).motifsID);
	}
/*
// calculate bayes ratio:	
	int newcount =0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end();allClustsIter ++)
	{
		newcount ++;
		(*allClustsIter).printClust();
		int a = (*allClustsIter).width;
	}
	cout <<newcount <<endl;
*/
// calculate posterior probabilities
	for(j=0;j<allMotifs.size();j++)
	{
		single = false;
		sum = factor_0;
		motifOrder = j;
//moved outside of loop		
		motifWidth = allMotifs[motifOrder].width;
		clustorder = 0;
//-		cout <<"j= "<<j<<endl;
//		if(j==34)
//			cout <<"here.";
		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end();allClustsIter ++)
		{
//			(*allClustsIter).printClust();
/*			if((*allClustsIter).motifsID.size()>1)
			{
				(*allClustsIter).minusMotif(motifOrder);
				clustplace = allClustsIter;
			}
*/
//			if(clustorder == 57)
//				cout <<"here.\n";
// decide shift
			clustWidth = (*allClustsIter).width;
			if((clustWidth >motifWidth)&&(abs(clustWidth -motifWidth)<=MAX_SHIFT))//cluster is longer
			{
//				cout <<"weird.\n";
//				clustlong.push_back(true);
				shift = clustWidth - motifWidth +1;
//				max = 0;
				for(k=0;k<shift;k++)
				{
//					motifStart = 0;
//					clustStart = k;
					ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,0,k));//,stage));
					sum = sum + ratio;
					//if((j==0)||(ratio >max))//removed 09/24/05
					//	max = ratio;			//removed 09/24/05
					//ratioVec.push_back(max); //removed 09/24/05
//					ratioVec.push_back(ratio);
				}//end of k
			}//end of if
			else if ((clustorder == allMotifs[motifOrder].clust)||((clustWidth <= motifWidth)&&(abs(clustWidth -motifWidth)<=MAX_SHIFT)))///motif is longer
			{
//				clustlong.push_back(false);
				if((clustorder == allMotifs[motifOrder].clust)&&((*allClustsIter).motifsID.size()==1))
				{
					single = true;
					ratio=0;
				}
				else 
				{
					int a = allMotifs[motifOrder].clust;
					if(clustorder == allMotifs[motifOrder].clust)
					{
					//	(*allClustsIter).minusMotif(motifOrder);
					//	clustplace = allClustsIter;
					}
					shift = motifWidth - clustWidth +1;
					for(k=0;k<shift;k++)
					{	
						if(clustorder == allMotifs[motifOrder].clust)
							ratio=exp(bayesMinusRatioSingle(motifOrder, allClustsIter,k,0));//,stage));
						else
							ratio=exp(bayesRatioSingle(motifOrder, allClustsIter,k,0));//,stage));
						sum = sum + ratio;
						int a = allMotifs[motifOrder].clust;
						int b = allMotifs[motifOrder].motifshift;
						if((clustorder == allMotifs[motifOrder].clust)&&(k == allMotifs[motifOrder].motifshift))
							numerator = ratio;
					}//end of k

				}//end of else
			}//end of else if
			clustorder ++;
		}//end of allClustsIter   
		if(single)
			allMotifs[motifOrder].postprob = factor_0/sum;
		else
		{
			allMotifs[motifOrder].postprob = numerator/sum;
//			(*clustplace).addMotif(motifOrder);
		}//end of else
	}//end of j
	return(totalclust);	
}//end of readResult

int readSeq(const char inputFileName[], int &widthmin)
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
			temMotif.clustshift = 0;
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
			genMotifProfile(seq,start,width,allMotifs[count].depth);
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

void genMotifProfile(vector <string> seq, const int start, const int width, const int depth)
{
	int j,k;
	
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
}//end of genMotifProfile

void genClustProfile(const int start, const int width, list <int> motifsID)
{
	int j,k;
	list <int>::iterator motifsIDIter;
	
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
		{
			for(motifsIDIter=motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter++)
			{
				int a = allClustProfiles[j][k + start];
				int b = allMotifs[(*motifsIDIter)].start;
				int c = allMotifs[(*motifsIDIter)].motifshift;
				int d = allMotifProfiles[j][k + allMotifs[(*motifsIDIter)].start + allMotifs[(*motifsIDIter)].motifshift];
				allClustProfiles[j][k + start] = allClustProfiles[j][k + start] 
					+ allMotifProfiles[j][k + allMotifs[(*motifsIDIter)].start + allMotifs[(*motifsIDIter)].motifshift];
			}//end of m
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
}//end of genMotifProfile

double finalRatioSingle(list <clustClass>::iterator allClustsIter)
{
    int j,k, begin;
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
		suminclust=suminclust + allClustProfiles[j][(*allClustsIter).start];
	sumex=0;
	difsum=0;
	motifsSize = (*allClustsIter).motifsID.size();
	if(motifsSize==1)
	  return 0;
	else
	{
//		for(k=0;k<4;k++)
//			cout <<(*allClustsIter).profile[k][0]<<endl;
		for(j=0;j<width;j++)
		{	
//-			if((*allClustsIter).power[j])//==1
//-			{
				sumex=0;
				for(k=0;k<4;k++)
				{	
					sumex=sumex+gammaln(allClustProfiles[k][(*allClustsIter).start + j]+beta[k]);
				}
				oldProfile=sumex-gammaln((double) suminclust+betasum);
				newProfile=0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
					begin = allMotifs[(*motifsIDIter)].motifshift;
					sumex=0;
					for(k=0;k<4;k++)
					{
						la = allMotifProfiles[k][allMotifs[(*motifsIDIter)].start + j +begin];
						sumex=sumex+gammaln((double) la+beta[k]);
					}//end of m
					newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
			    }//end of k
				difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
//-			}//end of if
		}//end of j
		difsum=difsum-(double) width/*remain*/*(double)(motifsSize-1)*tem;
		return(difsum);
	}//end of else
}//end of finalRatioSingle

//+void resultDisplay(const int motifCount, const double loglike,// const int bin, 
//+		   const double factor_0, list <clustClass> clusters, const char outPutFileName[],
//+		   const char cleanOutPutFileName[],const double cutoff, const double diftime)//-01/18/07, const char motifFileName[])
void resultDisplay(const int motifCount, const double loglike,// const int bin,
                   const double factor_0, list <clustClass> clusters, const char inputFileName[], 
                   const char outPutFileName[], const char cleanOutPutFileName[], const char cleanMemberFileName[],
                   const char memberFileName[], const char shiftFileName[],const char widthFileName[],
                   const double posprobacutoff, const double bayesratiocutoff, const double diftime)
{
	int j;//,k,m;
	list <clustClass>::iterator clustIter;
	list <int>::iterator motifsIDIter;
	ofstream outPutFile;
	ofstream cleanOutPutFile;
        ofstream cleanMemberFile;
//+	ofstream motifblockFile;// FASTA format motif sequence alignment file
	int finalClustNumber;
	int tot,count;
	double value =0;
	int motifsize;
	sortclustClass *clustarray, **sorted;
//	int start;
	int minWidth, maxWidth;
	bool begin;
	int *keeper;
	vector <int> leftover;
	int leftovercount;
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
	cleanOutPutFile.open(cleanOutPutFileName);
	if(!cleanOutPutFile)
	{
		cout << "ERROR: Unable to open file: " << cleanOutPutFileName<<endl;
	    exit(5);
	}//end of if
	cleanMemberFile.open(cleanMemberFileName);
        if(!cleanMemberFile)
        {
                cout << "ERROR: Unable to open file: " << cleanMemberFileName<<endl;
            exit(8);
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
	outPutFile << "*       BMC - VAR Result Display       * \n";
	outPutFile << "*                                      * \n";
	outPutFile << "**************************************** \n\n";
	
	cleanOutPutFile << "**************************************** \n";
        cleanOutPutFile << "*                                      * \n";
        cleanOutPutFile << "*     BMC -VAR Clean Result Display    * \n";
        cleanOutPutFile << "*                                      * \n";
        cleanOutPutFile << "**************************************** \n\n";

        outPutFile <<"***************************************************"<<endl;
        outPutFile <<"The local time now is "<< asctime(timeinfo) <<endl;
        outPutFile <<"The following parameters are in effect:"<<endl;
        outPutFile <<"Input file name:           " <<inputFileName<<endl;
        outPutFile <<"Output file name:          " <<outPutFileName<<endl;
        outPutFile <<"Cleaned Output file name:  " <<cleanOutPutFileName<<endl;
        outPutFile <<"Member file name:          " <<memberFileName<<endl;
        outPutFile <<"Shift position file name:  " <<shiftFileName<<endl;
        outPutFile <<"Motif width file name:     " <<widthFileName<<endl;
	outPutFile <<"Posteriro proba cutoff:    " <<posprobacutoff<<endl;
        outPutFile <<"Bayes ratio cutoff:	 " <<bayesratiocutoff<<endl;
        outPutFile <<"Values of alpha:           " <<factor_0<<endl;
        outPutFile <<"(Tunning parameter in DP model)."<<endl<<endl;
        outPutFile <<"It took BMC-VAR Result Display "<<diftime<<" minutes to finish this job."<<endl;
        outPutFile <<"***************************************************"<<endl<<endl;

        cleanOutPutFile <<"***************************************************"<<endl;
        cleanOutPutFile <<"The local time now is "<< asctime(timeinfo) <<endl;
        cleanOutPutFile <<"The following parameters are in effect:"<<endl;
        cleanOutPutFile <<"Input file name:           " <<inputFileName<<endl;
        cleanOutPutFile <<"Output file name:          " <<outPutFileName<<endl;
        cleanOutPutFile <<"Cleaned Output file name:  " <<cleanOutPutFileName<<endl;
        cleanOutPutFile <<"Member file name:          " <<memberFileName<<endl;
        cleanOutPutFile <<"Shift position file name:  " <<shiftFileName<<endl;
        cleanOutPutFile <<"Motif width file name:     " <<widthFileName<<endl;
        cleanOutPutFile <<"Posteriro proba cutoff:    " <<posprobacutoff<<endl;
        cleanOutPutFile <<"Bayes ratio cutoff:        " <<bayesratiocutoff<<endl;
        cleanOutPutFile <<"Values of alpha:           " <<factor_0<<endl;
        cleanOutPutFile <<"(Tunning parameter in DP model)."<<endl<<endl;
        cleanOutPutFile <<"It took BMC-VAR Result Display "<<diftime<<" minutes to finish this job."<<endl;
        cleanOutPutFile <<"***************************************************"<<endl<<endl;
	
	outPutFile << " Total number of all motifs: " << motifCount << "\n";
	cleanOutPutFile << " Total number of all motifs: " << motifCount << "\n";
    cout << "\n Total number of all motifs: " << motifCount << "\n";
	outPutFile << "\nThere are total of " << finalClustNumber << " clusters\n";
	cleanOutPutFile << "\nThere are total of " << finalClustNumber << " clusters\n";
    cout <<"\n There are total of of " << finalClustNumber << " clusters, see output file " << outPutFileName << " for details\n\n";
	outPutFile << "q=" << factor_0 << ", motif width range = " <<minWidth <<"--"<<maxWidth<<endl;
	cleanOutPutFile << "q=" << factor_0 << ", motif width range = " <<minWidth <<"--"<<maxWidth<<endl;
	outPutFile << "log likelihood = " << loglike <<endl<<endl;
    cleanOutPutFile << "log likelihood = " << loglike <<endl<<endl;
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
	keeper = new int[finalClustNumber];
    for (j=finalClustNumber;j>0;j--)
    {
//		outPutFile << "cluster " << finalClustNumber-j+1 << ", size= "<< (*sorted[j-1]).size <<"  ";// <<"Bayes Ratio=";
		outPutFile << "cluster " << finalClustNumber-j+1 << "["<<(*sorted[j-1]).order+1<<"], size= "<< (*sorted[j-1]).size <<"  ";
		outPutFile << "width = "<<(*sorted[j-1]).width;//02/18/06<<", number of On columns = " << (*sorted[j-1]).remain << "   ";
		outPutFile << " Average Bayes ratio = " << (*sorted[j-1]).postprob <<endl;
		tot=0;
		keeper[j] =0;
		for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
		{
			tot++;
			//if((tot%2==1)&&(tot!=1)) 
			//	outPutFile << endl;                      
			 outPutFile << tot << " <" <<(*motifsIDIter)<<"> "<< allMotifs[(*motifsIDIter)].name << " [" << allMotifs[(*motifsIDIter)].motifshift 
						<< "]  ("<< allMotifs[(*motifsIDIter)].postprob <<")"<<endl;
			if(allMotifs[(*motifsIDIter)].postprob > posprobacutoff)
				keeper[j] ++;
		}//end of motifsIDIter
		outPutFile <<endl<<endl;
	}//end of j
	// to find significant clusters.
    count =0;
	leftovercount =0;
    for (j=finalClustNumber;j>0;j--)
    {
		if((keeper[j]>1)&&((*sorted[j-1]).postprob > bayesratiocutoff))//>0
		{
                cleanOutPutFile << "cluster " << count+1 <<"("<<finalClustNumber-j+1 <<")"<< "["<<(*sorted[j-1]).order+1<<"], size= "<<keeper[j]<<"("<< (*sorted[j-1]).size<<")  ";
                cleanOutPutFile << "width = "<<(*sorted[j-1]).width;//02/18/06<<", number of On columns = " << (*sorted[j-1]).remain<<"  ";
                cleanOutPutFile << " Average Bayes ratio = " << (*sorted[j-1]).postprob <<endl;
                tot=0;
                for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
                {
					if(allMotifs[(*motifsIDIter)].postprob > posprobacutoff)
					{
	                        		tot++;
				            cleanOutPutFile << tot << " " << allMotifs[(*motifsIDIter)].name << " [" << allMotifs[(*motifsIDIter)].motifshift
                                                << "]  ("<< allMotifs[(*motifsIDIter)].postprob <<")"<<endl;
						cleanMemberFile << *motifsIDIter <<"  ";
					}//end of if
					else
					{
						leftover.push_back((*motifsIDIter));
						leftovercount ++;
					}//end of else
		         }//end of motifsIDIter
				 cleanOutPutFile <<endl<<endl;
				count ++;
			cleanMemberFile <<endl;
        }//end of if
		else
		{
			for(motifsIDIter=(*sorted[j-1]).motifsID.begin();motifsIDIter !=(*sorted[j-1]).motifsID.end();motifsIDIter++)
            {
				leftover.push_back((*motifsIDIter));
				leftovercount ++;
			}
		}//end of else
    }//end of j
	for(j=0;j<leftovercount;j++)
		cleanOutPutFile << j+1 <<" "<<leftover[j] <<" "<<allMotifs[leftover[j]].name <<endl;
	cleanOutPutFile << endl<<"There are " <<count<<" clusters left after clean up."<<endl;
	cleanOutPutFile << endl<<"There are " <<leftovercount<<" motifs left out after clean up."<<endl;
	cleanOutPutFile << endl<<"There are " <<motifCount - leftovercount<<" motifs passed clean up."<<endl;
	delete [] clustarray;
	delete [] sorted;
//	delete [] keeper;
	leftover.clear();
	outPutFile.close();
	cleanOutPutFile.close();
	cleanMemberFile.close();
}//end of resultDisplay

double bayesRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart)//, const int stage)
{
	int j,k;
	int Lold,Lz,Lnew,sum;
	double sumbayes,tem;
	double newProfile;
	double bayesRatio;
	int clustProfile,motifProfile;
	int width;
	
	int a = (*allClustsIter).width;
	int b = allMotifs[motifOrder].width;
//***** changed 03/01/07.
	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);//changed 9/13/05
//+	width = (*allClustsIter).width;	
// again changed on 05/08/07, since when motifwidth < cluster width, error will occur.
//***** changed 03/01/07.
	sum=0;

//	cout << "*****";
//	(*allClustsIter).printClust();
	for(j=0;j<4;j++)
	{
		sum=sum+allClustProfiles[j][(*allClustsIter).start];
	}
	Lold=sum;
	Lz=allMotifs[motifOrder].depth;
	Lnew=Lold+Lz;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
//-		if((stage !=1)||((*allClustsIter).power[j]))//removed ==1 08/10/05
//-		{
			for(k=0;k<4;k++)
			{
				clustProfile = allClustProfiles[k][(*allClustsIter).start + j + clustStart];
				motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
				newProfile = clustProfile + motifProfile;
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) clustProfile + beta[k])
					-gammaln((double) motifProfile + beta[k]);
			}//end of k
//-		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum)+gammaln((double) Lz+betasum)
		-gammaln(betasum)+sumbeta;
///	if(stage == 1)
///		bayesRatio = sumbayes+tem*(double) (*allClustsIter).remain;
///	else 
		bayesRatio = sumbayes+tem*(double) width;
	return(bayesRatio);
}//end of bayesRatioSingle

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter,int motifStart, int clustStart)//, const int stage)
{
	int j,k;
	int Lold,Lz,Lnew,sum;
	double sumbayes,tem;
	int newProfile,difProfile;
	double bayesRatio;
	int clustProfile,motifProfile;
	int width,shift;
	
	int a = (*allClustsIter).width;
	int b = allMotifs[motifOrder].width;
//***** changed 03/01/07.
//	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);//changed 9/13/05
	width = (*allClustsIter).width;
//***** changed 03/01/07.
/*
	list <int>::iterator motifsIDIter;
	for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
		cout <<(*motifsIDIter)<<endl;
	if(motifOrder ==0)
	{
		for(j=0;j<width;j++)
		{
			for(k=0;k<4;k++)
				cout <<allClustProfiles[k][(*allClustsIter).start + j+clustStart]<<" ";
			cout <<endl;
		}
	}
*/
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum + allClustProfiles[j][(*allClustsIter).start];
	}
	Lnew=sum;
	Lz=allMotifs[motifOrder].depth;
	Lold=Lnew - Lz;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<width;j++)
	{
//-		if((stage !=1)||((*allClustsIter).power[j]))//removed ==1 08/10/05
//-		{
			for(k=0;k<4;k++)
			{
//				clustProfile = (*allClustsIter).profile[k][j+clustStart];
//				motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
//				newProfile = clustProfile + motifProfile;
//				newProfile = (*allClustsIter).profile[k][j+clustStart];

				motifProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + motifStart];
				shift = allMotifs[motifOrder].motifshift;
				difProfile = allMotifProfiles[k][allMotifs[motifOrder].start + j + shift];
				clustProfile = allClustProfiles[k][(*allClustsIter).start + j+clustStart] - difProfile;
				newProfile = clustProfile + motifProfile;
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) clustProfile + beta[k])
					-gammaln((double) motifProfile + beta[k]);
			}//end of k
//-		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum)+gammaln((double) Lz+betasum)
		-gammaln(betasum)+sumbeta;
///	if(stage == 1)
///		bayesRatio = sumbayes+tem*(double) (*allClustsIter).remain;
///	else 
		bayesRatio = sumbayes+tem*(double) width;
	return(bayesRatio);
}//end of bayesMinusRatioSingle

void mqheap(HEAPSORT *p[], int total, int begin, int end)
{
	int j,k;
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

static void sift(HEAPSORT *p[], int k, int j, int m)
{
	int a,b;
	double t; //int
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


