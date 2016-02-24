// clustering program BMC 2.0
// origin from bmc by Steve Qin, 01/15/02
// modified on 03/16/04 using STL list and map
// modified 02/18/06 to work for long list
// for fixed width motifs only. 
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
//	short id;
	string name;
//	bool isEven;
	int start;
	short width;
	short depth;
//	short motifshift;
//	short clustshift;
	double postprob;
	double ratio;
//	vector <string> basepair;
//	vector <short> profile[4];
//	void printSeq();
//	void genProfile();
//	friend void readSeq(list <clustClass > &clustersList);
};

vector <motifClass> allMotifs;
int **allProfile;
short fixedwidth;

class clustClass
{
private:
public:
	list <short> motifsID;
//	vector <bool> power;
//-	short width;
//	short left;
//	short right;
//	double remain;
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
//	short left;
//	short right;
	short width;
//	vector <bool> power;
//	double remain;
	double postprob;
	~sortclustClass()
	{
		motifsID.clear();
	}
};

/*
void motifClass::printSeq()
{
	short j;

	for(j=0;j<depth;j++)
		cout << basepair[j] <<endl;
}
*/

void clustClass::printClust()
{
//	short j,k;
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
		cout << "Name " << allMotifs[(*motifsIDIter)].name; 
		//<< " " <<allMotifs[(*motifsIDIter)].width << 
		//	" " << allMotifs[(*motifsIDIter)].depth << " " << allMotifs[(*motifsIDIter)].isEven << endl;
		count++;
	}//end of motifsIDIter
	cout <<endl;
/*
	for(j=0;j<4;j++)
	{
		for(k=0;k<fixedwidth;k++)
		{
//			outputFile << profile[j][k] <<" ";
			cout << profile[j][k] << " ";
		}
		cout << endl;
//		outputFile << endl;
	}//end of j
*/
//	outputFile << endl;
//	outputFile.close();
}//end of clustClass::prshortClust()

//----- initialize the profile matrix when the first motif join this cluster

void clustClass::iniProfile()//short bin)
{
	short j,k;
	list <short>::iterator motifsIDIter;

	motifsIDIter=motifsID.begin();
//-	width=allMotifs[(*motifsIDIter)].width;

// first set profile matrix to 0

	for(j=0;j<4;j++)
	{
		for(k=0;k<fixedwidth;k++)
			profile[j].push_back(0);
	}
	for(j=0;j<4;j++)
	{
		for(k=0;k<fixedwidth;k++)
		{
			//short a = allMotifs[(*motifsIDIter)].profile[j][k];
			profile[j][k]=allProfile[j][(*motifsIDIter)*fixedwidth + k];//+(*motifsIDIter).profile[kstar[k]][motifWidth-1-j-start];
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
	bool begin = true;
	int startpos;

	startpos = motifOrder*fixedwidth;
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

	if((sign == -1)&&(motifsID.size()==1))
	{
		motifsIDIter = motifsID.begin();
		for(j=0;j<4;j++)
		{
			profile[j].clear();
			for(k=0;k<fixedwidth;k++)
			{	
				profile[j].push_back(allProfile[j][(*motifsIDIter)*fixedwidth + k]);
//				profile[j].push_back(allMotifs[(*motifsIDIter)].profile[j][k]);
//218				cout << profile[j][k] << " ";	
			}//end of k
//218			cout <<endl;
		}//end of j

	}//end of if
	else 
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<fixedwidth;k++)
			{
//				short e= allMotifs[motifOrder].profile[j][k];
				short e= allProfile[j][startpos + k];//03/03/07 "startpos" replace "motifOrder*fixedwidth".
				short f = profile[j][k];
				profile[j][k]=profile[j][k]+sign * allProfile[j][startpos + k];//03/03/07 "startpos" replace "motifOrder*fixedwidth".//allMotifs[motifOrder].profile[j][k+motifshift];
				if(profile[j][k] <0)
				{
					short e = profile[j][k];
					short f = allProfile[j][startpos + k];//03/03/07 "startpos" replace "motifOrder*fixedwidth".
					cout <<"error in cluster profile, motifOrder = " <<motifOrder <<"j = "<<j<<"k = " <<k<<
 						"profile+ " <<e<<"motif profile= "<<f<<endl;
					exit(15);
				}
			}//end of k
		}//end of j
	}//end of if
}//end of clustClass::updateProfile

double clustClass::loglikelihood()//const short stage)
{
	short j,k;
	short motifsSize;
	short la, lsum = 0;//, sumin,start;
	double logsum;//,tem;
	list <short>::iterator motifsIDIter;

//	tem=gammaln(betasum)-sumbeta;
	motifsSize = motifsID.size();
	lsum=0;
	logsum =0;
 	for(j=0;j<4;j++)
		lsum = lsum + profile[j][0];
	for(j=0;j<fixedwidth;j++)
	{
//-		if((stage!=1)||(power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				la = profile[k][j];
				logsum=logsum+gammaln( (double)la+beta[k]);
			}//end of k
			logsum=logsum-gammaln((double) lsum + betasum);
//-		}//end of if
/*	
		else if ((stage==1)&&(!power[j]))
		{
			for(motifsIDIter = motifsID.begin();motifsIDIter != motifsID.end();motifsIDIter ++)
			{
				sumin=allMotifs[(*motifsIDIter)].depth;
//-				start=allMotifs[(*motifsIDIter)].motifshift;// add (*motifsIDIter).shift 
				for(k=0;k<4;k++)
				{
//-					la = allMotifs[(*motifsIDIter)].profile[k][j + start];//+beta[k];
					la = allProfile[k][(*motifsIDIter)*fixedwidth + j];
					logsum=logsum+gammaln( (double) la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) sumin + betasum);
			}//end of 
			logsum=logsum+(double) (motifsSize-1)*tem;
		}//end of else if
*/
	}//end of j
	logsum = logsum + (double) fixedwidth *betasumdif;//tem;
	return(logsum);
}//end of clustClass::loglikelihood

list <clustClass> allClusts;
vector <list <clustClass>::iterator> motifHome;

short readSeq(const char inputFileName[]);//, short &widthmin)

void genProfile(vector <string> seq, const int start, const int motifOrder);

//void genProfile(vector <string> seq, const int start, const int width, const int depth);

void initialCluster(/*const short motifCount, */const short INITIAL_CLUSTER);

double bmc(const int initialclusternum, const short numChains, const short numItera, const double factor_0, list <clustClass> &finalClusters);

double loglikelihoodSingle();//const short stage);

double chain(const short NUM_CYCLE, const double FACTOR_0, list <clustClass> &bestClusters);

double operation(const short motifOrder, const double FACTOR_0);//, short &del);
				 //, const short stage);

//+ double clustInsertMotif(const short motifOrder, const short decision);//, const short stage);

double clustInsertMotif(const short motifOrder, const list <clustClass>::iterator clustIter);

double clustRemoveMotif(const short motifOrder, const list <clustClass>::iterator clustPoint);//, const short stage);

double loglikelihoodUpdateSingle(const short motifOrder, const list <clustClass>::iterator allClustsIter, const short stage);

//double diffWidthLoglikelihoodUpdateSingle (const short oldwidth, const short newwidth, const short motifOrder, 
//										   const list <clustClass>::iterator allClustsIter,const short stage);

//double onoffMetro(const list <clustClass>::iterator allClustsIter);

//double metSingle(const list <clustClass>::iterator allClustsIter, short indel, short position, short W);

//double fragment();

double bayesRatioSingle(short motifOrder, list <clustClass>::iterator allClustsIter);
						//,short motifStart, short clustStart, const short stage);

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter);

short decide(vector <double> &ratioVec, const double FACTOR_0, double &postprob);

short shiftdecide(vector <vector <double> > &ratioAll, const double FACTOR_0,short &shift, double &postprob);

double unran(int *na, int *nb, int *nc);

//double emptyLoglikelihoodupdateSingle(const short motifOrder);//, const list <clustClass>::iterator allClustsIter);
//										//, const short stage);

//void initWidth(const short bin, const short W);

//void resultDisplay(const short motifCount, const double loglike,
//		   const double factor_0, list <clustClass> clusters, const char outPutFileName[],
//                   const char memberFileName[]);

void resultDisplay(const short motifCount, const double loglike, const double factor_0, list <clustClass> clusters,
                   const char inputFileName[], const char outPutFileName[], const char memberFileName[],
                   const short numChains, const short numIter, const short initialclusternum, const double diftime);

void mqheap(HEAPSORT *p[], short total, short begin, short end);

static void sift(HEAPSORT *p[], short k, short j, short m);

double finalRatioEven(list <clustClass>::iterator allClustsIter);

double finalRatioOdd(list <clustClass>::iterator allClustsIter);

double finalRatioSingle(list <clustClass>::iterator allClustsIter);

//short checkWidth(	list <clustClass>::iterator clustPoint);

int main(int argc, char **argv)
{
	short j;
	list <clustClass> finalClusters;
//	vector <motifClass> fullMotifVector;
//-	short W,WMIN,WMAX;
//-	short  MAX_SHIFT;
	short numChains,numIter;
	double factor_0;
	short motifCount;
//	short bin=0;
	double loglike;
//	short IN_WIDTH;
	char inputFileName[MAX_LENGTH]="cisred7";//"60.txt";//"probfile2.txt";
        char outputFileName[MAX_LENGTH]="out.txt";
//-        char motifFileName[MAX_LENGTH]="motifblock.txt";
        char memberFileName[MAX_LENGTH]="member.txt";
	short initialclusternum =100;
//-	short widthmin;
	time_t start,end;
	double dif;

//	IN_WIDTH = 6;//5;//20;
//	W=IN_WIDTH;//6;
//	WMIN=6;//5;//16;//16;//4;//less than WMAX if consider fragmentation
//	WMAX=W;//16;//6;
//-	MAX_SHIFT=1;//5;
//	numChains=1;
//	numIter=5;
	factor_0=1.0;
//-	bin=2;//0;

//	srand((unsigned)time(NULL));

        if(argc<7)
        {
          printf("6 options need to be specified:\n\tinput filename,\n\toutput filename,\n\tcluster member filename,\n\tnumber of chains,\n\tnumber of iterations,\n\tinitial number of clusters\n");
            exit(0);
        }
        for(j=0;j<MAX_LENGTH;j++)
        {
            inputFileName[j]=argv[1][j];
            outputFileName[j]=argv[2][j];
//-            motifFileName[j]=argv[3][j];
            memberFileName[j]=argv[3][j];
        }
		
	numChains = atoi(argv[4]);
	numIter = atoi(argv[5]);
	initialclusternum = atoi(argv[6]);
// read in data, initial assign clusters and calculate cluster profiles
//	motifCount=readSeq(fullMotifVector,inputFileName,initialclusternum,widthmin);
	motifCount=readSeq(inputFileName);
//-	W = widthmin;
//-	WMIN = W;
//-	WMAX = W;
//+++	initialCluster(initialclusternum);//motifCount
	//exit(0);
// clustering
//-	cout << "cluster size= " << allClusts.size();
	time (&start);
	loglike = bmc(initialclusternum,numChains,numIter,factor_0,finalClusters);
	time(&end);
	dif = difftime (end,start) /60;
	allClusts.clear();
// output result
	resultDisplay(motifCount, loglike, factor_0, finalClusters, inputFileName, outputFileName, memberFileName, 
			numChains, numIter, initialclusternum,dif);
	finalClusters.clear();
	return(0);
}//end of main

void resultDisplay(const short motifCount, const double loglike, const double factor_0, list <clustClass> clusters, 
		   const char inputFileName[], const char outPutFileName[], const char memberFileName[], 
		   const short numChains, const short numIter, const short initialclusternum, const double diftime)
{
	short j;//,k,m;
	list <clustClass>::iterator clustIter;
	list <short>::iterator motifsIDIter;
	ofstream outPutFile;
//-	ofstream motifblockFile;// FASTA format motif sequence alignment file.
        ofstream memberFile;// cluster member motif information.
	short finalClustNumber;
	short tot,count;
	double value =0;
	short motifsize;
	sortclustClass *clustarray, **sorted;
//-	short start;
//-	short minWidth, maxWidth;
//-	bool begin;
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

	finalClustNumber = clusters.size();
// find out minimum and maximum cluster width
/*
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
*/
	outPutFile << "**************************************** \n";
	outPutFile << "*                                      * \n";
	outPutFile << "*       BMC 2.1 clustering Result      * \n";
	outPutFile << "*                                      * \n";
	outPutFile << "**************************************** \n\n";

        outPutFile <<"***************************************************"<<endl;
	outPutFile <<"The local time now is "<< asctime(timeinfo) <<endl;
	outPutFile <<"The following parameters are in effect:"<<endl;
	outPutFile <<"Input file name:           " <<inputFileName<<endl;
        outPutFile <<"Output file name:          " <<outPutFileName<<endl;
        outPutFile <<"Member file name:          " <<memberFileName<<endl;
        outPutFile <<"Number of chains:          " <<numChains<<endl;
        outPutFile <<"Number of cycles:          " <<numIter<<endl; 
        outPutFile <<"Number of initial clusters:" <<initialclusternum<<endl; 
        outPutFile <<"Values of alpha:           " <<factor_0<<endl;     
        outPutFile <<"(Tunning parameter in DP model)."<<endl<<endl;
	outPutFile <<"It took BMC-FIX "<<diftime<<" minutes to finish this job."<<endl;
	outPutFile <<"***************************************************"<<endl<<endl;

	outPutFile << " Total number of all motifs: " << motifCount << "\n";
	    cout << "\n Total number of all motifs: " << motifCount << "\n";
	outPutFile << "\nThere are total of " << finalClustNumber << " clusters\n";
    cout <<"\n There are total of of " << finalClustNumber << " clusters, see output file " << outPutFileName << " for details\n\n";
	outPutFile << "q=" << factor_0 <<", motif width range = " <<fixedwidth <<"--"<<fixedwidth<<endl;
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
//-		clustarray[count].left = (*clustIter).left;
//-		clustarray[count].right = (*clustIter).right;
//-		clustarray[count].remain = (*clustIter).remain;
		clustarray[count].motifsID = (*clustIter).motifsID;
//-		clustarray[count].power = (*clustIter).power;
//-		clustarray[count].width = (*clustIter).width;
		value=finalRatioSingle(clustIter);
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
		outPutFile << "width = "<<fixedwidth;//02/18/06<<", number of On columns = " << (*sorted[j-1]).remain << "   ";
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
		        outPutFile << tot << " " << allMotifs[(*motifsIDIter)].name << " ("<< 
						allMotifs[(*motifsIDIter)].postprob <<")"<<endl;
			memberFile << (*motifsIDIter) <<"  ";
		}//end of motifsIDIter
                outPutFile <<endl<<endl;
                memberFile <<endl;
	}//end of j
	delete [] clustarray;
	delete [] sorted;
	outPutFile.close();
//-	motifblockFile.close();
        memberFile.close();
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

double finalRatioSingle(list <clustClass>::iterator allClustsIter)
{
    short j,k;//, start;
	list <short>::iterator motifsIDIter;
//  	double oldProfile, newProfile;
	double sumin, oldsumex,newsumex;
	short suminclust,la;
	double difsum;
//	double tem;
	short motifsSize;
//-	short width;

//-	width = (*allClustsIter).width;
//	tem=gammaln(betasum)-sumbeta;
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
		for(j=0;j<fixedwidth;j++)
		{
//-			if((*allClustsIter).power[j])//==1
//-			{
				oldsumex=0;
				for(k=0;k<4;k++)
				{	
					oldsumex = oldsumex+gammaln((*allClustsIter).profile[k][j]+beta[k]);
				}
				//oldProfile=sumex-gammaln((double) suminclust+betasum);
				//newProfile=0;
				newsumex=0;
				for(motifsIDIter=(*allClustsIter).motifsID.begin();motifsIDIter !=(*allClustsIter).motifsID.end();motifsIDIter++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
//-					start = allMotifs[(*motifsIDIter)].motifshift;
					//sumex=0;
					for(k=0;k<4;k++)
					{
//						la = allMotifs[(*motifsIDIter)].profile[k][j+start];
						la = allProfile[k][(*motifsIDIter)*fixedwidth + j];
						newsumex = newsumex+gammaln((double) la+beta[k]);
					}//end of m
					newsumex = newsumex-gammaln((double) sumin+betasum);
			    }//end of k
//				difsum=difsum+oldProfile-newProfile;
				difsum=difsum+oldsumex-newsumex;
			//printf("%10.5f,%10.5f\n",oldProfile,newProfile);
//-			}//end of if
		}//end of j
//-		difsum=difsum-(double)((*allClustsIter).remain)*(double)(motifsSize-1)*tem;
		difsum=difsum-(double) fixedwidth*(gammaln((double) suminclust+betasum)+ (double)(motifsSize-1)*betasumdif);//tem;
		return(difsum);
	}//end of else
}//end of finalRatioSingle

/*
double finalRatioSingle(list <clustClass>::iterator allClustsIter)
{
    short j,k;//, start;
	list <short>::iterator motifsIDIter;
  	double oldProfile, newProfile;
	double sumin, sumex;
	short suminclust,la;
	double difsum;
//	double tem;
	short motifsSize;
//-	short width;

//-	width = (*allClustsIter).width;
//	tem=gammaln(betasum)-sumbeta;
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
		for(j=0;j<fixedwidth;j++)
		{
//-			if((*allClustsIter).power[j])//==1
//-			{
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
//-					start = allMotifs[(*motifsIDIter)].motifshift;
					sumex=0;
					for(k=0;k<4;k++)
					{
//						la = allMotifs[(*motifsIDIter)].profile[k][j+start];
						la = allProfile[k][(*motifsIDIter)*fixedwidth + j];
						sumex=sumex+gammaln((double) la+beta[k]);
					}//end of m
					newProfile=newProfile+sumex-gammaln((double) sumin+betasum);
			    }//end of k
				difsum=difsum+oldProfile-newProfile;
	    //printf("%10.5f,%10.5f\n",oldProfile,newProfile);
//-			}//end of if
		}//end of j
//-		difsum=difsum-(double)((*allClustsIter).remain)*(double)(motifsSize-1)*tem;
		difsum=difsum-(double) fixedwidth*(double)(motifsSize-1)*betasumdif;//tem;
		return(difsum);
	}//end of else
}//end of finalRatioSingle
*/

double bmc(const int initialclusternum, const short numChains, const short numItera, const double factor_0, list <clustClass> &finalClusters)
{
	short j;
	double like,large;
	list <clustClass> bestClusters;

	cout << "start iteration"<<endl;
	for(j=0;j<numChains;j++)
	{
		allClusts.clear();
		motifHome.clear();
	        initialCluster(initialclusternum);
		cout << "chain" << j+1 << endl;
		like=chain(numItera, factor_0, bestClusters);
		cout << "log likelihood =" << like <<", best cluster=" <<bestClusters.size() <<", cluster="<<allClusts.size()<<endl;
		if ((j==0) || (like > large))
		{
			finalClusters.clear();
			large=like;
			finalClusters.insert(finalClusters.end(),bestClusters.begin(),bestClusters.end());
			bestClusters.clear();
		}//end of if
	}//end of j
	bestClusters.clear();
	return(large);
}//end of bmc

double chain(const short NUM_CYCLE, const double FACTOR_0, list <clustClass> &bestClusters)
{
	short j;
	short cycle;
	short stage = 0;
	double like,like1,difference;
	list <clustClass>::iterator allClustsIter;
	list <motifClass>::iterator motifsIDIter;
	short currentClusterOrder;
//-	short del;
//-	vector <short> delVector;
	double diff =0;
	bool begin = true;
	double largest =0;

//*****
//  initialize
//*****
	cout <<"cluster= "<<allClusts.size()<<endl;
//*****
//  get initial likelihood
//*****
	like=loglikelihoodSingle();//0);
	cout <<"initial like="<<like<<endl;
//	double sum =0;
//	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//		sum = sum + (*allClustsIter).loglikelihood(0);
	//-cout << "original like=" << like <<endl;
//*****
//  start iteration
//*****
	stage = 0;///whether in annealing step
	begin = true;
	for(cycle = 0; cycle < NUM_CYCLE+ANNEALING;cycle ++ )
	{
		like=loglikelihoodSingle();//0);
		if(begin)
		{
			bestClusters.clear();
			largest = like;
			bestClusters.insert(bestClusters.end(),allClusts.begin(),allClusts.end());
			begin = false;
		}
//+		double sum =0;
//+		for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//+			sum = sum + (*allClustsIter).loglikelihood();//0);

//219		cout << "iteration=" << cycle+1 <<"\r";//"cluster=" << allClusts.size() << endl;
//***** special check on clusters
 		//-for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
		//-	(*allClustsIter).printClust(IN_WIDTH);		
		currentClusterOrder=0;
//		cout <<"cluster total = " << allClusts.size() <<endl;
		for(j=0;j<allMotifs.size();j++)// go through every motif
		{
//			cout <<"j="<<j;
//			if((cycle==0)&&(j==114))
//				cout<<"here.";
//219			cout <<"motif= "<<j<<" cluster size = "<<allClusts.size()<<endl;
//218			cout << "motif= " << j <<" " << "name=" <<allMotifs[j].name /*<< "clust=" << allMotifs[j].cluster*/ <<endl;
//			cout <<"shift= "<<allMotifs[23].motifshift<<endl;
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//			{
//				cout <<"width= "<<(*allClustsIter).width<<endl;
//				cout <<"remain= "<<(*allClustsIter).remain<<endl;
//218			}
			difference=operation(j, FACTOR_0);//,del);//,stage);
//03/03/07			like1=loglikelihoodSingle();//stage);
			like=like+difference;
//218			cout << "like=" <<like << "like1=" <<like1<<endl;
/* 03/03/07			if (fabs(like-like1)>0.00001)
			{
				cout << "problem here: motif = "<< j<<" dif=" << difference << "new like=" << like << " like1=" <<like1<<endl;
				exit(18);
			}
03/03/07 */
//***** special check on clusters
//			for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//				(*allClustsIter).printClust();
//added 02/15/07 cout <<"j="<<j<<" "<<allClusts.size()<<endl;
//***** end special check
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
//***** get the highest likelihood 02/21/06
		if(like > largest)
		{
			bestClusters.clear();
			largest = like;
			bestClusters.insert(bestClusters.end(),allClusts.begin(),allClusts.end());	
		}
	}//end of cycle
	//like1=loglikelihoodEven(clusters,WMIN);
	//cout << "like=" << like+difference << " like1=" <<like1<<endl;
	return(largest);//like+difference);changed 02/21/06	
}//end of chain

double operation(const short motifOrder, const double FACTOR_0)//, short &del)
				 //, const short stage)
{//for non shifted motifs only same width. 02/19/06
//	short j;
	short clustOrder = -1;
	vector <double> ratioVec;
	vector <list <clustClass>::iterator> clustsIterVec;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustPoint;
	short clusterNumber = 0;
	double ratio;
	short decision = -1, currentCluster = -1;//,shift = 0;
	clustClass tempClust;
	double likeminus,likeadd;
	bool empty = false;
	double postprob;

//	the pointer pointed to the cluster that the motif belongs to
	clustPoint = motifHome[motifOrder];
//  translate clustPoint to clustOrder 
	currentCluster = 0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		if (allClustsIter == clustPoint)
			break;
		else
			currentCluster ++;
	}//end of allClustsIter

	// remove motif from cluster 
	short a = (*clustPoint).motifsID.size();
	if((*clustPoint).motifsID.size() == 1) //for cluster with only 1 member, remove one exisitng cluster
	{
		empty = true;
//		cout <<"*****";
//		likeminus = - (*clustPoint).loglikelihood();//0);
//		allClusts.erase(clustPoint);
	}//end of if
	else
	{
	//	likeminus=clustRemoveMotif(motifOrder, clustPoint);//, stage);
		empty = false;
	}
//start trying to fit all clusters
	short count =0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		count ++;
		if(allClustsIter == clustPoint)
		{
			if(empty)
				ratioVec.push_back(0);
			else
			{
				ratio=exp(bayesMinusRatioSingle(motifOrder, allClustsIter));//,motifStart,clustStart,stage));
				ratioVec.push_back(ratio);
			}
		}
		else
		{
			ratio=exp(bayesRatioSingle(motifOrder, allClustsIter));//,motifStart,clustStart,stage));
			ratioVec.push_back(ratio);
		}
		clustsIterVec.push_back(allClustsIter);
	}//end of allClustiter
// decide which cluster fit the motif best
	decision=decide(ratioVec, FACTOR_0,postprob);
	allMotifs[motifOrder].postprob = postprob;
//	cout << "decision=" << decision;
// if the current cluster is the best, then do nothing
	if(empty)
	{
		if(decision ==0)
		{
//			cout <<"cluster size= "<<allClusts.size()<<endl;
			return(0);
		}
		else
		{
			likeminus = -(*clustPoint).loglikelihood();
//***** changed 03/03/07				
//+			likeadd=clustInsertMotif(motifOrder,decision-1);
			likeadd=clustInsertMotif(motifOrder,clustsIterVec[decision-1]);
			allClusts.erase(clustPoint);
//			cout <<"cluster size= "<<allClusts.size()<<endl;
			return(likeminus-likeadd);
		}
	}//end of if
	else if(decision==(currentCluster + 1))//added !empty 08/07/05
	{// same cluster
//		cout <<"cluster size= "<<allClusts.size()<<endl;
//		(*clustPoint).addMotif( motifOrder);
		(*clustPoint).motifsID.pop_front();
		(*clustPoint).motifsID.push_back(motifOrder);
		return(0);
	}//end of if
// if a different cluster is desired
	else 
	{
// remove it from original cluster
		//(*clustPoint).motifs.pop_front();
		likeminus=clustRemoveMotif(motifOrder, clustPoint);
		if(decision > 0)
		{
//***** changed 03/03/07			
//+			likeadd=clustInsertMotif(motifOrder,decision-1);//, stage);// moved down from before if 08/08/05
			likeadd=clustInsertMotif(motifOrder,clustsIterVec[decision-1]);
		}//end of if
// if a new cluster is desired
		else if (decision == 0)
		{
			tempClust.motifsID.push_back(motifOrder);
			//tempClust.width = allMotifs[motifOrder].width;
			tempClust.iniProfile();//bin);
			allClusts.push_back(tempClust);
			tempClust.motifsID.clear();
			clustPoint = allClusts.end();
			clustPoint --;
			motifHome[motifOrder] =clustPoint;
			likeadd = allMotifs[motifOrder].ratio;
		}//end of else if		
//		cout <<"cluster size= "<<allClusts.size()<<endl;	
		return(likeminus-likeadd);
	}//end of else
}//end of operation

/*
double operation(const short motifOrder, const double FACTOR_0)//, short &del)
				 //, const short stage)
{//for non shifted motifs only same width. 02/19/06
//	short j;
	short clustOrder = -1;
	vector <double> ratioVec;
//	vector < vector <double> > ratioAll;
	list <clustClass>::iterator allClustsIter;
	list <clustClass>::iterator clustPoint;
	short clusterNumber = 0;
	double ratio;
	short decision = -1, currentCluster = -1;//,shift = 0;
	clustClass tempClust;
	double likeminus,likeadd;
	bool empty = false;
	double postprob;
//-	vector <bool> clustlong;
//-	double likeold, likenew;

//	the pointer pointed to the cluster that the motif belongs to
	clustPoint = motifHome[motifOrder];
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
//-	del=-1;
//218	cout << "start of operation, motif = " <<motifOrder<<"size = " <<(*clustPoint).motifsID.size() <<endl;
	// remove motif from cluster 
	short a = (*clustPoint).motifsID.size();
	if((*clustPoint).motifsID.size() == 1) //for cluster with only 1 member, remove one exisitng cluster
	{
		empty = true;
		//likeminus=loglikelihoodUpdateSingle(motifOrder, clustPoint,stage);	
		likeminus = - (*clustPoint).loglikelihood();//0);
		allClusts.erase(clustPoint);
//218		cout <<"empty likeminus = "<<likeminus<<endl;
	}//end of if
	else
		likeminus=clustRemoveMotif(motifOrder, clustPoint);//, stage);
//start trying to fit all clusters
//-	short tmp = allClusts.size();
	short count =0;
	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
	{
		count ++;
// decide shift
//-		clustWidth = (*allClustsIter).width;
//-		motifWidth = allMotifs[motifOrder].width;
//-		clustlong.push_back(true);
//-		shift = 1;
//-		motifStart = 0;
//-		clustStart = 0;
		ratio=exp(bayesRatioSingle(motifOrder, allClustsIter));//,motifStart,clustStart,stage));
		ratioVec.push_back(ratio);
	}//end of allClustiter
// decide which cluster fit the motif best
	decision=decide(ratioVec, FACTOR_0,postprob);
	allMotifs[motifOrder].postprob = postprob;
//	cout << "decision=" << decision <<endl;
// if the current cluster is the best, then do nothing
	if((!empty)&&(decision==(currentCluster + 1)))//added !empty 08/07/05
	{// same cluster
		//(*clustPoint).motifs.pop_front();

		//(*clustPoint).updateProfile(1,allMotifs, motifOrder);
		//(*clustPoint).motifs.push_back(allMotifs[motifOrder]);
//-		clustlong.clear();
//-		if((oldmotifshift == allMotifs[motifOrder].motifshift)&&(oldclustshift == allMotifs[motifOrder].clustshift))
//-		{
			(*clustPoint).addMotif( motifOrder);
			return(0);
//-		}
//*****************
//		else 
//		{
//			likeold = (*clustPoint).loglikelihood(0);
//			(*clustPoint).addMotif( motifOrder);
//			likenew = (*clustPoint).loglikelihood(0);
//            likeadd = likeold - likenew;
// //218       cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
//            return(likeminus - likeadd) ;
//		}
//*********************
	}//end of if
// if a different cluster is desired
	else 
	{
// remove it from original cluster
		//(*clustPoint).motifs.pop_front();
		if(decision > 0)
		{
// insert it into the desired cluster
			//likeadd=clustInsertMotif(motifOrder,decision-1,clustWidth);//, stage);
			//if(MAX_SHIFT>1)
			//allMotifs[motifOrder].motifshift = shift;//update shift, added 08/03/05 //remove -MAX_SHIFT/2 08/06/05
//*************************
//			if(clustlong[decision-1])
//			{
//				allMotifs[motifOrder].motifshift = 0;
//				allMotifs[motifOrder].clustshift = shift;
//			}
//			else
//			{
//				allMotifs[motifOrder].motifshift = shift;
//				allMotifs[motifOrder].clustshift = 0;
//			}
//*****************************
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
//-			allMotifs[motifOrder].motifshift = 0;
//-			allMotifs[motifOrder].clustshift = 0;
//			tempMotifs.push_back(allMotifs[motifOrder]);//(*motifsIDIter);
			tempClust.motifsID.push_back(motifOrder);
			//tempClust.width = allMotifs[motifOrder].width;
			tempClust.iniProfile();//bin);
			allClusts.push_back(tempClust);
			tempClust.motifsID.clear();
			clustPoint = allClusts.end();
			clustPoint --;
			motifHome[motifOrder] =clustPoint;
//			tempMotifs.clear();
			//likeadd=loglikelihoodNewEven(clustWidth,motifsIDIter,stage);//***)
			if (empty)
				likeadd=likeminus;
			else
			{
				likeadd = allMotifs[motifOrder].ratio;
//-	deleted 02/15/07	likeadd=emptyLoglikelihoodupdateSingle(motifOrder);//,clustPoint);//,stage);
			}//end of else
		//-	cout << " likeadd2=" <<likeadd<<endl;
		}//end of else if		
// remove it from its current cluster
//-		del=motifOrder;//.first
//-		clustlong.clear();
//218		cout <<"likeminus= " <<likeminus <<"likeadd= "<<likeadd<<endl;
		return(likeminus-likeadd);
	}//end of else
}//end of operation
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

	//***** seed random number generator 
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
//-		shift=allMotifs[(*motifsIDIter)].motifshift;
		for(m=0;m<4;m++)
		{
//-			la = allMotifs[(*motifsIDIter)].profile[m][position+shift];
			la = allProfile[m][(*motifsIDIter)*width + position];
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

double onoffMetro(const list <clustClass>::iterator clustIter, short W, short WMIN, short bin)
{
	int na,nb,nc;
	double dif,ran,lower;
	short current, leftMargin,rghtMargin;

	//***** seed random number generator 
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

double clustRemoveMotif(const short motifOrder, const list <clustClass>::iterator clustPoint)//, const short stage)
{
//	short oldwidth = 0;
//	short newwidth = 0;
	double likenew,likeold;

//	oldwidth = (*clustPoint).width;
	likeold = (*clustPoint).loglikelihood();//0);
	(*clustPoint).minusMotif(motifOrder);
	likenew = (*clustPoint).loglikelihood();//0);
	return(likenew-likeold);
	cout <<"remove: old= "<<likeold<<"new= "<<likenew<<endl;
}// end of clustRemoveMotif

double clustInsertMotif(const short motifOrder, const list <clustClass>::iterator clustIter)//, const short stage)
{
//+	short count;
//+	list <clustClass>::iterator clustIter;
	double likenew,likeold;
//-	short clustWidth = 0;
//-	short newwidth = 0;
			
	likeold = (*clustIter).loglikelihood();//0);
	(*clustIter).addMotif( motifOrder);
	likenew = (*clustIter).loglikelihood();//0);
// assign this cluster to motifHome
	motifHome[motifOrder] = clustIter;
//218			cout <<"insert: old= "<<likeold<<"new= "<<likenew<<endl;
	return(likeold-likenew);		
}// end of clustInsertMotif

/*
double clustInsertMotif(const short motifOrder, const short decision)//, const short stage)
{
	short count;
	list <clustClass>::iterator clustIter;
	double likenew,likeold;
//-	short clustWidth = 0;
//-	short newwidth = 0;
			
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

/*
double emptyLoglikelihoodupdateSingle(const short motifOrder)//, const list <clustClass>::iterator allClustsIter)
										//, const short stage)
{
	short ** newProfile;
	short j,k;
	double logsum;//,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
//	short lsummotif;
	double difference;//,tem;
//-	short start;
//-	short width;

//-	width = (*allClustsIter).width;
	newProfile = new short*[4];
	for(j=0;j<4;j++)
		newProfile[j] = new short[fixedwidth];
//-	start= allMotifs[motifOrder].motifshift;
//	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	//lsum=0.0;
	motifsum=0;
//	lsummotif=allMotifs[motifOrder].depth;
	lsumnew=allMotifs[motifOrder].depth;//lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
	for (j=0;j<fixedwidth;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+start];
				newProfile[k][j] = allProfile[k][motifOrder*fixedwidth + j];//allMotifs[motifOrder].profile[k][j+start];//***shift to start//remove (*allClustsIter).profile[k][j]
			}//end of k
//-		}//end of if
	}//end of j
	for(j=0;j<fixedwidth;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
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
//-		}//end of if
/*
		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb = allProfile[k][motifOrder*width + j];//allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		}//end of else if
*//*
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	for(j=0;j<4;j++)
		delete [] newProfile[j];
	delete [] newProfile;
//	double x = allMotifs[motifOrder].ratio;
	return (-difference);
}//end of emptyLoglikelihoodupdateSingle
*/

//-----
// calculate likelihood for clusters of single motifs
//-----
double loglikelihoodSingle()//const short stage)
{
	short j,k;
//	double tem;
	short lsum;
	list <clustClass>::iterator allClustsIter;
//-	list <short>::iterator motifsIDIter;
	double logsum,la;
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
		for(j=0;j<fixedwidth;j++)
		{
//-			if((stage!=1)||((*allClustsIter).power[j]))
//-			{
				for(k=0;k<4;k++)
				{
					//+short a=(*allClustsIter).profile[k][j];
					la = (double) (*allClustsIter).profile[k][j];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
//03/03/07			logsum=logsum;// -gammaln((double) lsum + betasum);
//-			}//end of if
/*	
			else if ((stage==1)&&(!(*allClustsIter).power[j]))
			{
				for(motifsIDIter = (*allClustsIter).motifsID.begin();motifsIDIter != (*allClustsIter).motifsID.end();motifsIDIter ++)
				{
					sumin = allMotifs[(*motifsIDIter)].depth;
//-					start = allMotifs[(*motifsIDIter)].motifshift;// add (*motifsIDIter).shift 
					for(k=0;k<4;k++)
					{
//-						la = (double) allMotifs[(*motifsIDIter)].profile[k][j+start];//+beta[k];
						la = (double) allProfile[k][(*motifsIDIter)*width + j];
						logsum=logsum+gammaln(la+beta[k]);
					}//end of k
					logsum=logsum-gammaln((double) sumin + betasum);
				}//end of 
				logsum=logsum+(double) (motifsSize-1)*tem;
			}//end of else if
*/
		}//end of j
		logsum = logsum + (double) fixedwidth * (betasumdif - gammaln((double) lsum + betasum));//tem;
//218		cout <<"cluster like="<<logsum-logsumold<<endl;
//+        logsumold = logsum;
	}//end of allClustsIter
//218	cout <<endl;
	//logsum=logsum+allClusts.size()*(double) width *tem;
	return(logsum);
}//end of loglikelihoodSingle

double bayesRatioSingle(short motifOrder, list <clustClass>::iterator allClustsIter)
						//,short motifStart, short clustStart, const short stage)
{
	short j,k;
	short Lold,Lnew,sum;
	double sumbayes,tem;
	double newProfile;
	double bayesRatio;
//	short a,b;
//	short width;
        int startpos;

        startpos = motifOrder*fixedwidth;
//-	width = MIN((*allClustsIter).width,allMotifs[motifOrder].width);//changed 9/13/05
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum+(*allClustsIter).profile[j][0];
	}
	Lold=sum;
//	Lz=allMotifs[motifOrder].depth;
	Lnew=Lold+allMotifs[motifOrder].depth;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<fixedwidth;j++)
	{
//-		if((stage !=1)||((*allClustsIter).power[j]))//removed ==1 08/10/05
//-		{
			for(k=0;k<4;k++)
			{
				//a=(*allClustsIter).profile[k][j+clustStart];
				//b=allMotifs[motifOrder].profile[k][j+motifStart];
				newProfile=(*allClustsIter).profile[k][j] + allProfile[k][startpos + j];//[(j) + clustStart]//03/03/07 "startpos" replace "motifOrder*fixedwidth".
					//allMotifs[motifOrder].profile[k][j+motifStart];
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) (*allClustsIter).profile[k][j]+beta[k]);//[(j) + clustStart]
				//	-gammaln((double) allProfile[k][startpos + j] + beta[k]);//03/03/07 "startpos" replace "motifOrder*fixedwidth".
					//-gammaln((double) allMotifs[motifOrder].profile[k][j+motifStart]+beta[k]);
			}//end of k
//-		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum);
//	+gammaln((double) Lz+betasum)-gammaln(betasum)+sumbeta;
//-	if(stage == 1)
//-		bayesRatio = sumbayes+tem*(double) (*allClustsIter).remain;
//-	else 
//+	double x = allMotifs[motifOrder].ratio;
	bayesRatio = sumbayes+tem*(double) fixedwidth + allMotifs[motifOrder].ratio;
	return(bayesRatio);
}//end of bayesRatioSingle

double bayesMinusRatioSingle(int motifOrder, list <clustClass>::iterator allClustsIter)
{
	int j,k;
	int Lold,Lz,Lnew,sum;
	double sumbayes,tem;
	int oldProfile;
	double bayesRatio;
        int startpos;

        startpos = motifOrder*fixedwidth;	
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum + (*allClustsIter).profile[j][0];
	}
	Lnew=sum;
	Lz=allMotifs[motifOrder].depth;
	Lold=Lnew - Lz;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<fixedwidth;j++)
	{
		for(k=0;k<4;k++)
		{
			//a=(*allClustsIter).profile[k][j+clustStart];
			//b=allMotifs[motifOrder].profile[k][j+motifStart];
			oldProfile=(*allClustsIter).profile[k][j] - allProfile[k][startpos + j];//[(j) + clustStart]//03/03/07 "startpos" replace "motifOrder*fixedwidth".
				//allMotifs[motifOrder].profile[k][j+motifStart];
			sumbayes=sumbayes+gammaln((double) (*allClustsIter).profile[k][j] + beta[k])
				-gammaln((double) oldProfile + beta[k]);//[(j) + clustStart]
			//	-gammaln((double) allProfile[k][startpos + j] + beta[k]);//03/03/07 "startpos" replace "motifOrder*fixedwidth".
		}//end of k
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum);	
	//double x = allMotifs[motifOrder].ratio;
	bayesRatio = sumbayes+tem*(double) fixedwidth + allMotifs[motifOrder].ratio;
	return(bayesRatio);
}//end of bayesMinusRatioSingle

double loglikelihoodUpdateSingle (const short motifOrder, const list <clustClass>::iterator allClustsIter)
									//, const short stage)
{
//+	short ** newProfile;
	short j,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
//+	double motifsum;
	short lsummotif;
	double difference;//,tem;
//-	short width;
//-	short start;
        int startpos;

        startpos = motifOrder*fixedwidth;
//-	width = (*allClustsIter).width;
	// need to find out if the cluster width changed during the add/removal of the motif
//+	newProfile = new short*[4];
//+	for(j=0;j<4;j++)
//+		newProfile[j] = new short[fixedwidth];
//	start = allMotifs[motifOrder].motifshift;
//	tem=gammaln(betasum)-sumbeta;
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
	for (j=0;j<fixedwidth;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for (k=0;k<4;k++)
			{
				//a=allMotifs[motifOrder].profile[k][j+start];
				newProfile[k][j]=(*allClustsIter).profile[k][j]
							 + allProfile[k][motifOrder*fixedwidth + j];
//-							 +allMotifs[motifOrder].profile[k][j+start];//***shift to start
			}//end of k
//-		}//end of if
	}//end of j
+ */
	for(j=0;j<fixedwidth;j++)
	{
//-		if((stage!=1)||((*allClustsIter).power[j]))
//-		{
			for(k=0;k<4;k++)
			{
				la=(*allClustsIter).profile[k][j];
				logsum=logsum+gammaln((double)la+beta[k]);
//				lanew=newProfile[k][j];
				lanew = la + allProfile[k][startpos + j];//03/03/07 "startpos" replace "motifOrder*fixedwidth".
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum;//-gammaln(lsum+betasum);
			logsumnew=logsumnew;//-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
//-		}//end of if
/*		
		else if((stage==1)&&(!(*allClustsIter).power[j]))
		{
			for(k=0;k<4;k++)
			{
				lb=allMotifs[motifOrder].profile[k][j+start];//***shift to start//allMotifs[motifOrder].start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
//-		}//end of else if
*/
	}//end of j	
//+	difference=logsumnew-logsum+motifsum;
	difference=logsumnew-logsum -fixedwidth*(gammaln(lsumnew+betasum) - gammaln(lsum+betasum));
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
//+	for(j=0;j<4;j++)
//+		delete [] newProfile[j];
//+	delete [] newProfile;
	return difference;
}//end of loglikelihoodupdateSingle

//----- read in motif block data and put into map and list -----//
short readSeq(const char inputFileName[])//, short &widthmin)
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
//-			temMotif.motifshift = 0;
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
	fixedwidth = allMotifs[0].width;

	allProfile = new int*[4];
	for(j=0;j<4;j++)
	{
		allProfile[j] = new int[totalcol];
		for(k=0;k<totalcol;k++)
			allProfile[j][k]=0;
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
			genProfile(seq,start,count);//width,allMotifs[count].depth);
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
	if(minwidth != maxwidth)
//	if(minwidth == maxwidth)
//		fixedwidth = minwidth;
//	else
	{
		cout <<"motif widths are different, check."<<endl;
		exit(99);
	}
	return(motifCount);
}//end of readSeq

void genProfile(vector <string> seq, const int start, const int motifOrder)
//const int width, const int depth)
{
	int j,k;
	int depth, Lz;
	double sum, tem, ratio;
	
	depth = allMotifs[motifOrder].depth;
// get count
	for(j=0;j<depth;j++)
	{
		for(k=0;k<fixedwidth;k++)
		{
			switch (seq[j][k])
			{
			case 'a': allProfile[0][start + k]++; break;
			case 'c': allProfile[1][start + k]++; break;
			case 'g': allProfile[2][start + k]++; break;
			case 't': allProfile[3][start + k]++; break;
			case 'A': allProfile[0][start + k]++; break;
			case 'C': allProfile[1][start + k]++; break;
			case 'G': allProfile[2][start + k]++; break;
			case 'T': allProfile[3][start + k]++; break;
			default: cout << "error in motif profile" <<endl;
			}//end of switch
		}//end of k
	}//end of j
/*
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << allProfile[j][k + start] << " ";
		cout << endl;
	}//end of j
	cout << endl;
*/
	Lz = allMotifs[motifOrder].depth;
	sum = 0;
	for(j=0;j<fixedwidth;j++)
	{
		for(k=0;k<4;k++)
			sum = sum + gammaln((double) allProfile[k][start + j] + beta[k]);//03/02/07 "start" replaced "motifOrder*fixedwidth".
	}
	tem = gammaln((double) Lz+betasum) - betasumdif;//-gammaln(betasum)+sumbeta;
	ratio = -sum + tem*(double)fixedwidth;
	allMotifs[motifOrder].ratio = ratio;
}//end of genProfile

void initialCluster(/*const short motifCount,*/ const short INITIAL_CLUSTER)
{
	short j;
	list <motifClass> clust;
	list < clustClass >::iterator allClustsIter;
	list < motifClass >::iterator motifsIDIter;
	clustClass tempClust;
	short motifsSize;
//	short *width, group;

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
	for(j=0;j<INITIAL_CLUSTER;j++)// initialize clusters
	{
//-		allMotifs[j].clustshift =0;
//+		cout <<"here.*****"<<endl;
		tempClust.motifsID.push_back(j);
//		tempClust.width = width[j];
//-		tempClust.width = allMotifs[j].width;
		tempClust.iniProfile();//bin);
//+		short b = tempClust.profile[0][0];
//+		cout <<"*****b="<<b<<endl;//++
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
//-		allMotifs[j].clustshift =0;
		//group = j % INITIAL_CLUSTER;
		if ((j % INITIAL_CLUSTER) == 0)// set allClusts pointer to the beginning
		{
			allClustsIter = allClusts.begin();
//-			short a = (*allClustsIter).profile[0][0];
			(*allClustsIter).addMotif( j);
			motifHome.push_back(allClustsIter);
		}
		else // add motifs to existing clusters
		{
			allClustsIter ++;
//-			allMotifs[j].clustshift =0;		
			(*allClustsIter).addMotif(j);
			motifHome.push_back(allClustsIter);
		}
	}//end of if
//	delete []width;
//219	for(allClustsIter = allClusts.begin();allClustsIter != allClusts.end(); allClustsIter++)
//219		(*allClustsIter).printClust();
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

