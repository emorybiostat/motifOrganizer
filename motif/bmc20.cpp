// clustering program BMC 2.0
// origin from bmc by Steve Qin, 01/15/02
// modified on 03/16/04 using STL list and map
//
//
#include <map>
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
#define INITIAL_CLUSTER 2
#define IN_WIDTH 16//16//2
#define MIN(a,b) ((a) < (b) ? (a):(b))
#define MAX(a,b) ((a) > (b) ? (a):(b))

using namespace std;

double gammaln(double xx);
static double beta[4]={1.0,1.0,1.0,1.0};
static double betasum=beta[0]+beta[1]+beta[2]+beta[3];
static double sumbeta=gammaln(beta[0])+gammaln(beta[1])+gammaln(beta[2])+gammaln(beta[3]);
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
	int start;
//	int shift;
	vector <string> basepair;
	vector <int> profile[4];
//public:
	void printSeq();
	void genProfile();
	//void iniProfile();
//	friend void readSeq(list <clustClass > &clustersList);
};

class clustClass
{
private:
public:
	map <int,motifClass> cluster;
	vector <bool> power;
	double remain;
	//int left;
	//int right;
	//int width;
	vector <int> profile[4];
	void printClust(int clustWidth, ofstream &outputFile);
	//void genProfile(int wid);
	void clustClass::iniProfile(int clustWidth);
	void clustClass::addMotif(motifClass &motif, int clustWidth);
	int clustClass::minusMotif(int order, int clustWidth);
	void clustClass::updateProfile(int sign, motifClass &motif, int clustWidth);
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
			default: cout << "error in motif profile" <<endl;
			}//end of switch
		}//end of k
	}//end of j
	for(j=0;j<4;j++)
	{
		for(k=0;k<width;k++)
			cout << profile[j][k] << " ";
		cout << endl;
	}//end of j
	cout << endl;
}//end of motifClass::genProfile

/*
void clustClass::genProfile(int clustWidth)
{
	int j,k;
	int start;
	map <int, motifClass>::iterator motifIter;
	
	for(motifIter = cluster.begin(); motifIter != cluster.end(); motifIter++)
	{
		start=((*motifIter).second.width-clustWidth)/2;
		for(j=0;j<(*motifIter).second.depth;j++)
		{
			for(k=0;k<clustWidth;k++)
			{
				switch ((*motifIter).second.basepair[j][k+start])
				{
					case 'a': profile[0][k]++; break;
					case 'c': profile[1][k]++; break;
					case 'g': profile[2][k]++; break;
					case 't': profile[3][k]++; break;
					case 'A': profile[0][k]++; break;
					case 'C': profile[1][k]++; break;
					case 'G': profile[2][k]++; break;
					case 'T': profile[3][k]++; break;
					default: cout << "error in motif profile" <<endl;
				}//end of switch
			}//end of k
		}//end of j
	}//end of motifIter
	for(j=0;j<4;j++)
	{
		for(k=0;k<clustWidth;k++)
			cout << profile[j][k];
		cout << endl;
	}//end of j
}//end of clustClass::genProfile()
*/ 
void clustClass::printClust(int clustWidth, ofstream &outputFile)
{
//	int j,k;
	map <int, motifClass>::iterator motifIter;

	//clustWidth=profile[0].size();
	for(motifIter = cluster.begin(); motifIter != cluster.end(); motifIter++)
	{
		outputFile << "ID " << (*motifIter).first;
		outputFile << "Name " << (*motifIter).second.name << (*motifIter).second.width << 
			(*motifIter).second.depth << (*motifIter).second.isEven << endl;
		//(*motifIter).second.printSeq();
		cout << "ID " << (*motifIter).first;
		cout << "Name " << (*motifIter).second.name << (*motifIter).second.width << 
			(*motifIter).second.depth << (*motifIter).second.isEven << endl;
	}//end of motifIter
/*	for(j=0;j<4;j++)
	{
		for(k=0;k<clustWidth;k++)
		{
			outputFile << profile[j][k] <<" ";
			cout << profile[j][k] << " ";
		}
		cout << endl;
		outputFile << endl;
	}//end of j
	outputFile << endl;
*/
	outputFile.close();
}//end of clustClass::printClust()

//----- initialize the profile matrix when the first motif join this cluster

void clustClass::iniProfile(int clustWidth)
{
	int j,k;
	int start;
	map <int, motifClass>::iterator motifIter;
	int motifWidth;

	motifIter=cluster.begin();
	start=(*motifIter).second.start;
	motifWidth=(*motifIter).second.width;
// first set profile matrix to 0
	for(j=0;j<clustWidth;j++)
	{
		power.push_back(1);
		for(k=0;k<4;k++)
			profile[k].push_back(0);
	}
	for(j=0;j<clustWidth;j++)
	{
		power[j]=1;
		for(k=0;k<4;k++)
			profile[k][j]=(*motifIter).second.profile[k][j+start];//+(*motifIter).second.profile[kstar[k]][motifWidth-1-j-start];
	}
}//end of clustClass::iniProfile

/*
void clustClass::genProfile(int wid)
{
	int j,k;
	int offset;
	map <int, motifClass>::iterator motifIter;
	vector <int> profile[4];

// initialize
	for(j=0;j<wid;j++)
	{
		for(k=0;k<4;k++)
			profile[k].push_back(0);
	}
	for(motifIter = cluster.begin(); motifIter != cluster.end(); motifIter++)
	{
		offset=(*motifIter).second.width-width;
		for(j=0;j<wid;j++)
		{
			for(k=0;k<4;k++)
				profile[j][k]=profile[j][k]+(*motifIter).second.profile[j][k+offset];
		}	
	}
}//end of clustClass::genProfile()
*/

void clustClass::addMotif(motifClass &motif, int clustWidth)
{
	cluster[motif.id]=motif;
	updateProfile(1,motif,clustWidth);
}

int clustClass::minusMotif(int order, int clustWidth)
{
	if(cluster.size()>1)// there are ore han one element in the cluster
	{
		cluster.erase(order);
		updateProfile(-1,cluster[order],clustWidth);
		return(0);
	}
	else if (cluster.size() == 1)//there is only one element in this cluster
		return(1);
}

void clustClass::updateProfile(int sign, motifClass &motif, int clustWidth)
{
	int j,k;
	int start;

	start=motif.start;
/*
	for(j=0;j<4;j++)
	{
		for(k=0;k<clustWidth;k++)
			cout << profile[j][k] << " ";
		cout <<endl;
	}
	cout <<endl;

	for(j=0;j<4;j++)
	{
		for(k=0;k<clustWidth;k++)
			cout << motif.profile[j][k+start] << " ";
		cout <<endl;
	}
	cout <<endl;
*/
	for(j=0;j<4;j++)
	{
		for(k=0;k<clustWidth;k++)
		profile[j][k]=profile[j][k]+sign * motif.profile[j][k+start];
		//	profile[j][k]=profile[j][k]+sign * (motif.profile[j][k+start]+
		//				motif.profile[kstar[j]][motif.width-1-k-start]);
	}//end of j
}

class clustClass;

class motifClass;

int readSeq(int clustWidth, list <clustClass> &clustersList);

void bmc(const int numChains, const int numItera, const int W, const int WMIN, const int WMAX, const int MAX_SHIFT, 
		 const double factor_0, list <clustClass> &clusters, list <clustClass> &finalClusters);

double loglikelihoodEven(list <clustClass> clusters, const int width);

double chain(const int NUM_CYCLE, const int W, const int WMIN, const int WMAX, const int MAX_SHIFT, 
			 const double FACTOR_0, list <clustClass> &clusters);

double operation(const double FACTOR_0, const int clustWidth, const map <int, motifClass>::iterator motifIter, 
			     const list <clustClass>::iterator clustIter, const int order, list <clustClass> &clusters,
				 int &del);

double clustInsertMotif(list <clustClass> &clusters, const map <int, motifClass>::iterator motifIter,
					  const int decision, const int clustWidth);

double loglikelihoodUpdateEven(int width, const list <clustClass>::iterator clustIter,
						   const map <int, motifClass>::iterator motifIter);

double loglikelihoodNewEven(int width, /*const list <clustClass>::iterator clustIter,*/
						   const map <int, motifClass>::iterator motifIter);

double bayesRatio(motifClass motif, list <clustClass>::iterator clustIter,int W,int start);

int decide(vector <double> &ratioVec);

void output(int clustWidth, list <clustClass> clusters, char outputFilename[]);

double unran(int *na, int *nb, int *nc);

//void output(const list <clustClass> *clustersList);

double special(int width, const list <clustClass>::iterator clustIter,
						   const map <int, motifClass>::iterator motifIter);

int main()
{
	list <clustClass> clusters;
	list <clustClass> finalClusters;
	int W,WMIN,WMAX;
	int numChains,numIter,MAX_SHIFT;
	double factor_0;
	int motifCount;

	W=16;//16;//2;
	WMIN=16;//16;//2;
	WMAX=16;//16;//6;
	MAX_SHIFT=1;
	numChains=1;
	numIter=10;
	factor_0=1.0;

// read in data, initial assign clusters and calculate cluster profiles
	motifCount=readSeq(W,clusters);
// print out information
	//printClust(clusters);
// clustering
	cout << clusters.size();
	bmc(numChains,numIter,W, WMIN,WMAX,MAX_SHIFT,factor_0,clusters,finalClusters);
// output result
//	output(clusters);
	return(0);
}//end of main

void bmc(const int numChains, const int numItera, const int W, const int WMIN, const int WMAX, 
		 const int MAX_SHIFT, const double factor_0, list <clustClass> &clusters, list <clustClass> &finalClusters)
{
	int j;
	double like,large;

	for(j=0;j<numChains;j++)
	{
		cout << "chian" << j << endl;
		like=chain(numItera, W, WMIN,WMAX, MAX_SHIFT, factor_0,clusters);
		cout << "log likelihood =" << like <<endl;
		if ((j==0) || (like > large))
		{
			large=like;
			finalClusters=clusters;
		}//end of if
	}//end of j
}//end of bmc

double chain(const int NUM_CYCLE, const int W, const int WMIN, const int WMAX, const int MAX_SHIFT, 
			 const double FACTOR_0,list <clustClass> &clusters)
{
	int j,k;
	int cycle;
	bool stage;
	double like,like1,difference;
	list <clustClass>::iterator clustIter;
	map <int, motifClass>::iterator motifIter;
	int currentClusterOrder;
	int del;
	vector <int> delVector;
	int a;

//*****
//  initialize
//*****

//*****
//  get initial likelihood
//*****
	like=loglikelihoodEven(clusters,WMIN);
	cout << "original like=" << like <<endl;
//*****
//  start iteration
//*****
	stage=0;//whether in annealing step
	for(cycle = 0; cycle < NUM_CYCLE+ANNEALING;cycle ++ )
	{
		cout << "begin cycle=" << cycle <<"cluster=" << clusters.size() << endl;
		if(cycle == (NUM_CYCLE/3*2))
		{
			stage=1;
		}
		currentClusterOrder=0;
		cout <<"cluster total = " << clusters.size() <<endl;
		//for(clustIter = (clusters).begin();clustIter != (clusters).end();clustIter ++)
		clustIter=(clusters).begin();
		while (clustIter !=(clusters).end())
		{
		//	cout <<"motif total = " << (*clustIter).cluster.size() <<endl;
			a=clusters.size();
			for(motifIter = (*clustIter).cluster.begin();motifIter != (*clustIter).cluster.end();motifIter ++)
			{
				cout << "motif=" << (*motifIter).second.id <<" " << "clust=" <<currentClusterOrder <<endl;
				difference=operation(FACTOR_0, W,motifIter, clustIter, currentClusterOrder, clusters,del);
				like1=loglikelihoodEven(clusters,WMIN);
				like=like+difference;
//				output(W,clusters,"clust.txt");
				cout << "dif=" << difference << "new like=" << like << " like1=" <<like1<<endl;
				//exit(0);
				if(del !=-1)
					delVector.push_back(del);
			}
// remove it from its current cluster
			k=(*clustIter).cluster.size();
			k=delVector.size();
			if(delVector.size()==0)
			{
				++clustIter;
				currentClusterOrder++;
			}//end of if
			else
			{
				if((*clustIter).cluster.size()==delVector.size())//need to delete whole cluster
				{
					//(*clustIter).cluster.empty();
					//(*clustIter).cluster.clear();
					clustIter=clusters.erase(clustIter);
				}//end of if
				else if((*clustIter).cluster.size()>delVector.size())
				{
					for(j=0;j<delVector.size();j++)
					{
						(*clustIter).cluster.erase(delVector[j]);
					}
					++clustIter;
					currentClusterOrder++;
				}//end of else if
				else 
				{
					cout << "delete error" <<endl;
					exit(0);
				}//end of else
			}//end of else
			//output(W,clusters,"clust.txt");
			delVector.clear();
		}// end of clustIter

	}//end of cycle
	//like1=loglikelihoodEven(clusters,WMIN);
	//cout << "like=" << like+difference << " like1=" <<like1<<endl;
	return(like+difference);	
}//end of chain

void output(int clustWidth, list <clustClass> clusters, char outputFilename[])
{
	list <clustClass>::iterator clustIter;
	ofstream outputFile;
	int count;

	outputFile.open(outputFilename);
	if(!outputFile)
	{
		cout << "ERROR: Unable to open file: " << outputFilename << endl;
	    exit(2);
	}//end of if
	count=1;
	for(clustIter = clusters.begin();clustIter != clusters.end();clustIter ++)
	{
		cout << "cluster " << count << endl; 		
		count++;
		(*clustIter).printClust(clustWidth, outputFile);
	}
	outputFile.close();
}//end of output

double operation(const double FACTOR_0, const int clustWidth, const map <int, motifClass>::iterator motifIter, 
			     const list <clustClass>::iterator clustIter, const int order, list <clustClass> &clusters,
				 int &del)
{
	list <clustClass>::iterator clustIter0;
	vector <double> ratioVec;
	double ratio;
	int decision;
	int start;
	//motifClass temMotif;
	map <int,motifClass> tempClustMap; 
	clustClass tempClust;
	double likeminus,likeadd;
	int j,k;

	del=-1;
	// temporary remove motif from cluster 
/*	for(k=0;k<4;k++)
			{
				a=(*clustIter).profile[k][0];
				b=(*clustIter).profile[kstar[k]][1];
				}
	cout <<(*clustIter).cluster.size() <<endl;; 
*/
	(*clustIter).updateProfile(-1,(*motifIter).second, clustWidth);
/*
	for(k=0;k<4;k++)
			{
				a=(*clustIter).profile[k][0];
				b=(*clustIter).profile[kstar[k]][1];
				}
*/
// go over all current clusters
	if((*clustIter).cluster.size()>1)
		ratioVec.push_back(FACTOR_0);//represents new cluster
	else
		ratioVec.push_back(0);//for cluster with only 1 member, do not allow new cluster generating
	for(clustIter0 = (clusters).begin();clustIter0 != (clusters).end();clustIter0 ++)
	{
		start=0;//change
		ratio=exp(bayesRatio((*motifIter).second,clustIter0,clustWidth,start));
		ratioVec.push_back(ratio);
	}
// decide which cluster fit the motif best
	decision=decide(ratioVec);
	cout << "decision=" << decision <<endl;
// if the current cluster is the best, then do  nothing
	if(decision==(order+1))
	{
		(*clustIter).updateProfile(1,(*motifIter).second, clustWidth);
		return(0);
	}//end of if
// if a different cluster is desired
	else 
	{
/*		for(j=0;j<4;j++)
		{
			for(k=0;k<clustWidth;k++)
				cout << (*clustIter).profile[j][k] << " ";
			cout <<endl;
		}
*/
		likeminus=loglikelihoodUpdateEven(clustWidth,clustIter,motifIter);
		cout << "likeminus=" <<likeminus;
		if(decision > 0)
		{
// insert it into the desired cluster
			likeadd=clustInsertMotif(clusters, motifIter,decision-1,clustWidth);
			cout << " likeadd1=" <<likeadd<<endl;
		}//end of if
// if a new cluster is desired
		else if (decision == 0)
		{
// add a new cluster to list cluster, with this motif as its only member. 
			tempClustMap[(*motifIter).first]=(*motifIter).second;
			tempClust.cluster=tempClustMap;
			tempClust.iniProfile(clustWidth);
			clusters.push_back(tempClust);
			tempClustMap.clear();
			likeadd=loglikelihoodNewEven(clustWidth,motifIter);//***
			//likeadd=special(clustWidth,clustIter,motifIter);//***
			cout << " likeadd2=" <<likeadd<<endl;
		}//end of else if		
// remove it from its current cluster
		del=(*motifIter).first;
/*		
		a=(*clustIter).cluster.size();
		if((*clustIter).cluster.size()>1)
		{
			del=(*motifIter).first;
			//temMotif=(*motifIter).second;
			//(*clustIter).updateProfile(-1,(*motifIter).second,clustWidth);
			//(*clustIter).cluster.erase(del);
			//(*clustIter).cluster.erase(motifIter);
		}
		else if((*clustIter).cluster.size() == 1)//need to remove the whole cluster
		{
			//(*clustIter).cluster.empty();
			//(*clustIter).cluster.clear();
			clusters.erase(clustIter);
		}
*/
		return(likeadd-likeminus);
	}//end of else
}//end of operation

double loglikelihoodUpdateEven(int width, const list <clustClass>::iterator clustIter,
						   const map <int, motifClass>::iterator motifIter)
{
	short newProfile[4][IN_WIDTH];
	short j,jj,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
	short lsummotif,lb;
	double difference,tem;
	int start;

	start=(*motifIter).second.start;
/*	for (k=0;k<4;k++)
	{
		for (j=0;j<(*motifIter).second.width;j++)
			cout << (*motifIter).second.profile[k][j];
		cout <<endl;
	}
*/
	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
	motifsum=0;
	for (j=0;j<4;j++)
	{
		lsum=lsum+(*clustIter).profile[j][0];		   		
	}
	lsum=2*lsum;
	lsummotif=2*(*motifIter).second.depth;
	lsumnew=lsum+lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
	for (j=0;j<width;j++)
	{
		jj=MIN(j,width-1-j);
		//-if((stage!=1)||((*clustIter).cluster.power[jj]==1))
		//-{
			for (k=0;k<4;k++)
			{
				//a=(*motifIter).second.profile[k][j+start];
				newProfile[k][j]=(*clustIter).profile[k][j]
							 +(*motifIter).second.profile[k][j+start];//***shift to start
			}//end of k
		//-}//end of if
	}//end of j
	for(j=0;j<width/2;j++)
	{
		//-if((stage!=1)||((*clustIter).cluster.power[j]==1))
		//-{
			for(k=0;k<4;k++)
			{
				la=(*clustIter).profile[k][j]+(*clustIter).profile[kstar[k]][width-1-j];
				logsum=logsum+gammaln((double)la+beta[k]);
				lanew=newProfile[k][j]+newProfile[kstar[k]][width-1-j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		//-}//end of if
		//-else if((stage==1)&&((*clustIter).cluster.power[j]==-1))
		//-{
			for(k=0;k<4;k++)
			{
				lb=(*motifIter).second.profile[k][j+(*motifIter).second.start]+
					(*motifIter).second.profile[kstar[k]][width-1-j+(*motifIter).second.start];//***shift to start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		//-}//end of else if
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	return difference;
}//end of loglikelihoodupdateEven

double special(int width, const list <clustClass>::iterator clustIter,
						   const map <int, motifClass>::iterator motifIter)
{
	short newProfile[4][IN_WIDTH];
	short j,jj,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
	short lsummotif,lb;
	double difference,tem;
	int start;

	start=(*motifIter).second.start;
/*	for (k=0;k<4;k++)
	{
		for (j=0;j<(*motifIter).second.width;j++)
			cout << (*motifIter).second.profile[k][j];
		cout <<endl;
	}
*/
	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
	motifsum=0;
	lsum=0;
	lsummotif=2*(*motifIter).second.depth;
	lsumnew=lsummotif;
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
	for (j=0;j<width;j++)
	{
		jj=MIN(j,width-1-j);
		//-if((stage!=1)||((*clustIter).cluster.power[jj]==1))
		//-{
			for (k=0;k<4;k++)
			{
				//a=(*motifIter).second.profile[k][j+start];
				newProfile[k][j]=(*motifIter).second.profile[k][j+start];//***shift to start
			}//end of k
		//-}//end of if
	}//end of j
	for(j=0;j<width/2;j++)
	{
		//-if((stage!=1)||((*clustIter).cluster.power[j]==1))
		//-{
			for(k=0;k<4;k++)
			{
				la=0;
				logsum=logsum+gammaln((double)la+beta[k]);
				lanew=newProfile[k][j]+newProfile[kstar[k]][width-1-j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		//-}//end of if
		//-else if((stage==1)&&((*clustIter).cluster.power[j]==-1))
		//-{
			for(k=0;k<4;k++)
			{
				lb=(*motifIter).second.profile[k][j+(*motifIter).second.start]+
					(*motifIter).second.profile[kstar[k]][width-1-j+(*motifIter).second.start];//***shift to start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		//-}//end of else if
	}//end of j	
	difference=logsumnew-logsum+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	return difference;
}//end of loglikelihoodupdateEven
 
/*
double loglikelihoodUpdateEven(int width, const list <clustClass>::iterator clustIter,
						   const map <int, motifClass>::iterator motifIter)
{
	short newProfile[4][IN_WIDTH];
	short j,jj,k;
	double logsum,la,lsum;
	double logsumnew,lanew, lsumnew;
	double motifsum;
	short lsummotif,lb;
	double difference,tem;

	tem=gammaln(betasum)-sumbeta;
	logsum=0.0;
	logsumnew=0.0;
	lsum=0.0;
	motifsum=0;
	for (j=0;j<4;j++)
	{
		lsum=lsum+(*clustIter).profile[j][0];		   		
	}
	lsum=2*lsum;
	if(lsum==0)
	{
		printf("lsum=0\n");
		exit(0);
	}
	lsummotif=2*(*motifIter).second.depth;
	lsumnew=lsum-lsummotif;
	if(lsumnew<0)
	{
		printf("lsum<lsummotif\n");
		exit(0);
	}
//	printf("%10.5f, %5d, %10.5f\n",lsum,motif[rank].depth,lsumnew);
	for (j=0;j<width;j++)
	{
		jj=MIN(j,width-1-j);
		//-if((stage!=1)||((*clustIter).cluster.power[jj]==1))
		//-{
			for (k=0;k<4;k++)
			{
				newProfile[k][j]=(*clustIter).profile[k][j]
							 -(*motifIter).second.profile[k][j+(*motifIter).second.start];//***shift to start
				if(newProfile[k][j]<0)
				{
					printf("error=%d, %d, cluster size=%d\n",j,k,(*clustIter).cluster.size());
					exit(0);
				}//end of if
			}//end of k
		//-}//end of if
	}//end of j
	for(j=0;j<width/2;j++)
	{
		//-if((stage!=1)||((*clustIter).cluster.power[j]==1))
		//-{
			for(k=0;k<4;k++)
			{
				la=(*clustIter).profile[k][j]+(*clustIter).profile[kstar[k]][width-1-j];
				logsum=logsum+gammaln((double)la+beta[k]);
				lanew=newProfile[k][j]+newProfile[kstar[k]][width-1-j];
				logsumnew=logsumnew+gammaln((double) lanew+beta[k]);
//				printf("%10.5f %10.5f %10.5f %10.5f\n",la,lanew,logsum,logsumnew);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
			logsumnew=logsumnew-gammaln(lsumnew+betasum);
//			printf("**%10.5f %10.5f\n",lsum+betasum,lsumnew+betasum);
		//-}//end of if
		//-else if((stage==1)&&((*clustIter).cluster.power[j]==-1))
		//-{
			for(k=0;k<4;k++)
			{
				lb=(*motifIter).second.profile[k][j+(*motifIter).second.start]+
					(*motifIter).second.profile[kstar[k]][width-1-j+(*motifIter).second.start];//***shift to start
				motifsum=motifsum+gammaln((double) lb+beta[k]);
			}//end of k
			motifsum=motifsum-gammaln((double) lsummotif+betasum)+tem;
		//-}//end of else if
	}//end of j	
	difference=logsum-logsumnew+motifsum;
	//printf("%10.5f %10.5f  ",logsum,logsumnew);
	return difference;
}//end of loglikelihoodupdateEven
*/

double loglikelihoodNewEven(int width, /*const list <clustClass>::iterator clustIter,*/
						   const map <int, motifClass>::iterator motifIter)
{
	short j,k;
	double logsum,la,lsum;
	double difference,tem;
	int a,b;
	int start;

	start=(*motifIter).second.start;
	tem=gammaln(betasum);
	tem=tem-sumbeta;
	lsum=2*(*motifIter).second.depth;
	if(lsum==0)
	{
		printf("lsum=0\n");
		exit(0);
	}
	logsum=0.0;
	for(j=0;j<width/2;j++)
	{
		//-if((stage!=1)||((*clustIter).cluster.power[j]==1))
		//-{
			for(k=0;k<4;k++)
			{
				a=(*motifIter).second.profile[k][j+start];
				b=(*motifIter).second.profile[kstar[k]][width-1-j+start];
				la=(*motifIter).second.profile[k][j+start]+(*motifIter).second.profile[kstar[k]][width-1-j+start];
				logsum=logsum+gammaln((double)la+beta[k]);
			}//end of k
			logsum=logsum-gammaln(lsum+betasum);
		//-}//end of if
	}//end of j	
	difference=2*(logsum+tem*(double)width/2);
	return(difference);
}//end of loglikelihoodNewEven

double clustInsertMotif(list <clustClass> &clusters, const map <int, motifClass>::iterator motifIter,
					  const int decision, const int clustWidth)
{
	int count;
	list <clustClass>::iterator clustIter;
	double likeadd;
			
	count=0;
	//-output(clustWidth,clusters,"clust.txt");
	for(clustIter = (clusters).begin();clustIter != (clusters).end();clustIter ++)
	{
		if(count==decision)
		{
			likeadd=loglikelihoodUpdateEven(clustWidth,clustIter,motifIter);
			(*clustIter).addMotif((*motifIter).second,clustWidth);
			//(*clustIter).updateProfile(1,(*motifIter).second,clustWidth);
			return(likeadd);
		}
		else
			count++;
	}//end of clustIter
}// end of clustInsertMotif

int decide(vector <double> &ratioVec)
{
	int j;
	double sum;
	int number;
	int na,nb,nc;
	double ran;
	int decide;
	double compare, partialsum;

// seed random number generator 
//	srand((unsigned)time(NULL));
//	na=rand() +1;
//	nb=rand() -1;
//	nc=rand() ;
	na=nb=nc=18;

	sum=0.0;
	number=ratioVec.size();
	for(j=0;j<number;j++)
		sum=sum+ratioVec[j];
	//for(j=0;j<number;j++)
	//	ratioVec[j]=ratioVec[j]/sum;
	ran=unran(&na, &nb, &nc);
	compare=ran*sum;
	partialsum=ratioVec[0];
	decide=0;
	while (partialsum < compare)
	{
		decide++;
		partialsum=partialsum+ratioVec[decide];
	}//end of while
	return(decide);
}//end of decide

double bayesRatio(motifClass motif, list <clustClass>::iterator clustIter,int W,int start)
{
	int j,k;
	int Whalf;
	int Lold,Lz,Lnew;
	double sum,sumbayes,tem;
	double newProfile;
	double bayesRatio;
	int a,b;
	int c,d;

	Whalf=W/2;		
	sum=0;
	for(j=0;j<4;j++)
	{
		sum=sum+(*clustIter).profile[j][0];
	}
	Lold=sum*2;
	Lz=motif.depth*2;
	Lnew=Lold+Lz;
// calculate bayes factor
	sumbayes=0.0;
	for(j=0;j<W/2;j++)
	{
		if((*clustIter).power[j]==1)
		{
			for(k=0;k<4;k++)
			{
				a=(*clustIter).profile[k][j];
				b=(*clustIter).profile[kstar[k]][W-1-j];
				c=motif.profile[k][j+motif.start+start];
				d=motif.profile[kstar[k]][W-1-j+motif.start+start];
				newProfile=(*clustIter).profile[k][j]+(*clustIter).profile[kstar[k]][W-1-j]+
							motif.profile[k][j+motif.start+start]+motif.profile[kstar[k]][W-1-j+motif.start+start];
				sumbayes=sumbayes+gammaln((double) newProfile+beta[k])
					-gammaln((double) (*clustIter).profile[k][j]
					+(*clustIter).profile[kstar[k]][W-1-j]+beta[k])
					//-gammaln((double) newProfile+beta[k])
					-gammaln((double) motif.profile[k][j+motif.start+start]
					+motif.profile[kstar[k]][W-1-j+motif.start+start]+beta[k]);
			}//end of k
		}//end of if
	}//end of j
	tem=gammaln((double) Lold+betasum)-gammaln((double) Lnew+betasum)+gammaln((double) Lz+betasum)
		-gammaln(betasum)+sumbeta;
	bayesRatio=sumbayes+tem*(double) Whalf;
	return(bayesRatio);
}//end of bayesRatio

//-----
// calculate likelihood for clusters of even motifs
//-----
double loglikelihoodEven(list <clustClass> clusters, const int width)
{
	int j,k;
	double tem;
	int lsum;
	list <clustClass>::iterator clustIter;
	map <int,motifClass>::iterator motifIter;
	double logsum,la;
	int a,b;
	
	cout << "size=" << clusters.size() <<endl;
	tem=gammaln(betasum)-sumbeta;
	logsum=0;
	for(clustIter = (clusters).begin();clustIter != (clusters).end();clustIter ++)
	{
/*		for(j=0;j<4;j++)
		{
			for(k=0;k<width;k++)
				cout << (*clustIter).profile[j][k];
			cout <<"+++++"<<endl;
		}	
*/
		lsum=0;
		for(j=0;j<4;j++)
			lsum=lsum+(*clustIter).profile[j][0];
		lsum=2*lsum;
		for(j=0;j<width/2;j++)
		{
		//-	if((*clustIter).power[j]==1)
		//-	{
				for(k=0;k<4;k++)
				{
					a=(*clustIter).profile[k][j];
					b=(*clustIter).profile[kstar[k]][width-1-j];
					la = (double) (*clustIter).profile[k][j]+(*clustIter).profile[kstar[k]][width-1-j];
					logsum=logsum+gammaln(la+beta[k]);
				}//end of k
				logsum=logsum-gammaln((double) lsum + betasum);
		//-	}//end of if
		/*-	else if ((*clustIter).power[j]==0)
			{
				for(motifIter = (*clustIter).cluster.begin();motifIter != (*clustIter).cluster.end();motifIter ++)
				{
					sumin=2*(*motifIter).second.depth;
					for(k=0;k<4;k++)
					{
						la = (double) (*motifIter).second.profile[k][j]+(*motifIter).second.profile[kstar[k]][width-1-j]+beta[k];
						logsum=logsum+gammaln(la+beta[k]);
					}//end of k
					logsum=logsum-gammaln((double) sumin + betasum);
				}//end of motifIter
			}//end of else if
			*/
		}//end of j
	}//end of clustIter
	logsum=logsum+clusters.size()*(double) (width/2) *tem;
	return(logsum);
}//end of loglikelihoodEven

//----- perform clustering update on the current motif, remove it from the current cluster, -----//
//		and put it into the cluster which fits it the best									-----//
//_____																						-----//	

//----- read in motif block data and put into map and list -----//
int readSeq(int clustWidth, list <clustClass> &allCluster)
{
	map <int, motifClass> tempClustMap;
	map <int, motifClass>::iterator motifIter;
	//list < clustClass > allCluster;
	list < clustClass >::iterator allClustIter;
	vector <string> seq;
	//vector <motifClass> motifVector;
	motifClass temMotif;
	clustClass tempClust;
	const string inputFilename = "fileIn1.txt";//dat
	const string outputFilename = "out.txt";
	ifstream inFile;
	string lineString;
	istringstream iss;
	string motif;
	char c;
	int rowCount;
	int motifCount;
	int count;
	string tempName;
	ofstream outputFile;
	int signal=0;//signal for if it is the first line of the sfirst motif

	cout << "start reading in data" << endl;
	inFile.open(inputFilename.c_str());
	if(!inFile)
	{
		cout << "ERROR: Unable to open file: " << inputFilename << endl;
	    exit(2);
	}//end of if
// start reading in data
	motifCount=0;
	//temMotif.name=" ";
	signal=0;
	while (inFile)
	{
//-		output(clustWidth,allCluster,"clust.txt");
// get name
		getline(inFile,lineString);
// get DNA sequence
		iss.clear();
		iss.str(lineString + " ");
		iss >> c;
		if(!iss)
		{
			//cout << "blank line" <<endl;
			temMotif.width=motif.length();
			temMotif.depth=rowCount;
			temMotif.start=(temMotif.width-clustWidth)/2;
			if(temMotif.width %2 ==0)
				temMotif.isEven=0;
			else
				temMotif.isEven=1;
			temMotif.basepair=seq;
			seq.clear();
			temMotif.genProfile();
			//motifVector.push_back(temMotif);
			if(motifCount<INITIAL_CLUSTER)
			{
				tempClustMap[motifCount]=temMotif;
				tempClust.cluster=tempClustMap;
				tempClust.iniProfile(clustWidth);
				allCluster.push_back(tempClust);
				tempClustMap.clear();
			}
			else
			{
				if (motifCount % INITIAL_CLUSTER == 0)
					allClustIter = allCluster.begin();
				(*allClustIter).addMotif(temMotif,clustWidth);
				allClustIter ++;
			}
			motifCount++;
		}
		else if(c=='>')
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
		else
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
	outputFile.open(outputFilename.c_str());
	if(!outputFile)
	{
		cout << "ERROR: Unable to open file: " << outputFilename << endl;
	    exit(2);
	}//end of if
	count=0;
	for(allClustIter = allCluster.begin();allClustIter != allCluster.end();allClustIter ++)
	{
		cout << "cluster " << count << endl; 
		count++;		
		(*allClustIter).printClust(clustWidth, outputFile);
	}
	outputFile.close();
	return(motifCount);
}//end of readSeq

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