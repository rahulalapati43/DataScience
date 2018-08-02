/*************************************************************************
 *                           clustQ
 * clustQ is a highly efficient superposition-free approach to protein 
 * decoy clustering based on weighted internal distance comparisons.
 *
 * Developed by: Rahul Alapati
 * 
 * This program is for multi-model pairwise comparisons for rapid clustering 
 * of large number of protein decoys based on novel decoy-native similarity 
 * metric called WQ-score. Users need to provide the amino acid sequence of 
 * the target protein and an archive consisting of multiple structural decoys. 
 * All decoy structures must be in the PDB format. Users may obtain a brief 
 * instruction by running the program in help mode, with -h argument.
 * For comments, please email to bhattacharyad@auburn.edu.
 * 
 * Copyright (C) 2018 Debswapna Bhattacharya. All rights reserved.
 * Permission to use, copy, modify, and distribute this program for 
 * any purpose, with or without fee, is hereby granted, provided that
 * the notices on the head, the reference information, and this
 * copyright notice appear in all copies or substantial portions of 
 * the Software. It is provided "as is" without express or implied 
 * warranty.
 ************************** Change log ***********************************
 *     06/01/2018: The first version released.
 *************************************************************************/
 
#include <dirent.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <cfloat>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <bits/stdc++.h>
#include <time.h>

using namespace std;
using std::vector;

std::string line;
std::string prefix("ATOM");
std::string symbol(">");

const double PI = std::atan(1.0)*4;
const double RAD2DEG = 180/PI;
const double DEG2RAD = PI/180;

//global variables

vector<string> files = vector<string>();

const int NR_SEP = 1;
const int SR_SEP = 6;
const int MR_SEP = 12;
const int LR_SEP = 24;

double wt_nr = 1.0;
double wt_sr = 2.0;
double wt_mr = 4.0;
double wt_lr = 8.0;

char modFile[200] = "";
char natFile[200] = "";
char decoyFile[200] = "";
char fastaFile[200] = "";
char targetFile[200] = "";
char outputFile[200] = "";

// amino acid residue types
const int ALA = 0;
const int CYS = 1;
const int ASP = 2;
const int GLU = 3;
const int PHE = 4;
const int GLY = 5;
const int HIS = 6;
const int ILE = 7;
const int LYS = 8;
const int LEU = 9;
const int MET = 10;
const int ASN = 11;
const int PRO = 12;
const int GLN = 13;
const int ARG = 14;
const int SER = 15;
const int THR = 16;
const int VAL = 17;
const int TRP = 18;
const int TYR = 19;

// amino acide residue code to three letter format
string seq3[20] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

// amino acide residue code to one letter format
string seq[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

// point3d object
   struct point3d {
          double x;
          double y;
          double z;
   };

// pdbInfo object
   struct pdbInfo {
          int id;
          int aa;
          point3d ca;
   };

// poseInfo object
   struct poseInfo {
          int id;
          int aa;
          double r;
          double theta;
          double phi;
   };

//indicators
bool mFile = false;
bool nFile = false;
bool tFile = false;
bool oFile = false;
                                             
//function prototypes
void parseNextItem(int argc, char ** argv, int & i);
void parseCommandLine(int argc, char ** argv);
float calculateQScore(char *modFile, char *natFile);
int getAA(const char * aa);
void loadPdb(char *filename, vector<pdbInfo> &pdb);
double getDistance(point3d & p1, point3d &p2);
string getSequence(char *filename);
string getfastaSequence(char *filename);

/*************************************************************************
 * Name        : main
 * Purpose     : Compute CLUSTQ
 * Arguments   : int argc, char ** argv
 * Return Type : int
 *************************************************************************/
int main(int argc, char ** argv) {
					
  //display header info
   cout << endl;
   cout << "**********************************************************************" << endl;
   cout << "*                            clustQ                                  *" << endl;
   cout << "* Efficient protein decoy clustering using superposition-free        *" << endl;
   cout << "*            weighted internal distance comparisons.                 *" << endl;
   cout << "* For comments, please email to bhattacharyad@auburn.edu             *" << endl;
   cout << "**********************************************************************" << endl;
   
   vector< pair <float,string> > vect;
   time_t rawtime, rawtime2;
   struct tm * timeinfo;
   time(&rawtime);
   timeinfo = localtime(&rawtime);
 
   parseCommandLine(argc, argv);
   
   DIR* dirp = opendir(decoyFile);
   struct dirent * dp;

   if (dirp != NULL) {
       while ((dp = readdir(dirp)) != NULL){	
              if (strcmp(dp->d_name, ".") != 0 && strcmp(dp->d_name, "..") != 0 ) {
                  files.push_back(string(dp->d_name));
              }
       }
       closedir(dirp);
   }

   float QSCORE = 0.0;
   string decoys[files.size()];
   float averages[files.size()];
   string fastaseq;
   string pdbseq;
   string pdb2seq;
   
   printf("Processing started at: %s", asctime(timeinfo));
   cout << "Reading and Scoring Decoys..." << endl;
   fastaseq = getfastaSequence(fastaFile);
   for (unsigned int i = 0; i < files.size(); i++){
        double QScoreAverage = 0.0;
        double QScoreSum = 0.0;
        int count = 0;
        char model[200] = "";
	strcpy(modFile, files[i].c_str());
        strcat(model, decoyFile);
	strcat(model, modFile);
        decoys[i] = modFile;
        pdbseq = getSequence(model);
        if (pdbseq == fastaseq)
        {
	    for (unsigned int j = 0; j < files.size(); j++){
	         char native[200] = "";		
	         strcpy(natFile, files[j].c_str());
                 strcat(native, decoyFile);
                 strcat(native, natFile);
	         if (strcmp(model,native) != 0){
                     pdb2seq = getSequence(native);
                     if (pdb2seq == fastaseq)
                     { 
		     	QSCORE = calculateQScore(model,native);
                     	QScoreSum = QScoreSum + QSCORE;
		     	count = count + 1;
                     }
	         }	    
	    }
            QScoreAverage = QScoreSum/count;
            averages[i] = QScoreAverage;
        }
	else
	{
		averages[i] = 0.0;
	}  
   }
   
   int n = sizeof(averages)/sizeof(averages[0]);
   
   //Entering values in vector of pairs
   for (int r = 0; r < n; r++) 
	vect.push_back( make_pair(averages[r],decoys[r]));
   
   sort(vect.begin(), vect.end());

   cout << "done" << endl;
   printf("Writing output to %s...\n",outputFile);
   ofstream outFile;
   outFile.open(outputFile);
   
   for (int s = n-1; s >= 0; s--) {
        if (vect[s].first == 0.0)
		outFile << vect[s].second << " " << "X" << endl;
        else
        	outFile << vect[s].second << " " << fixed << vect[s].first << endl;
   }

   printf("done \n") ;
   time(&rawtime2);
   timeinfo = localtime(&rawtime2);
   printf("Processing finished at: %s",asctime(timeinfo));
   cout << "Total processing time: " << rawtime2-rawtime << " second(s)" << endl;
}
 
/*************************************************************************
 * Name        : parseCommandLine
 * Purpose     : parse command line arguments
 * Arguments   : int argc, char ** argv
 * Return Type : void
 *************************************************************************/  
void parseCommandLine(int argc, char ** argv) {
   int i = 1;
    while (i < argc)
        parseNextItem(argc, argv, i);
    if (!mFile) {
        cout << endl;
        cout << "Error! Decoy file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                    cout << "   -t target : target name" << endl;
                    cout << "   -f fasta  : fasta file" << endl;
		    cout << "   -d decoy  : decoy directory" << endl;
		    cout << "   -o output : output file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
    else {
        ifstream fin(decoyFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Decoy file not present" << endl << endl;
                                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                                    cout << "   -t target : target name" << endl;
                                    cout << "   -f fasta  : fasta file" << endl;
				    cout << " 	-d decoy  : decoy directory" << endl;
				    cout << " 	-o output : output file" << endl;
                                    cout << "   -h help   : this message" << endl;
                                    exit(0);
        }
    }
    if (!nFile) {
        cout << endl;
        cout << "Error! Fasta file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                    cout << "   -t target : target name" << endl;
                    cout << "   -f fasta  : fasta file" << endl;
		    cout << "   -d decoy  : decoy directory" << endl;
		    cout << "   -o output : output file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
    else {
        ifstream fin( fastaFile);
        if( fin.fail() ) {
            cout << endl;
            cout << "Error! Fasta file not present" << endl << endl;
                                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                                    cout << "   -t target : target name" << endl;
                                    cout << "   -f fasta  : fasta file" << endl;
				    cout << "   -d decoy  : decoy directory" << endl;
				    cout << "   -o output : output file" << endl;
                                    cout << "   -h help   : this message" << endl;
                                    exit(0);
        }
    }
    if (!tFile) {
        cout << endl;
        cout << "Error! Target file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                    cout << "   -t target : target name" << endl;
                    cout << "   -f fasta  : fasta file" << endl;
                    cout << "   -d decoy  : decoy directory" << endl;
                    cout << "   -o output : output file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
    if (!oFile) {
        cout << endl;
        cout << "Error! Output file must be provided" << endl << endl;
                    cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
                    cout << "   -t target : target name" << endl;
                    cout << "   -f fasta  : fasta file" << endl;
                    cout << "   -d decoy  : decoy directory" << endl;
                    cout << "   -o output : output file" << endl;
                    cout << "   -h help   : this message" << endl;
                    exit(0);
    }
}

/*************************************************************************
 * Name        : parseNextItem
 * Purpose     : parse next item in command line argument
 * Arguments   : int argc, char ** argv, int & i
 * Return Type : void
 **************************************************************************/
 void parseNextItem(int argc, char ** argv, int & i) {
        if (strncmp(argv[i], "-d", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No decoy file provided" << endl << endl;
            cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
            cout << "   -d decoy  : decoy directory" << endl;
  	    cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        mFile = true;
        strcpy(decoyFile, argv[++i]);
        strcat(decoyFile, "/");
    }
    else if (strncmp(argv[i], "-f", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No fasta file provided" << endl;
            cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
	    cout << "   -d decoy  : decoy directory" << endl;
 	    cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        nFile = true;
        strcpy(fastaFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-h", 2) == 0) {
            cout << endl;
            cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
 	    cout << "   -d decoy  : decoy directory" << endl;
	    cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
    }
    else if (strncmp(argv[i], "-t", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
	    cout << "Error! No target file provided" << endl;
            cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
            cout << "   -d decoy  : decoy directory" << endl;
            cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
        tFile = true;
	strcpy(targetFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-o", 2) == 0) {
	if (argc < i + 2) {
            cout << endl;
	    cout << "Error! No output file provided" << endl;
            cout << "Usage: " << argv[0] << " -t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
            cout << "   -d decoy  : decoy directory" << endl;
            cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);
        }
	oFile = true;
	strcpy(outputFile, argv[++i]);
    }
     else {
            cout << endl;
            cout << "Error! Invalid option" << endl << endl;
            cout << "Usage: " << argv[0] << "-t target -f fasta -d decoy -o output" << endl;
            cout << "   -t target : target name" << endl;
            cout << "   -f fasta  : fasta file" << endl;
	    cout << "   -d decoy  : decoy directory" << endl;
	    cout << "   -o output : output file" << endl;
            cout << "   -h help   : this message" << endl;
            exit(0);

    }
    i++;
} 

/*************************************************************************
 * Name        : calculateQScore
 * Purpose     : calculate QSCORE
 * Arguments   : char *modFile, char *natFile
 * Return Type : float
 **************************************************************************/
float calculateQScore(char *modFile, char *natFile) { 
      // load model
         vector<pdbInfo> modPdb;
         loadPdb(modFile, modPdb);
      
      //load native
        vector<pdbInfo> natPdb;
        loadPdb(natFile, natPdb);

      	if (modPdb.size() != natPdb.size()) {
	    cout << modFile <<endl;
            cout <<natFile << endl;
      	    cout << "Error! size mismatch between model and native" << endl;
            return 0.0;
      	}

      	double num_nr = 0.0;
      	double sum_nr = 0.0;
      	for (int i = 0; i < modPdb.size(); i++) {
             for (int j = i + NR_SEP; j < modPdb.size(); j++) {
                  if (abs(i - j) < SR_SEP) {
                        sum_nr += exp ( - pow((getDistance(modPdb[i].ca, modPdb[j].ca) - getDistance(natPdb[i].ca, natPdb[j].ca)), 2.0));
                        num_nr += 1.0;
                }
            }
          }

      	double QNR = 0.0;
    	if (num_nr > 0.0){
            QNR = sum_nr / num_nr;
    	}
    	else {
        	wt_nr = 0.0;
    	}

    	double num_sr = 0.0;
    	double sum_sr = 0.0;

    	for (int i = 0; i < modPdb.size(); i++) {
             for (int j = i + SR_SEP; j < modPdb.size(); j++) {
                  if (abs(i - j) < MR_SEP) {
                        sum_sr += exp ( - pow((getDistance(modPdb[i].ca, modPdb[j].ca) - getDistance(natPdb[i].ca, natPdb[j].ca)), 2.0));
                        num_sr += 1.0;
                }
            }
          }
	
	double QSR = 0.0;
    	if (num_sr > 0.0){
            QSR = sum_sr / num_sr;
    	}	
    	else{
        	wt_sr = 0.0;
    	}

    	double num_mr = 0.0;
    	double sum_mr = 0.0;

    	for (int i = 0; i < modPdb.size(); i++) {
             for (int j = i + MR_SEP; j < modPdb.size(); j++) {
                  if (abs(i - j) < LR_SEP) {
                        sum_mr += exp ( - pow((getDistance(modPdb[i].ca, modPdb[j].ca) - getDistance(natPdb[i].ca, natPdb[j].ca)), 2.0));
                        num_mr += 1.0;
                }
            }
          }

    	double QMR = 0.0;
    	if (num_mr > 0.0){
            QMR = sum_mr / num_mr;
    	}
    	else{
        	wt_mr = 0.0;
    	}

    	double num_lr = 0.0;
    	double sum_lr = 0.0;
    	for (int i = 0; i < modPdb.size(); i++) {
             for (int j = i + LR_SEP; j < modPdb.size(); j++) {
                	sum_lr += exp ( - pow((getDistance(modPdb[i].ca, modPdb[j].ca) - getDistance(natPdb[i].ca, natPdb[j].ca)), 2.0));
                	num_lr += 1.0;
            }
          }

        double QLR = 0.0;
    	if (num_lr > 0.0){
           QLR = sum_lr / num_lr;
    	}
    	else {
        	wt_lr = 0.0;
    	}

        double QSCORE = ((wt_nr * QNR) + (wt_sr * QSR) + (wt_mr * QMR) + + (wt_lr * QLR)) / (wt_nr + wt_sr + wt_mr + wt_lr);

 return QSCORE;
}

/*************************************************************************
 * Name        : getAA
 * Purpose     : convert AA name to a numerical code
 * Arguments   : const char * aa
 * Return Type : int
 **************************************************************************/

int getAA(const char * aa) {
    if (strlen(aa) == 3) {
        if (strcmp(aa, "ALA") == 0)
            return (ALA);
        else if (strcmp(aa, "ARG") == 0)
            return (ARG);
        else if (strcmp(aa, "ASN") == 0)
            return (ASN);
        else if (strcmp(aa, "ASP") == 0)
            return (ASP);
        else if (strcmp(aa, "CYS") == 0)
            return (CYS);
        else if (strcmp(aa, "GLN") == 0)
            return (GLN);
        else if (strcmp(aa, "GLU") == 0)
            return (GLU);
        else if (strcmp(aa, "GLY") == 0)
            return (GLY);
        else if (strcmp(aa, "HIS") == 0)
            return (HIS);
        else if (strcmp(aa, "ILE") == 0)
            return (ILE);
        else if (strcmp(aa, "LEU") == 0)
            return (LEU);
        else if (strcmp(aa, "LYS") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
        else if (strcmp(aa, "PHE") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
        else if (strcmp(aa, "PHE") == 0)
            return (LEU);
        else if (strcmp(aa, "LYS") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
        else if (strcmp(aa, "PHE") == 0)
            return (PHE);
        else if (strcmp(aa, "PRO") == 0)
     	    return (PRO);
        else if (strcmp(aa, "SER") == 0)
            return (SER);
        else if (strcmp(aa, "THR") == 0)
            return (THR);
        else if (strcmp(aa, "TRP") == 0)
            return (TRP);
        else if (strcmp(aa, "TYR") == 0)
            return (TYR);
        else if (strcmp(aa, "VAL") == 0)
            return (VAL);
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else if (strlen(aa) == 1) {
        if (aa[0] == 'A')
            return ALA;
        else if (aa[0] == 'C')
            return CYS;
        else if (aa[0] == 'D')
            return ASP;
        else if (aa[0] == 'E')
            return GLU;
        else if (aa[0] == 'F')
            return PHE;
        else if (aa[0] == 'G')
            return GLY;
        else if (aa[0] == 'H')
            return HIS;
        else if (aa[0] == 'I')
            return ILE;
        else if (aa[0] == 'K')
            return LYS;
        else if (aa[0] == 'L')
            return LEU;
        else if (aa[0] == 'M')
            return MET;
        else if (aa[0] == 'N')
            return ASN;
        else if (aa[0] == 'P')
   	    return PRO;
        else if (aa[0] == 'Q')
            return GLN;
        else if (aa[0] == 'R')
            return ARG;
        else if (aa[0] == 'S')
            return SER;
        else if (aa[0] == 'T')
            return THR;
        else if (aa[0] == 'V')
            return VAL;
        else if (aa[0] == 'W')
            return TRP;
        else if (aa[0] == 'Y')
            return TYR;
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid amino acid " << aa << endl;
        exit(0);
    }
}


/**************************************************************************
 * Name        : loadPdb
 * Purpose     : loads a pdb file into a vector of pdbInfo object
 * Arguments   : char *filename, vector<pdbInfo> &pdb
 * Return Type : void
 **************************************************************************/

void loadPdb(char *filename, vector<pdbInfo> &pdb) {
    string line, str;
    string atom ("ATOM ");
    int prevRes = -999999;
    point3d caAtom;
    pdbInfo pdbData;
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            if(line.compare(0, atom.length(), atom)==0) {
                int res = atoi(line.substr(22, 4).c_str());
                int anmae = getAA(line.substr(17, 3).c_str());
		// seek for the next residue
		if (res != prevRes) {
                    if (prevRes != -999999) {
			// insert the data collected so far
			 pdbData.ca = caAtom;
                        pdb.push_back(pdbData);
                    }
                    prevRes = res;
                    pdbData.id = res;
                    pdbData.aa = anmae;
                }
		// consider the first alternate location id (i.e. A) if present
		if ( line.compare(16, 1, " ") == 0 || line.compare(16, 1, "A") == 0 ) {
		   // get the CA atom coordinate
		   if( line.compare(12, 4, "CA  ") == 0 || line.compare(12, 4, " CA ") == 0 || line.compare(12, 4, "  CA") == 0 ) {
                        caAtom.x = atof(line.substr(30, 8).c_str());
                        caAtom.y = atof(line.substr(38, 8).c_str());
                        caAtom.z = atof(line.substr(46, 8).c_str());
                    }
                }
            }
        }
        fin.close();
	// for the last residue
	 pdbData.ca = caAtom;
        pdb.push_back(pdbData);
    }
    else {
        cout << "Error! pdb file can not open " << filename << endl;
        exit(0); 
    }

}
 
/*************************************************************************
 * Name        : getDistance
 * Purpose     : gets distance between two points
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 **************************************************************************/
double getDistance(point3d & p1, point3d &p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/*************************************************************************
 * Name        : getSequence
 * Purpose     : gets sequence from pdb file
 * Arguments   : char *filename
 * Return Type : string
 **************************************************************************/

string getSequence(char *filename) {
    ifstream inFile;
    inFile.open(filename);
    string aminoacid;
    string res;
    string sequence;
    string prev = "-1";
    string pos;

    if (inFile) {
        for(std::string line; getline(inFile,line); )
        {
          if (std::equal(prefix.begin(), prefix.end(), line.begin())) {
              aminoacid = line.substr(17,3);
              pos = line.substr(22,4);
              if (aminoacid.length() == 3) {
                if (aminoacid == "ALA")
                        res = "A";
                else if (aminoacid == "CYS")
                        res = "C";
                else if (aminoacid == "ASP")
                        res = "D";
                else if (aminoacid == "GLU")
                        res = "E";
                else if (aminoacid == "PHE")
                        res = "F";
                else if (aminoacid == "GLY")
                        res = "G";
                else if (aminoacid == "HIS")
                        res = "H";
                else if (aminoacid == "ILE")
                        res = "I";
                else if (aminoacid == "LYS")
                        res = "K";
                else if (aminoacid == "LEU")
                        res = "L";
                else if (aminoacid == "MET")
                        res = "M";
                else if (aminoacid == "ASN")
                        res = "N";
                else if (aminoacid == "PRO")
                        res = "P";
                else if (aminoacid == "GLN")
                        res = "Q";
                else if (aminoacid == "ARG")
                        res = "R";
                else if (aminoacid == "SER")
                        res = "S";
                else if (aminoacid == "THR")
		         res = "T";
                else if (aminoacid == "VAL")
                        res = "V";
                else if (aminoacid == "TRP")
                        res = "W";
                else if (aminoacid == "TYR")
                        res = "Y";
                else
                        res = "X";
            }
            if (pos != prev) {
                sequence += res;
                prev = pos;
            }
          }
        }
    }
  return sequence;
}

/*************************************************************************
 * Name        : getfastaSequence
 * Purpose     : gets fasta sequence from fasta file
 * Arguments   : char *filename
 * Return Type : string
 ***************************************************************************/

string getfastaSequence(char *filename) {
    ifstream fastaFile;
    fastaFile.open(filename);
    string fastaseq;

    if (fastaFile) {
        for(std::string line; getline(fastaFile,line); )
        {
          if (!(std::equal(symbol.begin(), symbol.end(), line.begin()))) { 
              fastaseq += line;
          }
        }
    }
    else {
	printf("Unable to open: %s.\n",filename);
	exit(0);
    }
 return fastaseq;
}	
