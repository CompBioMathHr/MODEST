#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <algorithm>

const int AL=21; // alfabet length
const int MQL=80; // maximal query length
const int MNM=10000; // maximal number of motifs

int mw; // motif width
int nos; // number of sequences
int noi; // number of iterations
int nkeys=0; // number of keys
float th; // threshold

float scale; 

char alfabet[21]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','x'}; 
float background[20]={0.078, 0.051, 0.043, 0.053, 0.019, 0.043, 0.063, 0.072, 0.023, 0.053, 0.091, 0.059, 0.022, 0.039, 0.052, 0.068, 0.059, 0.014, 0.032, 0.066}; 

float PAM120[20][20]={ {0.7128871,   0.0029970,   0.0101898,   0.0138861,   0.0027972,   0.0071928,   0.0207792,   0.0453546,   0.0027972,   0.0055944,   0.0076923,   0.0063936,   0.0023976,   0.0025974,   0.0276723,   0.0548452,   0.0443556,   0.0000999,   0.0024975,   0.0269730},
                       {0.0061944,   0.8008792,   0.0043960,   0.0009991,   0.0024978,   0.0193826,   0.0015986,   0.0033969,   0.0175842,   0.0045959,   0.0031971,   0.0795284,   0.0025977,   0.0024978,   0.0116895,   0.0227795,   0.0061944,   0.0045959,   0.0002997,   0.0050954},
                       {0.0223754,   0.0042953,   0.6336030,   0.0757167,   0.0005993,   0.0100889,   0.0192788,   0.0271701,   0.0345620,   0.0062931,   0.0069923,   0.0499451,   0.0005993,   0.0026970,   0.0060933,   0.0619319,   0.0273699,   0.0001998,   0.0063930,   0.0037958},  
                       {0.0235764,   0.0008991,   0.0649351,   0.7016983,   0.0001998,   0.0138861,   0.1050949,   0.0254745,   0.0080919,   0.0027972,   0.0009990,   0.0157842,   0.0002997,   0.0002997,   0.0035964,   0.0180819,   0.0103896,   0.0000999,   0.0005994,   0.0031968},  
                       {0.0078921,   0.0025974,   0.0007992,   0.0002997,   0.9314685,   0.0002997,   0.0003996,   0.0032967,   0.0023976,   0.0046953,   0.0003996,   0.0004995,   0.0001998,   0.0003996,   0.0028971,   0.0233766,   0.0034965,   0.0000999,   0.0071928,   0.0072927},  
                       {0.0185814,   0.0214785,   0.0102897,   0.0167832,   0.0002997,   0.7265734,   0.0675325,   0.0081918,   0.0407592,   0.0027972,   0.0135864,   0.0269730,   0.0041958,   0.0003996,   0.0175824,   0.0103896,   0.0075924,   0.0000999,   0.0004995,   0.0053946},  
                       {0.0352612,   0.0013985,   0.0162821,   0.0994906,   0.0001998,   0.0521426,   0.7117171,   0.0175807,   0.0042953,   0.0044951,   0.0030966,   0.0165818,   0.0003996,   0.0002997,   0.0076915,   0.0146838,   0.0065927,   0.0000999,   0.0022975,   0.0053941},  
                       {0.0452502,   0.0005993,   0.0133853,   0.0139846,   0.0003996,   0.0028968,   0.0100889,   0.8458695,   0.0005993,   0.0006992,   0.0027969,   0.0055938,   0.0002997,   0.0024973,   0.0058935,   0.0341624,   0.0070922,   0.0000999,   0.0002997,   0.0074918},  
                       {0.0061944,   0.0218803,   0.0403637,   0.0108902,   0.0024978,   0.0468578,   0.0072934,   0.0036967,   0.7985813,   0.0007993,   0.0096913,   0.0076931,   0.0003996,   0.0048956,   0.0116895,   0.0068938,   0.0035968,   0.0000999,   0.0090918,   0.0068938},  
                       {0.0156827,   0.0066926,   0.0063930,   0.0030966,   0.0025971,   0.0028968,   0.0066926,   0.0017980,   0.0006992,   0.7186095,   0.0490460,   0.0095895,   0.0106882,   0.0171811,   0.0030966,   0.0061932,   0.0230746,   0.0000999,   0.0028968,   0.1129757},  
                       {0.0097883,   0.0026968,   0.0024970,   0.0005993,   0.0001998,   0.0066920,   0.0027966,   0.0029964,   0.0025969,   0.0201758,   0.8723532,   0.0032960,   0.0169796,   0.0140831,   0.0049940,   0.0031962,   0.0053935,   0.0000999,   0.0027966,   0.0257691},  
                       {0.0066933,   0.0408591,   0.0258741,   0.0082917,   0.0002997,   0.0133866,   0.0095904,   0.0057942,   0.0035964,   0.0047952,   0.0055944,   0.8230769,   0.0083916,   0.0003996,   0.0052947,   0.0165834,   0.0176823,   0.0001998,   0.0002997,   0.0032967},  
                       {0.0140873,   0.0094915,   0.0014987,   0.0007993,   0.0002997,   0.0086922,   0.0030972,   0.0032970,   0.0006994,   0.0255770,   0.0956139,   0.0416625,   0.7219502,   0.0093915,   0.0031971,   0.0092916,   0.0135878,   0.0000999,   0.0005995,   0.0370666},  
                       {0.0051948,   0.0024975,   0.0026973,   0.0002997,   0.0002997,   0.0003996,   0.0003996,   0.0026973,   0.0047952,   0.0150849,   0.0303696,   0.0005994,   0.0024975,   0.8655345,   0.0025974,   0.0067932,   0.0029970,   0.0024975,   0.0477522,   0.0039960},  
                       {0.0468531,   0.0092907,   0.0055944,   0.0035964,   0.0026973,   0.0131868,   0.0076923,   0.0091908,   0.0070929,   0.0007992,   0.0072927,   0.0080919,   0.0003996,   0.0002997,   0.8210789,   0.0357642,   0.0130869,   0.0000999,   0.0002997,   0.0075924},  
                       {0.0685246,   0.0124863,   0.0364599,   0.0126860,   0.0106882,   0.0052942,   0.0099890,   0.0446509,   0.0034962,   0.0031965,   0.0033963,   0.0189791,   0.0022975,   0.0044951,   0.0253721,   0.6655679,   0.0611328,   0.0021976,   0.0025971,   0.0064929},  
                       {0.0646354,   0.0034965,   0.0189810,   0.0080919,   0.0028971,   0.0049950,   0.0060939,   0.0106893,   0.0028971,   0.0146853,   0.0079920,   0.0242757,   0.0043956,   0.0027972,   0.0106893,   0.0725275,   0.7149850,   0.0001998,   0.0024975,   0.0221778},  
                       {0.0006993,   0.0182817,   0.0024975,   0.0002997,   0.0001998,   0.0003996,   0.0001998,   0.0003996,   0.0025974,   0.0002997,   0.0096903,   0.0010989,   0.0001998,   0.0073926,   0.0003996,   0.0107892,   0.0005994,   0.9386613,   0.0049950,   0.0002997},  
                       {0.0050949,   0.0003996,   0.0083916,   0.0007992,   0.0071928,   0.0004995,   0.0023976,   0.0005994,   0.0091908,   0.0029970,   0.0059940,   0.0026973,   0.0002997,   0.0636364,   0.0003996,   0.0051948,   0.0048951,   0.0024975,   0.8719281,   0.0048951},  
                       {0.0374625,   0.0026973,   0.0030969,   0.0029970,   0.0046953,   0.0027972,   0.0051948,   0.0120879,   0.0023976,   0.0654346,   0.0350649,   0.0034965,   0.0087912,   0.0011988,   0.0053946,   0.0067932,   0.0200799,   0.0000999,   0.0024975,   0.7777223}
                     }; 

struct atom {
	char sequence[80]; 
	int row; 
	int rowpos; // row position
	float Evalue; 
	long double postScore; 
	}; 
	
float f(int x) {
	float y; 
	if (x<36) y=0.505514-0.00551429*((float)x); 
	else y=-9.05/((float)x*x) +9.5/((float)x)+0.05; 
	return y; 
	}

int findindex(char letter) { // find index of a letter in an alfabet
	for(int i=0; i<21; i++) if (letter==alfabet[i]) return i; 
	return -1; 
	}

int main(int argc, char *argv[]) {
	
	struct atom	motif[MNM];
	
	if (argc==6) {
				
		FILE *file; 
		
		// open query 
		int ok=1; 
		int keys[MQL]; 
		file=fopen(argv[2],"r"); 
		nos=0; // set the number of sequences to zero
 		if (file==0) printf("Can't open file %s!\n", argv[2]); 
		else {
			char tempquery[MQL]; // MQL = maximal query length
			int i=0;   
			int isint=0;  
			int ischar=0; 
			int nuintrows=0; 
			char c; 
			while ( (c=getc(file))!=EOF ) { 
				if (c=='X') c='x'; 
				if (findindex(c)!=-1) { tempquery[i]=c; i++; ischar=1; }
				else { 
					if (c=='0' || c=='1') {
						isint=1;
						if (c=='1') {keys[nkeys]=i; nkeys++;}
						i++; 
						} 
					else {
						if (isint==1) { nuintrows++; isint=0; }  
						else {
							tempquery[i]='\0'; 
							strcpy(motif[nos].sequence,tempquery); 
							i=0; 
							nos++; 
							ischar=0;
							isint=0;
							}  
					    }
					 }
				if (isint==1 && ischar==1) { printf("Query is not properly formated!\n"); ok=0; }  
				if (nuintrows>1) { printf("Query is not properly formated!\n"); ok=0; } 
				}
			mw=strlen(motif[0].sequence); // motif width
			for(int j=0; j<nos; j++) if (mw!=strlen(motif[j].sequence)) ok=0; 
			}
		fclose(file); 
		
		// how many (real) rows does fasta or text proteome have
		file=fopen(argv[1],"r"); 
		
		int countnames=0; 
		int countinnerrows=0;
		int nurows; // number of rows
		 
		if (file==0) { printf("Can't open file %s!\n", argv[1]); ok=0; } 
		else {
			char c; 
			char cold='\n'; 
			while ( (c=getc(file))!=EOF ) { 
				if (c=='>' && cold=='\n') {
					countnames++;
					countinnerrows=0; 
					} 
				if (c=='\n') {
					countinnerrows++; 
					}
				cold=c; 
				} 
			} 
		fclose(file); 
		
		if (countnames>0) nurows=countnames; 
		else nurows=countinnerrows;

		
		int *rowlength; // row length array
		rowlength=new int[nurows]; 
		int *rownamelength; // rowname length array
		rownamelength=new int[nurows]; 
		
		for(int i=0; i<nurows; i++) {
			rowlength[i]=0; 
			rownamelength[i]=0; 
			}
		
		file=fopen(argv[1],"r");
		if (file==0) printf("Can't open file %s!\n", argv[1]); 
		else {
			if (countnames>0) { // fasta format
				int i=-1; 
				int state=0;
				char c;  
				char cold='\n'; 
				while ( (c=getc(file))!=EOF ) {
					if (c=='>' && cold=='\n') state=1;   
					if (state==1 && c!='\n') rownamelength[i+1]++;
					if (state==1 && c=='\n') {state=0; i++;}
					if (state==0 && findindex(c)!=-1 && c!='x') rowlength[i]++;					
					cold=c; 
					}
				} else { // plain text format
					int i=0; 
					char c; 
					while ( (c=getc(file))!=EOF ) {
						if (c=='\n') i++; 
						if (findindex(c)!=-1 && c!='x') rowlength[i]++; 
						}
					}
			}	
		 
		// make space for the proteome, names, aminoacidnumber and haming vector
		char **proteome;
		proteome=new char* [nurows];  
		for(int i=0; i<nurows; i++) proteome[i]=new char [rowlength[i]+1]; // +1 is for \0
			
		char **names; 
		names=   new char* [nurows];
		for(int i=0; i<nurows; i++) names[i]=new char [rownamelength[i]+1]; // +1 is for \0 
		
		float *ls; 
		ls=new float [nurows]; 
		
		float **haming; 
	    haming=new float* [nurows]; 
	    for(int i=0; i<nurows; i++) {
			haming[i]=new float [rowlength[i]+mw-1];
			for(int k=0; k<rowlength[i]+mw-1; k++) haming[i][k]=0; // set haming vector to be nulvector
			} 
		
		int **aminoacidnu; // for every row how many is there amino acid A, R, ... 
		aminoacidnu= new int* [nurows];  // Here +5 is needed, otherwise segmentation fault can occur
		for(int i=0; i<nurows; i++) aminoacidnu[i]=new int [20]; 
		
		
		for(int i=0; i<nurows; i++) 
			for(int k=0; k<20; k++) aminoacidnu[i][k]=0;
		 
		
		// open file and save it in an array of strings: proteome, also fill the aminoacidnu array
		file=fopen(argv[1],"r"); 
		if (file==0) printf("Can't open file %s!\n",argv[1]);
		else {
			if (countnames>0) { // fasta format
				int i=-1; 
				int j=0; 
				int l=0; 
				int state=0;
				char c; 
				char cold='\n'; 
				int indx;  
				
				while ( (c=getc(file))!=EOF ) {
					if (c=='>' && cold=='\n') state=1;   
					if (state==1 && c=='\n') {
						state=0; 
						if (i>=0) proteome[i][j]='\0'; 
						names[i+1][l]='\0'; 
						i++; 
						j=0;
						l=0; 
						}
					if (state==1 && c!='\n') {
						names[i+1][l]=c;
						l++;
						} 
					indx=findindex(c); 
					if (state==0 && indx!=-1 && c!='x') {
						proteome[i][j]=c;
						aminoacidnu[i][indx]++;  
						j++;
						}					
					cold=c; 
					}
				proteome[i][j]='\0';
			
			} else { // text format
				char c; 
				int i=0;
				int j=0;
				int indx;   
				while ( (c=getc(file))!=EOF ) { 
					if (c!='\n') {
						indx=findindex(c);
						if (indx!=-1 && c!='x') {
							proteome[i][j]=c;
						    aminoacidnu[i][indx]++;
						    j++; 
						    }
						}
					else {
						proteome[i][j]='\0';
						names[i][0]='\0';  
						i++;
						j=0;  
					}
				} 
			}  
		} 
		fclose(file);  			
		
		
		// read scale and number_of_iterations from command line
		scale=atof(argv[3]);
		noi=atoi(argv[4]); 
		th=atof(argv[5]); 

		
		// make room for position of every amino acid in proteome
		int **pos; 
		pos=new int* [nurows*20]; 
		for (int i=0; i<nurows; i++) 
		     for (int k=0; k<20; k++) pos[i+nurows*k]=new int [aminoacidnu[i][k]+1];
		
		
			
		file=fopen(argv[1],"r"); 
		if (file==0) printf("Can't open file %s!\n",argv[1]);
		else {
			if (countnames>0) { // fasta format
				int hpos[20]; // help positions
				char c;
				char cold='\n'; 
				int state=0; 
				int i=-1;  
				int j=0; 
				int indx; 
				while ( (c=getc(file))!=EOF ) {
					
					if (c=='>' && cold=='\n') { 
						state=1;  
						for(int k=0; k<20; k++) hpos[k]=0;   
						i++; 
						j=0; 			 						 
						} 
	
					if (state==1 && c=='\n') state=0; 
					
					indx=findindex(c); 
					
					if (state==0 && indx!=-1 && c!='x') {
						pos[i+indx*nurows][hpos[indx]]=j; 
						hpos[indx]++; 
						j++; 
					    }
					    
				    cold=c; 
				    }
				 
				
				} else { // text format
					int hpos[20]; // help positions
					int i=0;
					int j=0;
					int indx; 
					char c;  
					for(int k=0; k<20; k++)  hpos[k]=0; 								
					while ( (c=getc(file))!=EOF ) { 
						if (c!='\n') { 
							indx=findindex(c);
							if (indx != -1 && c!='x') {
								pos[i+nurows*indx][hpos[indx]]=j; 
								hpos[indx]++; 
								j++; 
								}
						    } else { 
								i++; 
								for(int k=0; k<20; k++) hpos[k]=0;  
								j=0;  
								}
						 }	
				}
		}
			
		
		
	
		if (ok==1) { // if all sequences of the query are of equal length and the query is properly formated
						
			int mof[mw][20]; // matrix of frequencies
			int nuaa[mw]; // number of amino acids
			int nux[mw]; // number of x (any amino acid symbol)
			int nudaa[mw]; // number of distinct amino acids
			float wm[mw][20]; // weigted model 
			float sw[MNM]; // sequece weights 
			float lom[mw][20]; // logaritmic model
			float lomold[mw][20]; // old logaritmic model
			float comb[mw][20]; // combination model
			float nm[mw][20]; // natural model
			
			float mean; 
			float oldmean=-1000000; 
			float scale2;  
			float threshold;
			 		
		
			
			for (int iter=0; iter<noi; iter++) { // iteration begins here  
	
				for(int i=0; i<mw; i++) // set mof to be zero matrix
					for(int j=0; j<AL; j++) mof[i][j]=0;   

				for(int i=0; i<mw; i++) {nuaa[i]=0; nux[i]=0;} // set number of amino acids and number of x to zero

				for(int i=0; i<mw; i++) // fill the mof with empirical frequencies, and compute nuaa and nux
					for(int l=0; l<nos; l++) {
						int a=findindex(motif[l].sequence[i]); 
						if (a<20) { mof[i][a]=mof[i][a]+1; nuaa[i]++; }
						if (a==20) nux[i]++; 
						}
				
				
				
				for(int i=0; i<mw; i++) { 
					nudaa[i]=0; 
					for(int j=0; j<20; j++) if (mof[i][j]!=0) nudaa[i]++;  // compute the number of distinct amino acids 
					if (nux[i]>0) nudaa[i]=20; // if there is a x in a column then number of distinct amino acids is 20 
					if (nuaa[i]==0) for(int b=0; b<20; b++) mof[i][b]=1;  // what if there is no amino acids in some imput column i.e. it consits only of symbols x 
					}
	
				for(int l=0; l<nos; l++) { // here we compute sequence weights
					float ssw=0; 
					for(int i=0; i<mw; i++) {
						int b=findindex(motif[l].sequence[i]); 
						if (b<20) ssw+=1.0/(float)(nudaa[i]*mof[i][b]); 
						}
					if (ssw==0) ssw=1; 
					sw[l]=ssw; 
					}	
	
				// here we normalize sequence weights
				float ssw=0; 
				for(int l=0; l<nos; l++) ssw=ssw+sw[l]; 
				for(int l=0; l<nos; l++) sw[l]=sw[l]/ssw; 
    
				for(int i=0; i<mw; i++) // set wm to be zero matrix
					for(int j=0; j<20; j++) wm[i][j]=0; 
    
				for(int i=0; i<mw; i++) // compute wm from sw
					for(int l=0; l<nos; l++) { 
						int b=findindex(motif[l].sequence[i]); 
						if (b<20) wm[i][b] +=sw[l]; 
						if (b==20) for(int a=0; a<20; a++) wm[i][a] +=0.05*sw[l]; 
						} 
    
				for(int i=0; i<mw; i++) // add pseudo count
					for(int j=0; j<20; j++) wm[i][j]=(wm[i][j]+0.01/(float)nos)/(1.0+0.2/(float)nos); 
	
				for(int i=0; i<mw; i++)	// set nm to be zero matrix
					for(int j=0; j<20; j++) nm[i][j]=0; 
			
				for(int i=0; i<mw; i++) // compute nm i.e. wm + evolution given by PAM120 matrix
					for(int j=0; j<20; j++) 
						for(int k=0; k<20; k++) nm[i][j]+=wm[i][k]*PAM120[k][j]; 
    
				if (iter==1) for(int i=0; i<mw; i++) // save lom in the first iteration
								for(int j=0; j<20; j++)
									lomold[i][j]=lom[i][j];
				
				float al=f(nos); // compute factor of convex mixing of wm and nm
											
				for(int i=0; i<mw; i++) {
					for(int j=0; j<20; j++) {
						comb[i][j]=(1-al)*wm[i][j]+al*nm[i][j]; // convex combination of wm and nw
					    if (nkeys>0) for(int k=0; k<nkeys; k++) if (keys[k]==i) comb[i][j]=0.01*comb[i][j]+0.99*wm[i][j]; // conserved positions 	
						if (nux[i]>0) comb[i][j]=0.9*background[j]+0.1*comb[i][j]; // if there is x somewhere in the i-th column
						lom[i][j]=log(comb[i][j])-log(background[j]); // compute PSSM lom
						if (iter>0) lom[i][j]=log(0.5*(exp(lom[i][j])+exp(lomold[i][j]))); // combine lom from this and first iteration
						}
					}
				
			
		
				for(int i=0; i<nurows; i++) // set haming vector to zero
					for(int k=0; k<rowlength[i]+mw-1; k++) haming[i][k]=0; 
					
				for(int i=0; i<nurows; i++)	 // compute haming vector using PSSM lom
					for(int j=0; j<20; j++)
						for(int z=0; z<mw; z++) 
							for(int k=0; k<aminoacidnu[i][j]; k++) 
								haming[i][pos[i+j*nurows][k]+mw-1-z]+=lom[z][j]; 
				
				file=fopen("list.txt","w"); // save maxima of haming vectors in a file list.txt so python can be called 
				for(int i=0; i<nurows; i++) {
					float ms=haming[i][mw-1]; 
					for(int j=mw-1; j<rowlength[i]+mw-2; j++)
					if (haming[i][j]>ms) ms=haming[i][j]; 
					fprintf(file,"%f\n", ms); 
					}
				fclose(file); 
				
				char text[3][30];   
	
				file = popen("python logistic.py", "r"); // call python to make a logistic fit
				int it=0; 
				int jt=0; 
				char c; 
				while ( (c=getc(file))!=EOF ) {
					if (c=='\n') {text[it][jt]='\0'; it++; jt=0; } 
					else {text[it][jt]=c; jt++; } 
					}
				pclose(file);
  
				mean=atof(text[0]); // mean paremeter of logistic fit
				scale2=atof(text[1]); // scale paremeter of logistic fit
     
    
				if (mean==oldmean) break; // if there is no change in mean of the logistic distribution break the iteration loop
				oldmean=mean; 
				
				threshold=mean+scale*scale2; // set the threshold
				
				nos=0;
				for(int i=0; i<nurows; i++) { 
					for(int j=mw-1; j<rowlength[i]; j++) {
						if (haming[i][j]>=threshold) { 
							motif[nos].Evalue=(0.5 - 0.5 * tanh((haming[i][j]-mean)/(2*scale2)))*nurows; 
							
							for(int ii=0; ii<mw; ii++) motif[nos].sequence[ii]=proteome[i][j-mw+1+ii];  
							motif[nos].sequence[mw]='\0'; 
							motif[nos].row=i+1;  
							motif[nos].rowpos=j-(mw-1)+1; 
							nos++;  
							if (nos>MNM) {
								nos=MNM; 
								iter=iter+noi; // end iteration by adding noi = number of iterations
								}
							}
						}
					}
				if (nos==0) iter=iter+noi; // end iteration 
				else {
					if (iter<=4) {
		
						float mof2[mw][20];  // matrix of frequencies 2
						float lom2[mw][20];  // logaritmic model 2
			 		
	
						for(int i=0; i<mw; i++) // set mof2 to be zero matrix
							for(int j=0; j<AL; j++) mof2[i][j]=0;   


						for(int i=0; i<mw; i++) { // fill the mof2 with empirical frequencies
							for(int l=0; l<nos; l++) {
								int a=findindex(motif[l].sequence[i]); 
									mof2[i][a]=mof2[i][a]+1;  
								}
							}
				
						for(int i=0; i<mw; i++) // add pseudo count
							for(int j=0; j<20; j++) {
								mof2[i][j]=mof2[i][j]/(float)nos; 
								mof2[i][j]=(mof2[i][j]+0.01/(float)nos)/(1.0+0.2/(float)nos); 
							}				
				
						for(int i=0; i<mw; i++) {
							for(int j=0; j<20; j++) {
								lom2[i][j]=log(mof2[i][j])-log(background[j]); // compute PSSM lom
								}
							}
						for(int l=0; l<nos; l++) {
							motif[l].postScore=0; 
							for(int i=0; i<mw; i++) 
								motif[l].postScore=motif[l].postScore+lom2[i][findindex(motif[l].sequence[i])];  
								motif[l].postScore=exp(motif[l].postScore); 
							}
		
						
						long double list[nos];
						long double listmean=0;   
						
						
						for(int l=0; l<nos; l++) {
							list[l]=motif[l].postScore; 
							listmean=listmean+list[l]/(float)nos; 								
							}
						
						float temp; 
						for(int l1=0; l1<nos; l1++)
							for(int l2=l1+1; l2<nos; l2++)
								if (list[l1]>list[l2]) {
									temp=list[l2]; 
									list[l2]=list[l1]; 
									list[l1]=temp; 
									}

						int din=nos/4; // down index 
						float down=0; 
						if (din>0) down=list[din]; 
						
  
						
						int newnos=0; 
						for(int l=0; l<nos; l++) {
							if ( ( motif[l].postScore>=th*listmean) || ( motif[l].postScore>=down ) ) {
								for(int ii=0; ii<=mw; ii++) motif[newnos].sequence[ii]=motif[l].sequence[ii]; 
								motif[newnos].postScore=motif[l].postScore; 
								motif[newnos].row=motif[l].row; 
								motif[newnos].rowpos=motif[l].rowpos; 
								motif[newnos].Evalue=motif[l].Evalue; 
								newnos++; 
								} 
							}
						
 	

						printf("Thrown out %d\n", nos-newnos);  
						nos=newnos; 
						}// iter<4
					}
				} // end iteration for loop 
			
			
			if (nos<MNM) printf("Final number of sequences: %d  (scale=%f, number_of_iterations=%d)\n",nos,scale,noi);
			if (nos>=MNM) printf("Maximal output size reached. Try increasing scale!\n");
			for(int i=0; i<nos; i++) {
				printf("sequence%d\t%d-%d\t(Evalue = %f)\n", motif[i].row, motif[i].rowpos,motif[i].rowpos+mw-1,motif[i].Evalue);
				if (names[motif[i].row-1][0]!='\0') printf("%s \n", names[motif[i].row-1]);  
				printf("%s\n", motif[i].sequence);
				}
			
			system("rm list.txt"); 
				
			} else printf("All sequences of the query must be of equal length. Query can consist only of amino acids ARNDCQEGHILKMFPSTWYV and symbol X or x.\n"); 
				
		} else printf("Four arguments expected: proteome query scale number_of_iterations\n"); 
	
	return 0; 
	}; 
