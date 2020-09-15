/*#############################################################################
  
  # BayClone_C : 
  # Copyright (C) <2018>  <Subhajit Sengupta>
  # This file is part of BayClone_C
  # This is a helper file for BayClone_C
  # This file is responsible for parsing input from VCF and CN file

  # by Subhajit Sengupta (subhajit@uchicago.edu); April, 2017  
  
  # This file is used to get input data for CN and allele count 
  # from Clonal CN consensus Sanger format input and VCf file

  #  BayClone_C is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

  #  BayClone_C is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  #  along with BayClone_C. If not, see <http://www.gnu.org/licenses/>.
  
#############################################################################
# To run this module:
#  ./parseInputData <input vcf file> <input CN file> <output file>
#############################################################################*/

#include "stdio.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define MAX_LINE_LEN	50000
#define MAX_FLD_LEN		500


// Sanger format VCF file
//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOUR
//#FORMAT = GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM

#define VCF_CHR_NAME_INDX	0
#define VCF_POS_INDX	1
#define VCF_REF_INDX	3	
#define VCF_ALT_INDX	4	
#define VCF_PASS_INDX	6
#define VCF_INFO_INDX	7
#define VCF_FORMAT_INDX		8
#define VCF_SAMPLE_OFFSET	9
#define VCF_NORMAL_INDX	9
#define VCF_TUMOR_INDX	10

/// Clonal consensus CN calls
#define CC_CHR_NAME_INDX	0
#define CC_CHR_POS_STRT_INDX 1
#define CC_CHR_POS_END_INDX 2
#define CC_CN_INDX	3
#define CC_MAJ_INDX	4
#define CC_MIN_INDX	5
#define	CC_CFRAC_INDX	6


#define BB_CHR_NAME_INDX	0
#define BB_CHR_POS_STRT_INDX 1
#define BB_CHR_POS_END_INDX 2
#define BB_MAJ1A_INDX	7
#define BB_MIN1A_INDX	8
#define	BB_FRAC1A_INDX	9
#define BB_MAJ2A_INDX	10
#define BB_MIN2A_INDX	11
#define	BB_FRAC2A_INDX	12

typedef struct _SNV_CN
{
	// location info	
	char chrName[10];
	int startPos;
	int endPos;

	// from mutect vcf file
	int N_tot_T; // tumor total cnt
	int n_alt_T; // tumor alt cnt
	int N_tot_N; // normal total cnt 
	int n_alt_N; // tumor alt cnt

	// from BB file
	int nMaj1_A;
	int nMin1_A;
	int totCN;
	double frac1_A;
	int nMaj2_A;
	int nMin2_A;
	double frac2_A;
	double sampleCN; // our input M for BayClone2
	double maxAllelCnt;	

}SNV_CN;


int readNextLine (char *l, FILE *f)
{
	if (fgets (l, MAX_LINE_LEN, f) == NULL) 
	{
		/* NULL.  EOF was read.  reset l, return */
		l[0] = '\0';
	   	return 1;
	}
	else
	{ 
		if (strlen (l) == MAX_LINE_LEN -1) 
		{
			/* Line too long - do some malloc stuff.
	 		 * For now, just terminate and return. */
 	 		//error("Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
 	 		fprintf(stdout,"Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
	     	l[0] = '\0';
	     	return 2;
		}
		else    
		{
			/* Proper line read. Return 0 */
	    	return 0;
	   	}
	}
}

int file_line_count(char *file_name)
{
	char c; 
	int lines = 0;
	FILE *fp;
	//fprintf(stdout,"fileName: %s\n",file_name);
	fp = fopen(file_name,"r");
	if(NULL==fp)
	{	
		//fprintf(stdout,"BB file open error !!\n");
		return 0;
	}
	//fprintf(stdout,"reading files...\n");
	while ((c = fgetc(fp)) != EOF)
	{
		if (c == '\n') 
			++lines;
	}

	fclose(fp);
	return lines;
}

int export_N_n(char* vcf_name, SNV_CN **SC1, int* VCF_n_rcd)
{

	//Rprintf("%s:%s:%s\n",*R_vcf_name,*R_N_tot_file,*R_n_alt_file);

	//char vcf_name[100];
	//char N_tot_file[100];
	//char n_alt_file[100];

	FILE* fp;
	FILE* ft1,* ft2;
	char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];
	char* pch,*pch1;
	//char *pch2;
	int cnt1,s_cnt;
	//int GT_fld;
	//int ALT_CNT_fld;
	//int REF_CNT_fld;
	char fmt_fld[MAX_FLD_LEN];
	char pass_str[MAX_FLD_LEN];
	int i,j;
	int e_flag;
	//int N_tot;
	int n_alt,N_minus_n;
	char N_n_cnt_str[MAX_FLD_LEN];
	int N_tot_samples,n_alt_samples;
	int total_line_cnt;
	int header_line_cnt = 0;
	int snv_line_cnt;

	char chrName[10];
	int posSNV;
	char refStr[10];
	char altStr[10];
	
	int refID,altID;
	int fCnt[4],rCnt[4];
	char nucleo_arr[4] = {'A','C','G','T'};

	char* pass_pch;

    char delimit[]=" \t;";

	//strcpy(vcf_name,R_vcf_name);
	//strcpy(N_tot_file,R_N_tot_file);
	//strcpy(n_alt_file,R_n_alt_file);

	total_line_cnt = file_line_count(vcf_name);
	fprintf(stdout,"%d lines\n",total_line_cnt);

	fp = fopen(vcf_name,"r");		

	//ft1=fopen(N_tot_file,"w");
	//ft2=fopen(n_alt_file,"w");
	
	int vcf_line_cnt = 0;
	if(NULL==fp)
		//||(NULL==ft1)||(NULL==ft2))
	{
		//error("file open error !\n");
		fprintf(stdout,"VCF file open error !!\n");
		return -1;
	}
	else
	{
		while (readNextLine(vcf_line, fp) == 0)
		{
    		vcf_line_cnt++;
			if(vcf_line_cnt == 1)
			{		
				snv_line_cnt = total_line_cnt - 1;
			 	*SC1 = (SNV_CN*)(calloc(snv_line_cnt,sizeof(SNV_CN)));
			}
			else
			{
				if(vcf_line[0] != '#')
				{
					e_flag = 0;
					
 					cnt1=0;                 
					pch = strtok(vcf_line,delimit);
	                n_alt = 0;
					N_minus_n = 0;
					double vaf = 0.0;
					while (pch != NULL)
	                {
						if(cnt1 == VCF_CHR_NAME_INDX)
							strcpy(chrName,pch);
						 
						if(cnt1 == VCF_POS_INDX)
							sscanf(pch,"%d",&posSNV);	
						
						// for alt count
						if(strstr(pch, "t_alt_count") != NULL)
						{
							//fprintf(stdout,"%s\n",pch);
							sscanf(pch,"%*[^=]=%d",&n_alt);
							
						}
						
						// for ref count		
						if(strstr(pch, "t_ref_count") != NULL)
						{
							//fprintf(stdout,"%s\n",pch);
							sscanf(pch,"%*[^=]=%d",&N_minus_n);
						}

						if(strstr(pch, "VAF") != NULL)
						{
							//fprintf(stdout,"%s\n",pch);
							sscanf(pch,"%*[^=]=%lf",&vaf);
						}
											
						pch = strtok (NULL, delimit);
						cnt1++;
					}	


					N_tot_samples = n_alt+N_minus_n;
					n_alt_samples = n_alt;
					if((N_tot_samples == 0)||(vaf==0.0))
						e_flag = 1;
					//printf(stdout,"%d\t%d\n",N_tot_samples,e_flag);	
					if(e_flag==0)
					{
						int k = (*VCF_n_rcd);
						//fprintf(stderr,"%s\n",chrName);
						strcpy((*SC1)[k].chrName,chrName);
						(*SC1)[k].startPos = posSNV;
						(*SC1)[k].endPos = posSNV;
						(*SC1)[k].N_tot_T = N_tot_samples;
						(*SC1)[k].n_alt_T = n_alt_samples;
						(*SC1)[k].N_tot_N = 0;
						(*SC1)[k].n_alt_N = 0;

						(*VCF_n_rcd)++;		
					}
				}
			}
		}
	
	}// else file open success

	fclose(fp);
	//fclose(ft1);
	//fclose(ft2);
	return 0;
	//Rprintf("files are generated!\n");
	//fprintf(stdout,"files are generated!\n");

}
int get_sample_copy_val_for_SNVs(char* BB_file_name,SNV_CN** SC2,int* BB_n_rcd)
{
	FILE* fp;
	char BB_line[MAX_LINE_LEN];
	//fprintf(stdout,"BB file Name: %s\n",BB_file_name);	

	fp = fopen(BB_file_name,"r");
	if(NULL == fp)
	{
		fprintf(stdout,"BB input file open error!!\n");
		return -1;
	}
	int total_line_cnt = file_line_count(BB_file_name);
	//fprintf(stdout,"BB line cnt: %d\n",total_line_cnt);
	
	*SC2 = (SNV_CN*)(calloc(total_line_cnt,sizeof(SNV_CN)));
	int cnt1,i=0;
	char* pch;
	int nLine = 0;
	while (readNextLine(BB_line, fp) == 0)
	{
			nLine++;
			//fprintf(stdout,"%s\n",BB_line);
			/* If it is a comment, keep track as header line */
	     	if (BB_line[0] == '#' || nLine == 1)
				continue;
			else
			{
				cnt1=0;
				pch = strtok(BB_line,"\t");
				while (pch != NULL)
				{
					//Rprintf("%s\n",pch);
					switch(cnt1)
					{
						case BB_CHR_NAME_INDX:
							strcpy((*SC2)[i].chrName,pch);break;
						case BB_CHR_POS_STRT_INDX:
							sscanf(pch,"%d",&((*SC2)[i].startPos));break;
						case BB_CHR_POS_END_INDX:
							sscanf(pch,"%d",&((*SC2)[i].endPos));break;
						case BB_MAJ1A_INDX:
							sscanf(pch,"%d",&((*SC2)[i].nMaj1_A));break;
						case BB_MIN1A_INDX:
							sscanf(pch,"%d",&((*SC2)[i].nMin1_A));break;
						case BB_FRAC1A_INDX:
							sscanf(pch,"%lf",&((*SC2)[i].frac1_A));break;
						//default: fprintf(stdout,"BB file parsing error !!\n");break;	
					}
					if((*SC2)[i].frac1_A < 1.0)
					{
						switch(cnt1)
						{	
							case BB_MAJ2A_INDX:
								sscanf(pch,"%d",&((*SC2)[i].nMaj2_A));break;
							case BB_MIN2A_INDX:
								sscanf(pch,"%d",&((*SC2)[i].nMin2_A));break;
							case BB_FRAC2A_INDX:
								sscanf(pch,"%lf",&((*SC2)[i].frac2_A));break;
						}
					}			
					else
					{
						(*SC2)[i].nMaj2_A = 0;
						(*SC2)[i].nMin2_A = 0;
						(*SC2)[i].frac2_A = 0.0;
							
					}
					pch = strtok (NULL, "\t");
					cnt1++;
				}
				i++;
			}

	}
	
	*BB_n_rcd = i;
	
	/*
	for (int j=0;j<(*BB_n_rcd);j++)
	{
			fprintf(stdout,"%s:%d-%d\t", (*SC2)[j].chrName,(*SC2)[j].startPos,(*SC2)[j].endPos);

			fprintf(stdout,"(%d+%d)*%lf+(%d+%d)*%lf = ",SC2[j].nMaj1_A,SC2[j].nMin1_A,SC2[j].frac1_A,
													  SC2[j].nMaj2_A,SC2[j].nMin2_A,SC2[j].frac2_A);
			SC2[j].sampleCN =  (SC2[j].nMaj1_A+SC2[j].nMin1_A)*SC2[j].frac1_A + 
							(SC2[j].nMaj2_A+SC2[j].nMin2_A)*SC2[j].frac2_A; 
			fprintf(stdout,"%lf\n",SC2[j].sampleCN);					

	}	
	*/
	fclose(fp);
	return 0;
}


int main(int argc,char* argv[])
{
	SNV_CN *SC1, *SC2;
	int BB_n_rcd = 0;
	int VCF_n_rcd = 0;
	FILE* ft;

	if(argc < 4)
	{
		fprintf(stdout,"Please enter: <prog_name> <input vcf file> <input CCC file> <output file>\n");
		return -1;
	}	
	else
	{	
		ft = fopen(argv[3],"w");
		if(NULL == ft)
		{
			fprintf(stdout,"output file open error!!\n");
			return -1;
		}		
		
		int iRet = get_sample_copy_val_for_SNVs(argv[2],&SC2,&BB_n_rcd);
		if(iRet == -1)
			return -1;
		
		fprintf(stdout,"%d BB records\n",BB_n_rcd);
		
		for (int j=0;j<BB_n_rcd;j++)
		{
			
			//fprintf(stdout,"%s:%d-%d\t", SC2[j].chrName,SC2[j].startPos,SC2[j].endPos);
			
			//fprintf(stdout,"(%d+%d)*%lf+(%d+%d)*%lf = ",SC2[j].nMaj1_A,SC2[j].nMin1_A,SC2[j].frac1_A,
			//										  SC2[j].nMaj2_A,SC2[j].nMin2_A,SC2[j].frac2_A);
			SC2[j].sampleCN =  (SC2[j].nMaj1_A+SC2[j].nMin1_A)*SC2[j].frac1_A + 
							   (SC2[j].nMaj2_A+SC2[j].nMin2_A)*SC2[j].frac2_A; 
			///actually this involves CCF of that CNA region; so this is sampleCN in the tumor fraction
			//fprintf(stdout,"%lf\n",SC2[j].sampleCN);					
			SC2[j].maxAllelCnt =  SC2[j].nMaj1_A*SC2[j].frac1_A + SC2[j].nMaj2_A*SC2[j].frac2_A;
			if(SC2[j].nMin1_A*SC2[j].frac1_A + SC2[j].nMin2_A*SC2[j].frac2_A > SC2[j].maxAllelCnt)
				SC2[j].maxAllelCnt = SC2[j].nMin1_A*SC2[j].frac1_A + SC2[j].nMin2_A*SC2[j].frac2_A; 
			
		}	
		

		iRet = export_N_n(argv[1],&SC1,&VCF_n_rcd);
		if(iRet == -1)
			return -1;
		
		fprintf(stdout,"%d VCF records\n",VCF_n_rcd);

		/*
		for(int j=0;j<VCF_n_rcd;j++)
		{	
			fprintf(stdout,"%s:%d-%d\t",SC1[j].chrName,SC1[j].startPos,SC1[j].endPos);	
			fprintf(stdout,"%d:%d\t%d:%d\n",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);		
		}
		*/

		//fprintf(ft,"#chr\tpos\ttumor_N\ttumor_n\tsampleCN_tumor\n");
		int bFound = 0;
		for(int j=0;j<VCF_n_rcd;j++)
		{
			SC1[j].sampleCN = 2.0;
			bFound = 0;	
			//SC1[j].sampleCN = -1.0; // when only clonal CN info is avalable
			int loc1 = SC1[j].startPos;
			for(int k=0;k<BB_n_rcd;k++)
			{
				if( (strcmp(SC1[j].chrName,SC2[k].chrName)==0) && (loc1 >= SC2[k].startPos) && (loc1 <= SC2[k].endPos) )
				{
					SC1[j].sampleCN = SC2[k].sampleCN;
					SC1[j].maxAllelCnt = SC2[k].maxAllelCnt;

					SC1[j].nMaj1_A = SC2[k].nMaj1_A;
					SC1[j].nMin1_A = SC2[k].nMin1_A;
					SC1[j].frac1_A = SC2[k].frac1_A;

					SC1[j].nMaj2_A = SC2[k].nMaj2_A;
					SC1[j].nMin2_A = SC2[k].nMin2_A;
					SC1[j].frac2_A = SC2[k].frac2_A;
					bFound = 1;	
				}	
			}
			//fprintf(stdout,"%s:%d\t",SC1[j].chrName,SC1[j].startPos);	
			//fprintf(stdout,"%d\t%d\t%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);
			//fprintf(stdout,"%lf\n",SC1[j].sampleCN);
			if((bFound ==1) && (SC1[j].sampleCN > 0))
			{
				fprintf(ft,"%s\t%d\t",SC1[j].chrName,SC1[j].startPos);	
				//fprintf(ft,"%d\t%d\t%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);
				fprintf(ft,"%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T);

				fprintf(ft,"%lf\t",SC1[j].sampleCN);
				fprintf(ft,"%d\t%d\t%lf\t",SC1[j].nMaj1_A,SC1[j].nMin1_A,SC1[j].frac1_A); 
				if(SC1[j].frac1_A < 1.0){				
					fprintf(ft,"%d\t%d\t%lf\n",SC1[j].nMaj2_A,SC1[j].nMin2_A,SC1[j].frac2_A);
				}else{
					fprintf(ft,"NA\tNA\tNA\n");
				}
				//fprintf(stdout,"%s\t%d done!!\n",SC1[j].chrName,SC1[j].startPos); 
				//fprintf(ft,"%lf\n",SC1[j].maxAllelCnt);
			}
				
		}		

		fprintf(stdout,"%s outout file is written !!\n",argv[3]);
	}


	fclose(ft);

	free(SC1);
	free(SC2);

	return 0;
}

