#############################################################################
  
  # BayClone_C : 
  # Copyright (C) <2018>  <Subhajit Sengupta>

  # This file is part of BayClone_C
  # This is the main R file for BayClone_C
  # This file is responsible for running the algorithm and producing the
  # output in the desired file format specified by PCAWG (consensus)

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
  
############################################################################## 
# Requires: mclust, fpc package
# Please install those packages before running this code
######################################################################################
# To run this module:
#  R CMD BATCH --no-save --no-restore '--args 
#   <input_file_created_in_previous_step> <purity_ploidy_file>' BayCloneC.R run.log
######################################################################################

rm(list=ls(all=TRUE))


is.nonint <-function(a){
    return(!(a==round(a)))
}

clonality_test_for_sample <-function(data,pur,config_param){
    
    
    TOTAL_RD_INDX = config_param$TOTAL_RD_INDX
    VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
    SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
    MAJOR1_CN_INDX = config_param$MAJOR1_CN_INDX
    MINOR1_CN_INDX = config_param$MINOR1_CN_INDX
    
    res = NULL
    
    if(any(is.nonint(data[,SAMPLE_T_CN_INDX]))){
        res$flag = 0
    }else{
        res$flag = 1
    }
        
    VAF = data[,VARIANT_RD_INDX]/data[,TOTAL_RD_INDX]

    pur_adjust = 2*(1-pur)/pur
   
    Z_arr = rep(0,length(VAF))
        
    for(j in 1:length(VAF)){
        max_allele_cnt = max(data[j,MAJOR1_CN_INDX],data[j,MINOR1_CN_INDX])
        z_val_possible = 1:max_allele_cnt
        #cat(data[j,SAMPLE_T_CN_INDX],"\t",pur_adjust,"\t",z_val_possible,"\t",VAF[j],"\n")
        Z_vec_arr = abs(VAF[j] - z_val_possible/(data[j,SAMPLE_T_CN_INDX]+pur_adjust))
        Z_arr[j] = z_val_possible[which.min(Z_vec_arr)]
    }
        
    E_VAF = Z_arr/(data[,SAMPLE_T_CN_INDX]+pur_adjust)
        
    residual_VAF = VAF - E_VAF
    sd_res_VAF = NULL
    
    if(length(residual_VAF) > 1)
        sd_res_VAF = sd(residual_VAF)
    
    res$sd = sd_res_VAF
    
    if(res$flag == 0 | is.null(sd_res_VAF)){
        return(res)
    }
        
    std_res_VAF = (residual_VAF-mean(residual_VAF))/sd_res_VAF
    
    MC_std_res_VAF = Mclust(std_res_VAF)

    T_cluster_from_MC = as.vector(table(MC_std_res_VAF$classification))

    c_mean = MC_std_res_VAF$parameters$mean
    
    cnt1 = 0
    S = sum(T_cluster_from_MC)
    
    del.ct = NULL
    for(k in 1:length(T_cluster_from_MC)){
        if(T_cluster_from_MC[k]/S < 0.05){  ### smaller cluster size
            cnt1 = cnt1+1
            del.ct = c(del.ct,k)
        } 
    }
    
    if(!is.null(del.ct)){
        c_mean = c_mean[-del.ct]
    }
    
    nClust = dim(MC_std_res_VAF$z)[2] -  cnt1
        
    if(nClust > 1){
        data_sd = sd(std_res_VAF)
        cat("***HERE***\n")
        if(sd(c_mean) < 0.5*data_sd)  ### if standard deviation of cluster mean is less than the half of the data std then it is clonal
            nClust = 1
    }
        
    if(nClust > 1){
        res$flag = 0
    }else{
        res$flag = 1
    }
    
    return(res)
}

find_ccf_for_each_SNV <-function(data,pur,consensus_ploidy,clonal_flag,dir_path,config_param){
    
    CHR_NAME_INDX = config_param$CHR_NAME_INDX
    CHR_POS_INDX = config_param$CHR_POS_INDX
    
    TOTAL_RD_INDX = config_param$TOTAL_RD_INDX
    VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
    SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
    MAJOR1_CN_INDX = config_param$MAJOR1_CN_INDX
    MINOR1_CN_INDX = config_param$MINOR1_CN_INDX
    FRAC1_INDX = config_param$FRAC1_INDX
    MAJOR2_CN_INDX = config_param$MAJOR2_CN_INDX
    MINOR2_CN_INDX = config_param$MINOR2_CN_INDX
    FRAC2_INDX = config_param$FRAC2_INDX
    
    MAX_CN_CUTOFF = config_param$MAX_CN_CUTOFF
    PLOIDY_MAX_CN_CUTOFF = config_param$PLOIDY_MAX_CN_CUTOFF
        
    
    # f1 = sprintf("./%s/%s_mult_ccf.dat",dir_path,sample_id)
    
    VAF = data[,VARIANT_RD_INDX]/data[,TOTAL_RD_INDX]
    
    vaf_adj = (pur*data[,SAMPLE_T_CN_INDX] + (1-pur)*2)/pur
    
    #### mutation CN
    mutation_CN_for_SNV = VAF*vaf_adj
    
    nSNV = length(data[,VARIANT_RD_INDX])
    
    mult = rep(0,nSNV)
    ccf = rep(0,nSNV)
        
    ##### this part is not needed for PCAWG @Jan21 
    ##### calculates ploidy from sample_avg_tumor_CN #####
    id = which(data[,SAMPLE_T_CN_INDX] < PLOIDY_MAX_CN_CUTOFF)
    S = length(id)
    ploidy = sum(data[id,SAMPLE_T_CN_INDX] )/S
    
    ##############
    
    
    for(i in 1:nSNV){
        bFlag = 0
        mult_option = NULL
        if(data[i,FRAC1_INDX] == 1.0){
            m1 = max(data[i,MAJOR1_CN_INDX],data[i,MINOR1_CN_INDX])
            mult_space = seq(1,m1)
        }else{
            m1 = max(data[i,MAJOR1_CN_INDX],data[i,MINOR1_CN_INDX])
            m2 = max(data[i,MAJOR2_CN_INDX],data[i,MINOR2_CN_INDX])
            
            set1 = seq(0,m1)*data[i,FRAC1_INDX]
            set2 = seq(0,m2)*data[i,FRAC2_INDX]
            
            v0 = as.vector(outer(set1,set2,'+'))
            mult_space = v0[v0!=0]
            bFlag = 1 
        }
        
        if(mutation_CN_for_SNV[i] < MAX_CN_CUTOFF)
        {
            dist_CCF_to_one = abs(1.0-(mutation_CN_for_SNV[i]/mult_space))
            
            #dist_CCF_to_one = abs(CCF_init_est[i] - mult_option)
            
            idx = which.min(dist_CCF_to_one)
            mult[i] = mult_space[idx]
        }else{
            #mult[i] = floor(CCF_init_est[i])
            mult[i] = round(mutation_CN_for_SNV[i])
        }
        
        if(clonal_flag == 1){
            ccf[i] = 1.0
        }else{
            ccf[i] = mutation_CN_for_SNV[i]/mult[i]
        	if(bFlag == 1){
				if(mult[i] %in% set1){
					### x+0 type configuration then CCF should be lower than 1.0
					ccf[i] = ccf[i]*data[i,FRAC1_INDX]
				}
				if(mult[i] %in% set2){
					### 0+x type configuration then CCF should be lower than 1.0
					ccf[i] = ccf[i]*data[i,FRAC2_INDX]
				}
			}
		}
    }

    data1 = cbind(data,mutation_CN_for_SNV,mult,ccf)
    names(data1) = c("chr","pos","N","n","M_B","nMaj1_A","nMin1_A","frac1_A","nMaj2_A","nMin2_A","frac2_A","mutation_CN","multiplicity","ccf")
    # write.table(data1,row.names=F,sep="\t",quote=F,file=f1)

    ###############
    
    f2 = sprintf("./%s/%s_multiplicity.txt",dir_path,sample_id)
    chr_name = as.character(data[,CHR_NAME_INDX])
    chr_loc = data[,CHR_POS_INDX]
    tumor_CN = data[,SAMPLE_T_CN_INDX]
    
    col4=rep(NA,length(mult))
    col5=rep(NA,length(mult))
    T1 = cbind(chr_name,chr_loc,tumor_CN,mult,col4,col5)
    write.table(T1,col.names=c("chr","pos","tumour_copynumber","multiplicity","multiplicity_options","probabilities"),row.names=F,sep="\t",quote=F,file=f2)
    
    
    f3 = sprintf("./%s/%s_purity_ploidy.txt",dir_path,sample_id)
    T1 = cbind(sample_id,pur,consensus_ploidy)  ### ploidy is our estimate
    write.table(T1,col.names=F,row.names=F,sep="\t",quote=F,file=f3)
    
    res = NULL
    res$ccf = ccf
    res$ploidy = ploidy
    
    return(res)
}

find_nclust_ccf_of_SNV <- function(ccf,clonal_flag,config_param){
    
    MAX_CLUSTER = config_param$MAX_CLUSTER  ### 7 for pcawg
    MIN_CLUST_DIST = config_param$MIN_CLUST_DIST
    nClust = -1
    if(clonal_flag == 1){            
        #nClust = 0
		CCF_Info$mean_ccf = 1.0
		CCF_Info$old_label = 1
		CCF_Info$mut_assign_vec = rep(1,length(ccf))
		CCF_Info$old_size = length(ccf)
    }else{
        #MC=Mclust(ccf,1:MAX_CLUSTER)

        MC_BIC = mclustBIC(ccf, 1:MAX_CLUSTER)
        S_MC = summary(MC_BIC, ccf)
		CCF_Info = NULL
		
		if(S_MC$G > 1){
			merged_MC = mergenormals(ccf,S_MC,method="ridge.uni")
       
        	#uniq_class = sort(unique(merged_MC$clustering))
        	#uniq_class_len =  length(uniq_class)

        	#cat("means = ",MC$parameters$mean,"\n")
        	#cat("assignment 1 = ",MC$classification[1:20],"\n")

        	#mean_vec = rep(0,uniq_class_len)
        	#cnt_id = 0 
        	#for(m_id in uniq_class){
            #	cnt_id = cnt_id+1
           # 	data_id = which(MC$classification == m_id)
            #	mean_vec[cnt_id] = MC$parameters$mean[m_id]
            #	MC$classification[data_id] = cnt_id
        	#}   
        	#cat("median = ",median_from_data,"\n")

        	#cat("assignment 2 = ",MC$classification[1:20],"\n")
        	#nMC = uniq_class_len

        #CCF_Info = NULL

        #CCF_Info$mean_ccf =  MC$parameters$mean
        #CCF_Info$mean_ccf = median_from_data
        #CCF_Info$mean_ccf = mean_vec
        #cat("mean_ccf from MC = ",MC$parameters$mean,"\n")
        #cat("median from data = ",median_from_data,"\n")
        #CCF_Info$old_label = 1:nMC
        #CCF_Info$old_label = 1:cnt_id

 
		
			cat("all means ",merged_MC$muarray,"\n")
			cat("kept cluster", merged_MC$clusternumbers,"\n")

			selected_mean_ccf = merged_MC$muarray[merged_MC$clusternumbers]
			cat("new mean_ccf ",selected_mean_ccf,"\n")

			uniq_class = sort(unique(merged_MC$clustering))
        	uniq_class_len =  length(uniq_class)
			cat("uniq_class ", uniq_class, "\n")
	
			filtered_selected_mean_ccf = merged_MC$muarray[merged_MC$clusternumbers[uniq_class]]
			cat("filtered_selected_mean_ccf = ",filtered_selected_mean_ccf,"\n")
	
        	mean_vec = rep(0,uniq_class_len)
        	cnt_id = 0 
        	for(m_id in uniq_class){
            	cnt_id = cnt_id+1
            	data_id = which(merged_MC$clustering == m_id)
            	mean_vec[cnt_id] = filtered_selected_mean_ccf[m_id]
        	}   

			CCF_Info$old_label = 1:cnt_id
        	CCF_Info$mean_ccf = mean_vec
			CCF_Info$mut_assign_vec = merged_MC$clustering
			CCF_Info$old_size = as.vector(table(merged_MC$clustering))
		
		}else{
			CCF_Info$mean_ccf = 1.0
			CCF_Info$old_label = 1
			CCF_Info$mut_assign_vec = rep(1,length(ccf)) 
			CCF_Info$old_size = length(ccf)	
		}
	
	}
	
	return(CCF_Info) ### return nClust-1
}

find_cluster_center_by_kmean <- function(data,ccf,clonal_flag,pur,nc_to_run_with,dir_path,config_param,CCF_Info){
    
    cat("nclust for KK = ",nc_to_run_with,"\n")
    
    CHR_NAME_INDX = config_param$CHR_NAME_INDX
    CHR_POS_INDX = config_param$CHR_POS_INDX
    
    TOTAL_RD_INDX = config_param$TOTAL_RD_INDX
    VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
    
    SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
    
    MAJOR1_CN_INDX = config_param$MAJOR1_CN_INDX
    MINOR1_CN_INDX = config_param$MINOR1_CN_INDX
    FRAC1_INDX = config_param$FRAC1_INDX
    
    MAJOR2_CN_INDX = config_param$MAJOR2_CN_INDX
    MINOR2_CN_INDX = config_param$MINOR2_CN_INDX
    FRAC2_INDX = config_param$FRAC2_INDX
    KMEAN_CENTER_CUTOFF = config_param$KMEAN_CENTER_CUTOFF
    
    #nc_to_run_with = nc_mclust[i,2]+1 ### add one (clonal population)
    
    new_number_of_cluster = nc_to_run_with-1
   
    KK=kmeans(ccf,nc_to_run_with,nstart=200)  ### CHK !!

    mean_ccf = KK$centers
    old_label = sort(unique(KK$cluster))
    
    mut_assign_vec = KK$cluster
    old_size = KK$size
    
    I0=which(KK$centers >= 1-((1-KMEAN_CENTER_CUTOFF)/pur))
    #I0=which(KK$centers >= KMEAN_CENTER_CUTOFF) ### or 0.95 ?
    cat("KK_centers = ",KK$centers,"\n")
    cat("I0 = ",I0,"\n")
    
    mean_ccf[I0] = 1.0
    I1=old_label[I0]
    cat("I1 = ",I1,"\n")
    cat("old_label = ",old_label,"\n")
    temp_label = setdiff(old_label,I1)
    cat("temp_label = ",temp_label,"\n")
    
    temp_size = old_size[temp_label]
    
    cat("temp_size",temp_size,"\n")
    cat("old_size",old_size,"\n")
    new_mean_ccf = mean_ccf
    new_mut_assign_vec = mut_assign_vec
    new_size = old_size
    new_label = old_label ###
    
    L1 = length(I1)
    if(L1>1){
        
        min_id1 = min(I1)
        size_min_id1 = sum(old_size[I1])
        
        snv_indx_to_change = which(mut_assign_vec %in% I1 )
        new_mut_assign_vec[snv_indx_to_change] = min_id1
        
        new_label = c(temp_label,min_id1)    
        
        new_size = c(temp_size,size_min_id1)
        
        new_mean_ccf = mean_ccf[new_label]
        
        new_number_of_cluster = length(new_label) - 1
    
    
        cat("new label = ",new_label,"\n")
        cat("new_size = ",new_size,"\n")
        cat("new_number_of_cluster = ",new_number_of_cluster,"\n")
        
        #sorted_label = sort(new_mean_ccf,index.return=T)$ix
        sorted_label = rank(new_mean_ccf)
        sorted_mut_assign_vec = new_mut_assign_vec
        
        for(j in 1:length(new_label)){
            mut_id = which(new_mut_assign_vec == new_label[j])
            sorted_mut_assign_vec[mut_id] = sorted_label[j] 
        }
        cat("new_mean_ccf = ",new_mean_ccf,"\n")
        cat("sorted_label = ",sorted_label,"\n")
        
        display_id = sort(new_mean_ccf,index.return=T)$ix
        cat("display_id = ",display_id,"\n")
    
    }else{
        sorted_label = rank(new_mean_ccf)
        sorted_mut_assign_vec = new_mut_assign_vec
        for(j in 1:length(new_label)){
            mut_id = which(new_mut_assign_vec == new_label[j])
            sorted_mut_assign_vec[mut_id] = sorted_label[j] 
        }
        display_id = sort(new_mean_ccf,index.return=T)$ix
        cat("display_id = ",display_id,"\n")
        cat("else part\n")
    }
    
    f2 = sprintf("./%s/%s_subclonal_structure.txt",dir_path,sample_id)
    T2 = cbind(sorted_label[display_id],new_size[display_id],new_mean_ccf[display_id]*pur) ### cellular prevalence by multiplying with purity
    
    write.table(T2,col.names = c("cluster","n_ssms","proportion"),row.names=F,quote=F,sep="\t",file=f2)
    
    f3 = sprintf("./%s/%s_mutation_assignments.txt",dir_path,sample_id)
    chr_name = as.character(data[,CHR_NAME_INDX])
    chr_loc = data[,CHR_POS_INDX]
    T3 = cbind(chr_name,chr_loc,sorted_mut_assign_vec)
    write.table(T3,col.names=c("chr","pos","cluster"),row.names=F,quote=F,sep="\t",file=f3)
    ################
    
    f4 = sprintf("./%s/%s_number_of_clusters.txt",dir_path,sample_id)
    T4 = cbind(sample_id,new_number_of_cluster)
    write.table(T4,col.names=F,row.names=F,quote=F,sep="\t",file=f4)
    
    nSNV = length(ccf)
    
    cnt_clonal = new_size[which.max(new_mean_ccf)]
    cat("cnt_clonal = ",cnt_clonal,"\n")
    cnt_subclonal = nSNV - cnt_clonal
    clonal_snv_prop = cnt_clonal/nSNV
    
    f5 = sprintf("./%s/%s_clonal_fraction.txt",dir_path,sample_id)
    T5 = cbind(sample_id,cnt_clonal,cnt_subclonal,clonal_snv_prop)
    write.table(T5,col.names=F,row.names=F,quote=F,sep="\t",file=f5)
    
}

generate_result_based_on_mclust <- function(data,CCF_Info,clonal_flag,pur,dir_path,config_param){
    
    
    CHR_NAME_INDX = config_param$CHR_NAME_INDX
    CHR_POS_INDX = config_param$CHR_POS_INDX
    
    TOTAL_RD_INDX = config_param$TOTAL_RD_INDX
    VARIANT_RD_INDX = config_param$VARIANT_RD_INDX
    
    SAMPLE_T_CN_INDX = config_param$SAMPLE_T_CN_INDX
    
    MAJOR1_CN_INDX = config_param$MAJOR1_CN_INDX
    MINOR1_CN_INDX = config_param$MINOR1_CN_INDX
    FRAC1_INDX = config_param$FRAC1_INDX
    
    MAJOR2_CN_INDX = config_param$MAJOR2_CN_INDX
    MINOR2_CN_INDX = config_param$MINOR2_CN_INDX
    FRAC2_INDX = config_param$FRAC2_INDX
    KMEAN_CENTER_CUTOFF = config_param$KMEAN_CENTER_CUTOFF
    
    #nc_to_run_with = nc_mclust[i,2]+1 ### add one (clonal population)
    
    new_number_of_cluster = length(CCF_Info$mean_ccf)-1
    
    mean_ccf = CCF_Info$mean_ccf
    old_label = CCF_Info$old_label
    
    mut_assign_vec = CCF_Info$mut_assign_vec
    old_size = CCF_Info$old_size
    
    I0=which(mean_ccf > KMEAN_CENTER_CUTOFF)
    #I0=which(KK$centers >= KMEAN_CENTER_CUTOFF) ### or 0.95 ?
    #I0 = which(mean_ccf >= 1-((1-KMEAN_CENTER_CUTOFF)/pur))
	cat("mean_ccf from mclust and merge = ",mean_ccf,"\n")
    cat("I0 = ",I0,"\n")
    
    mean_ccf[I0] = 1.0
    I1=old_label[I0]
    cat("I1 = ",I1,"\n")
    cat("old_label = ",old_label,"\n")
    temp_label = setdiff(old_label,I1)
    cat("temp_label = ",temp_label,"\n")
    
    temp_size = old_size[temp_label]
    
    cat("temp_size",temp_size,"\n")
    cat("old_size",old_size,"\n")
    new_mean_ccf = mean_ccf
    new_mut_assign_vec = mut_assign_vec
    new_size = old_size
    new_label = old_label ###
    
    L1 = length(I1)
    if(L1>1){
        
        min_id1 = min(I1)
        size_min_id1 = sum(old_size[I1])
        
        snv_indx_to_change = which(mut_assign_vec %in% I1 )
        new_mut_assign_vec[snv_indx_to_change] = min_id1
        
        new_label = c(temp_label,min_id1)    
        
        new_size = c(temp_size,size_min_id1)
        
        new_mean_ccf = mean_ccf[new_label]
        
        new_number_of_cluster = length(new_label) - 1
        
        
        cat("new label = ",new_label,"\n")
        cat("new_size = ",new_size,"\n")
        cat("new_number_of_cluster = ",new_number_of_cluster,"\n")
        
        #sorted_label = sort(new_mean_ccf,index.return=T)$ix
        sorted_label = rank(new_mean_ccf)
        sorted_mut_assign_vec = new_mut_assign_vec
        
        for(j in 1:length(new_label)){
            mut_id = which(new_mut_assign_vec == new_label[j])
            sorted_mut_assign_vec[mut_id] = sorted_label[j] 
        }
        cat("new_mean_ccf = ",new_mean_ccf,"\n")
        cat("sorted_label = ",sorted_label,"\n")
        
        display_id = sort(new_mean_ccf,index.return=T)$ix
        cat("display_id = ",display_id,"\n")
        
    }else{
        sorted_label = rank(new_mean_ccf)
        sorted_mut_assign_vec = new_mut_assign_vec
        for(j in 1:length(new_label)){
            mut_id = which(new_mut_assign_vec == new_label[j])
            sorted_mut_assign_vec[mut_id] = sorted_label[j] 
        }
        display_id = sort(new_mean_ccf,index.return=T)$ix
        cat("display_id = ",display_id,"\n")
        cat("else part\n")
    }
    
    f2 = sprintf("./%s/%s_subclonal_structure.txt",dir_path,sample_id)
    T2 = cbind(sorted_label[display_id],new_size[display_id],new_mean_ccf[display_id]*pur) ### cellular prevalence by multiplying with purity
    
    write.table(T2,col.names = c("cluster","n_ssms","proportion"),row.names=F,quote=F,sep="\t",file=f2)
    
    f3 = sprintf("./%s/%s_mutation_assignments.txt",dir_path,sample_id)
    chr_name = as.character(data[,CHR_NAME_INDX])
    chr_loc = data[,CHR_POS_INDX]
    T3 = cbind(chr_name,chr_loc,sorted_mut_assign_vec)
    write.table(T3,col.names=c("chr","pos","cluster"),row.names=F,quote=F,sep="\t",file=f3)
    ################
    
    f4 = sprintf("./%s/%s_number_of_clusters.txt",dir_path,sample_id)
    T4 = cbind(sample_id,new_number_of_cluster)
    write.table(T4,col.names=F,row.names=F,quote=F,sep="\t",file=f4)
    
    nSNV = length(mut_assign_vec)
    
    cnt_clonal = new_size[which.max(new_mean_ccf)]
    cat("cnt_clonal = ",cnt_clonal,"\n")
    cnt_subclonal = nSNV - cnt_clonal
    clonal_snv_prop = cnt_clonal/nSNV
    
    f5 = sprintf("./%s/%s_clonal_fraction.txt",dir_path,sample_id)
    T5 = cbind(sample_id,cnt_clonal,cnt_subclonal,clonal_snv_prop)
    write.table(T5,col.names=F,row.names=F,quote=F,sep="\t",file=f5)
}

###############################

require(mclust)
require(fpc)
args <- commandArgs(TRUE)
set.seed(17345)

###############################

config_param = NULL

config_param$CHR_NAME_INDX = 1
config_param$CHR_POS_INDX = 2

config_param$TOTAL_RD_INDX = 3
config_param$VARIANT_RD_INDX = 4

config_param$SAMPLE_T_CN_INDX = 5

config_param$MAJOR1_CN_INDX = 6
config_param$MINOR1_CN_INDX = 7
config_param$FRAC1_INDX = 8
config_param$MAJOR2_CN_INDX = 9
config_param$MINOR2_CN_INDX = 10
config_param$FRAC2_INDX = 11

config_param$MAX_CLUSTER = 6 ### 7 for PCAWG
config_param$MIN_CLUST_DIST = 0.1

config_param$MAX_CN_CUTOFF = 6.0  ### for multiplicty calculation
config_param$PLOIDY_MAX_CN_CUTOFF = 15.0  ### for ploidy calculation previously 20
config_param$KMEAN_CENTER_CUTOFF = 0.9 #0.8 0.95

config_param$SAMPLE_NAME_INDX = 1
config_param$PURITY_INDX = 2
config_param$PLOIDY_INDX = 3
##############################
#sample_id = args[1]
sample_id = "sample"
input.file = args[1] 
pur_ploidy_file_name = args[2]
#R_code_src_dir = args[3]
#rda_output_dir = args[4]
#out_fig_folder = args[5]
out_txt_folder = "./"
##############################

#consensus_purity_ploidy_table = read.table("~/ICGC_PROJECT/consensus.20170119.purity.ploidy.annotated.txt",header=T)
#consensus_purity_ploidy_table = read.table("./pp_table.txt",header=T) ### true purity
#consensus_purity_ploidy_table = read.table("./pp_table_est.txt",header=T)  ### est purity
#consensus_purity_ploidy_table = read.table("./consensus.20170119.purity.ploidy.annotated.txt",header=T) ### true purity

SAMPLE_NAME_INDX = config_param$SAMPLE_NAME_INDX
PURITY_INDX = config_param$PURITY_INDX
PLOIDY_INDX = config_param$PLOIDY_INDX
#id = which(consensus_purity_ploidy_table[,SAMPLE_NAME_INDX] == sample_id)
#consensus_purity = consensus_purity_ploidy_table[id,PURITY_INDX]
#consensus_ploidy = consensus_purity_ploidy_table[id,PLOIDY_INDX]

##############################
data = data.frame(read.table(input.file))

#### Jan 28 update to read from individual purity file
dir_input_file = dirname(input.file)
#pur_ploidy_file_name = sprintf("%s/%s.pur.ploidy.txt",dir_input_file,sample_id)
pur_ploidy_val_str = read.table(pur_ploidy_file_name,header=T)
pur_val_from_file = as.numeric(unlist(pur_ploidy_val_str)[PURITY_INDX])
ploidy_val_from_file = as.numeric(unlist(pur_ploidy_val_str)[PLOIDY_INDX])
##########

#pur_val_from_file = consensus_purity
#ploidy_val_from_file = consensus_ploidy

cat("purity and ploidy = ",pur_val_from_file,ploidy_val_from_file,"\n")
#res_clonal = clonality_test_for_sample(data,consensus_purity,config_param)
#res_clonal = clonality_test_for_sample(data,pur_val_from_file,config_param)
#clonal_flag = res_clonal$f

clonal_flag = 0
cat("clonal_flag = ",clonal_flag,"\n")

#res_ccf = find_ccf_for_each_SNV(data,consensus_purity,consensus_ploidy, clonal_flag,out_txt_folder,config_param)
res_ccf = find_ccf_for_each_SNV(data,pur_val_from_file,ploidy_val_from_file, clonal_flag,out_txt_folder,config_param)

#nc_from_mclust = find_nclust_ccf_of_SNV(res_ccf$ccf,clonal_flag,config_param)
CCF_Info = NULL
CCF_Info = find_nclust_ccf_of_SNV(res_ccf$ccf,clonal_flag,config_param)
#cat("nc_from_mclust = ",nc_from_mclust,"\n")

#nc_to_run_with = nc_from_mclust+1
#find_cluster_center_by_kmean(data,res_ccf$ccf,clonal_flag,consensus_purity,nc_to_run_with,out_txt_folder,config_param)
#find_cluster_center_by_kmean(data,res_ccf$ccf,clonal_flag,pur_val_from_file,nc_to_run_with,out_txt_folder,config_param)

generate_result_based_on_mclust(data,CCF_Info,clonal_flag,pur_val_from_file,out_txt_folder,config_param)

##############################
cat("\n analysis finished and generated files for sample ",sample_id,"!!!\n")
##############################

