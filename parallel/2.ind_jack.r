library(data.table)
library(parallel)
library(MASS)

options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

pheno_y1 = args[1]
pheno_y2 = args[2]
k_ind = as.integer(args[3])


###
matrixV = function(GRM,v1_n,N){
    # '''
    # generate V based on GRM and individual blocks & snp blocks
    # Vs are symmetric matrices
    # GRM: 2N*2N
    # return: dictionary (or other data structure) that contains 14 Vs, each of dimension 2N*2N (N = the number of siblings)
    # '''
    V = list()
    V[[1]] = 0; V[[2]] = 0; V[[3]] = 0; V[[4]] = 0
    if (v1_n==1) {
        ## GRM top left block
        index = 1:N
        V[[1]] = GRM[index,index]-diag(1,N)
    }else if (v1_n==2){
        ## GRM bottem right block
        index = (N+1):(2*N)
        V[[1]] = 4*GRM[index,index]-diag(1,N)*2
    }else if (v1_n==3){
        M = matrix(c(0,1,1,0),2,2)
        V[[1]] = kronecker(diag(1,N/2), M)
    }else if ( v1_n==12){
        index1 = 1:N
        index2 = (N+1):(2*N)
        V[[1]] = 2*(GRM[index1,index2]+ GRM[index2,index1])-diag(1,N)*2
    }else if (v1_n==4){
        index = 1:N
        V[[2]] = GRM[index,index]  
        V[[3]] = V[[2]]
    }else if ( v1_n==5){
        index = (N+1):(2*N)
        V[[2]] = GRM[index,index]*4
        V[[3]] = V[[2]]
    }else if ( v1_n==6){
        V[[2]] = diag(1,N)
        V[[3]] = V[[2]]
    }else if (v1_n==7){
        M = matrix(c(0,1,1,0),2,2)
        V[[2]] = kronecker(diag(1,N/2), M)
        V[[3]] = V[[2]]
    }else if (v1_n==11){
        index1 = 1:N
        index2 = (N+1):(2*N)
        V[[2]] = 2*GRM[index1,index2]
        V[[3]] = 2*GRM[index2,index1]
    }else if (v1_n==14){
        index1 = 1:N
        index2 = (N+1):(2*N)
        V[[3]] = 2*GRM[index1,index2]
        V[[2]] = 2*GRM[index2,index1]
    }else if (v1_n==8){
        index = 1:N
        V[[4]] = GRM[index,index]-diag(1,N)
    }else if (v1_n==9){
        index = (N+1):(2*N)
        V[[4]] = 4*GRM[index,index]-diag(1,N)*2
    }else if (v1_n==10){
        M = matrix(c(0,1,1,0),2,2)
        V[[4]] = kronecker(diag(1,N/2), M)
    }else if (v1_n==13){
        index1 = 1:N
        index2 = (N+1):(2*N)
        V[[4]] = 2*(GRM[index1,index2] + GRM[index2,index1])-diag(1,N)*2
    }
    return(V)
}



traceAB = function(V1,V2,n_cores=1){
    # '''
    # A and B are both blockwise matrices with four blocks like [[A1,A2],[A3,A4]]
    # trace(AB) = trace(A1B1+A2B3+A3B2+A4B4)
    # V1,V2 contains GRM index for each block submatrix
    # '''
  
    i_v <- c(1,2,3,4)
    j_v <- c(1,3,2,4)
    
    # result = mcmapply(i=i_v,j=j_v,function(i,j){
    result = c()
    for(k in 1:4){
      i=i_v[k]
      j=j_v[k]
      if((!is.null(dim(V1[[i]])))&(!is.null(dim(V2[[j]])))){
        result = c(result,sum(V1[[i]] * V2[[j]]))
      }else{
        result = c(result,0)
      }
    }

    tr = sum(result)
    return(tr)
}



matrixA = function(n_para,GRM,N,n_cores){
    # about 21 min for one matrix A
    # '''
    # matrixA is symmetric matrix, so just need to calculate lower triangle
    # grm : grm matrix
    # n_para : number of parameters
    # '''
    v1 = c(1,1,1,1,2,2,2,3,3,12,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,11,11,14)
    v2 = c(1,2,3,12,2,3,12,3,12,12,4,5,6,7,11,14,5,6,7,11,14,6,7,11,14,7,11,14,11,14,14)
    A = matrix(0,n_para,n_para)
    for (row in 1:length(v1)){
        a=Sys.time()
        V1 = matrixV(GRM,v1[row],N)
        V2 = matrixV(GRM,v2[row],N)
        A[v1[row],v2[row]] = traceAB(V1,V2,n_cores) 
        print(Sys.time()-a)
        #    user  system elapsed 
        # 28.914  10.293  40.120 
    }
    A[8,8] = A[1,1]
    A[9,9] = A[2,2]
    A[10,10] = A[3,3]
    A[8,9] = A[1,2]
    A[8,10] = A[1,3]
    A[9,10] = A[2,3]
    A[8,13] = A[1,12]
    A[9,13] = A[2,12]
    A[10,13] = A[3,12]
    A[13,13] = A[12,12]
    A = A+ t(A)
    diag(A) = diag(A)/2
    return(A)
}

matrixB = function(GRM,y1,y2,n_para,N,n_cores){
    # '''
    # matrix B = yVy-trace(V)
    # '''
    A = matrix(0,n_para,1)
    rowlist= c(1:n_para)
    # B = mcmapply(row=rowlist,function(row){
    B = c()
    for(row in 1:n_para){
        print(row)
      
      V = matrixV(GRM,row,N)
      term1 = 0; term2 = 0
      for(i in 1:4){
        if(is.null(dim(V[[i]]))){
            next
        }
        if(i==1){
            term1 = term1 + t(y1) %*% V[[1]] %*% y1
        } else if(i==2){
            term1 = term1 + t(y1) %*% V[[2]] %*% y2
        } else if(i==3){
            term1 = term1 + t(y2) %*% V[[3]] %*% y1
        } else if(i==4){
            term1 = term1 + t(y2) %*% V[[4]] %*% y2
        }
        term2 = term2 + sum(diag(V[[i]]))
      }
      B = c(B,term1 - term2)
      
    }
    
    return(B)
}


# calculate = function(K_ind,K_snp,N_grm,M_grm,fn_y1,fn_y2,dout){
    # """
    # N: number of siblings
    # M: number of SNPs
    # K_ind: number of individual blocks
    # K_snp: number of SNP blocks
    # """

    ## number of parameters = 14
    npara = 14
    N_grm = 34162
    M_grm = 4736711
    numCores <- detectCores()
    n_cores = floor(numCores/2)
    ## number of individual blocks
    K_ind = 200
    ## number of snp blocks
    # K_snp = 900

    ## read phenotypes
    # dpheno = '/z/Comp/lu_group/Members/jsong/Minorjobs/clarazou/SiblingGC/5000_mom/sibling_phenotype'
    # id_sub = '/z/Comp/lu_group/Members/jsong/Sibling/UKB/Traits/jackknife_snp/codes/indlist.txt'
    # id = '/z/Comp/lu_group/Members/jsong/Sibling/fid_2sib.txt'

    # fn = c('phenotype_Height.txt','phenotype_BMI.txt','phenotype_HouseholdIncome.txt','phenotype_EduYears.txt','phenotype_OverallHealth.txt')
    # n = c(34104,34077,29501,33886,34048)
    # h = read.table('phenotype_OverallHealth.txt')
    # summary(h)
    # sum(!is.na(h))
    # h=cbind(id,h)
    # h[h[,3]<0,3]=NA
    # write.table(h,paste0(dout,'OverallHealth.txt'),quote=F,row.names=F,col.names=F,sep='\t')


    fn_y1 = paste0('/z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/phenotypes/',pheno_y1,'.txt')
    fn_y2 = paste0('/z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/phenotypes/',pheno_y2,'.txt')
    y1 = as.data.frame(fread(fn_y1)) ## FID, IID, y
    y2 = as.data.frame(fread(fn_y2))

    ## read GRMs 
    grm_all = as.matrix(fread('/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/line1_rmed/chr_all_grm_zs.rel',fill=TRUE,sep='\t')) ## takes a while


    ## read ind lists
    ## individual order in GRM file
    # fam.f = as.data.frame(fread('/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/par_all.fam'))
    grm_id = as.data.frame(fread('/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/line1_rmed/chr_all_grm_zs.rel.id'))
    colnames(grm_id)[1]='FID'
    grm_sib_id = grm_id[1:N_grm,]
    ## cleaning phenotypes
    pheno_complete = complete.cases(cbind(y1[,3],y2[,3]))
    y_all = cbind(y1[pheno_complete,3],y2[pheno_complete,3])
    y_all_id = y1[pheno_complete,2] #IID

    grm_id_left = grm_sib_id[grm_sib_id$IID %in% y_all_id,]
    ta = table(grm_id_left$FID)
    count = attr(ta, "dimnames")[[1]][ta==2]
    grm_id_left_families = grm_id_left[grm_id_left$FID %in% count,]
    

    N = dim(grm_id_left_families)[1]
    blocks_IID = grm_id_left_families[order(grm_id_left_families$FID),'IID'] ## reorder
    index = match(blocks_IID,grm_id$IID)
    index2 = c(index,index+N_grm)
    grm_all = grm_all[index2,index2]
    grm_id = grm_id[index,]


    y_all = y_all[match(blocks_IID,y_all_id),]
    y_all_id = y_all_id[match(blocks_IID,y_all_id)]


    ## Point est for the whole dataset
    # y1_std = (y_all[,1]-mean(y_all[,1]))/sd(y_all[,1])
    # y2_std = (y_all[,2]-mean(y_all[,2]))/sd(y_all[,2])
    # A = matrixA(npara,grm_all,N,n_cores) 
    # B = matrixB(grm_all,y1_std,y2_std,npara,N,n_cores)
    # est_all = solve(A,B)
  
    ## Ind Jackknife
    # partition indlist into different blocks
    # note: every individual have two IDs: family ID, i.e. FID, and individual ID, i.e. IID. Each family has two siblings, so that means there will be two individuals that share the same FID. We want to split samples based on FID, so the siblings in the same family will be partitioned into the same block. And, we want to organize the dataset in the way that two siblings are continent in the dataset. So the rows of the data will be like: sib1_fam1, sib2_fam1, sib1_fam2, sib2_fam2, sib1_fam3, sib2_fam3,...
    set.seed(1)
    spliting = sample(K_ind,length(unique(grm_id_left_families$FID)),replace=T)
    # est_indblock = matrix(0,K_ind,npara)
    
    # K_ind_list= c(1:K_ind)
    # est_indblock = mcmapply(k_ind=K_ind_list,function(row){
      
      blocks_FID = unique(grm_id_left_families$FID)[spliting==k_ind]
      blocks_tmp = grm_id[grm_id$FID %in% blocks_FID,]
      blocks_IID = blocks_tmp[order(blocks_tmp$FID),'IID'] ## reorder
      N_indblock = length(blocks_IID)
      index = match(blocks_IID,grm_id$IID)
      index_grm = c(index,index+N)
      grm_indblock = grm_all[-index_grm, -index_grm]
      A_indblock = matrixA(npara,grm_indblock,N-N_indblock,n_cores)

      index_y = match(blocks_IID,y_all_id)
      y1_block = y_all[-index_y,1]
      y1_block_std = (y1_block-mean(y1_block))/sd(y1_block)
      y2_block = y_all[-index_y,2]
      y2_block_std = (y2_block-mean(y2_block))/sd(y2_block)
    #   y_indblock = c(y1[-index],y2[-index])
      B_indblock = matrixB(grm_indblock,y1_block_std,y2_block_std,npara,N-N_indblock,n_cores)
      est_indblock = solve(A_indblock,B_indblock)


    dout = paste0('/z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/2.indjack/',pheno_y1,'_',pheno_y2)
    if(!dir.exists(dout)){
        dir.create(file.path(dout))
    }
    setwd(dout)
    write.table(A_indblock,paste0('A_',k_ind,'.txt'),quote=F,row.names=F,col.names=F,sep='\t')
    write.table(B_indblock,paste0('B_',k_ind,'.txt'),quote=F,row.names=F,col.names=F,sep='\t')
    write.table(est_indblock,paste0('est_',k_ind,'.txt'),quote=F,row.names=F,col.names=F,sep='\t')
  
    # },mc.cores=n_cores)
    
    # var_indjack = sum((est_indblock-rowMeans(est_indblock))^2)*(K_ind-1)/K_ind


    # ## SNP Jackknife
    # fam.f = '/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/par_all.fam'
    # map.f = '/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/line1_rmed/block1.map'
    # dosage.f = '/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/line1_rmed/chr1.dosage.gz'
    # dout = '/z/Comp/lu_group/Members/jsong/Sibling/coding/dosage/all_ordered/line1_rmed'

    # est_snpblock = matrix(0,K_snp,npara)
    
    # K_snp_list= c(1:K_snp)
    # est_indblock = mcmapply(k_snp=K_snp_list,function(row){
      
    #   blocks_snp = blocks_snp_all[k_snp]
    #   map.f = paste0(dout,'/map/block',k_snp,'.map') ## pre-written
    #   ## call PLINK
    #   m_block = dim(read.table(map.f))[1]
    #   command = paste0('/z/Comp/lu_group/Software/plink/plink2_linux_x86_64_20190708/plink2 --fam ',fam.f,' --map ',map.f,' --import-dosage ',dosage.f,' noheader --make-rel square --out ',dout,'/blockgrm/grm_block',k_snp)
    #   system(command)
    #   ## read block GRM 
    #   grm_snpblock = as.matrix(fread(paste0(dout,'/blockgrm/grm_block',k_snp,'.rel'))) 
    #   grm_snpblock = (grm_all*M_grm-grm_snpblock*m_block)/(M_grm-m_block)
    #   system(paste0('rm ',dout,'/blockgrm/grm_block',k_snp,'.rel'))
    #   A_snpblock = matrixA(npara,grm_snpblock,N,n_cores)
    #   y = np.concatenate((y1,y2))
    #   B_snpblock = matrixB( grm_snpblock,y,npara,N,n_cores)
    #   solve(A_snpblock,B_snpblock)
      
    # },mc.cores=n_cores)

    # var_snpjack = sum((est_snpblock-rowMeans(est_snpblock))^2)*(K_snp-1)/K_snp

    ## final output
    # return (data.frame(est_all,var_indjack,var_snpjack))
# }


# calculate(K_ind,K_snp,N,M,fn_y1,fn_y2,dout)


# dout = '/z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/2.indjack'
# # dout = paste0('/z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/2.indjack/',pheno_y1,'_',pheno_y2)
# dsh = paste0(dout,'/sh/')
# if(!dir.exists(dsh)){
#         dir.create(file.path(dsh))
# }
# setwd(dsh)

# id = c('height','bmi','OverallHealth','EA','income')

# for(i in 1:4){
#     for(j in (i+1):5){
#         dsh_ij = paste0(dsh,id[i],'_',id[j],'/')
#         if(!dir.exists(dsh_ij)){
#                 dir.create(file.path(dsh_ij))
#         }
#         for(k in 1:200){
#             pipeline_d = c('#!/bin/bash\n',paste0('/s/bin/Rscript /z/Comp/lu_group/Members/jsong/Sibling/coding/UKB/2.ind_jack.r ',id[i],' ',id[j],' ',k))
#             write.table(pipeline_d,paste0(dsh_ij,k,".sh"),row.names=F, col.names=F, quote=F, sep="\n")
#         }

#         t = paste0(dsh_ij, "condor_submit_settings.sh")

#         write("Universe = vanilla", file = t)
#         # cat("Environment = PYSPARK_SUBMIT_ARGS=\"\"--driver-memory 8g pyspark-shell\"\"\n", file = t, append = T)
#         cat("plusone = $(Process) + 1\n", file = t, append = T)
#         cat("NewProcess = $INT(plusone)\n", file = t, append = T)
#         cat("executable = $(NewProcess).sh\n", file = t, append = T)
#         cat("Request_memory = 120G\n", file = t, append = T)
#         cat("Request_cpus = 1\n", file = t, append = T)
#         cat("Output = simple$(NewProcess).out\n", file = t, append = T)
#         cat("Log = simple$(NewProcess).log\n", file = t, append = T)
#         cat("error = simple$(NewProcess).error\n", file = t, append = T)
#         cat("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n", file = t, append = T)
#         cat("periodic_release = (NumJobStarts < 10) && ((CurrentTime - EnteredCurrentStatus) > 300)\n", file = t, append = T)
#         cat("Queue", 200, "\n", file = t, append = T)

#         # make .sh files executable:
#         system(paste0("chmod a+x ", dsh_ij, "*.sh"))
#         setwd(dsh_ij)
#         system(paste0('condor_submit ',"condor_submit_settings.sh"))
#     }
# }

