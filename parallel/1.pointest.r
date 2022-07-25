library(data.table)
library(parallel)
library(MASS)

options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

fn_y1 = args[1]
fn_y2 = args[2]
grm_path = args[3]
grm_id_path = args[4]
dout = args[5]


###
matrixV = function(GRM,v1_n,N){
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
    i_v <- c(1,2,3,4)
    j_v <- c(1,3,2,4)
    
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
    v1 = c(1,1,1,1,2,2,2,3,3,12,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,11,11,14)
    v2 = c(1,2,3,12,2,3,12,3,12,12,4,5,6,7,11,14,5,6,7,11,14,6,7,11,14,7,11,14,11,14,14)
    A = matrix(0,n_para,n_para)
    for (row in 1:length(v1)){
        V1 = matrixV(GRM,v1[row],N)
        V2 = matrixV(GRM,v2[row],N)
        A[v1[row],v2[row]] = traceAB(V1,V2,n_cores) 
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
    A = matrix(0,n_para,1)
    rowlist= c(1:n_para)
    B = c()
    for(row in 1:n_para){
        
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


npara = 14
N_grm = 34162
M_grm = 4736711
numCores <- detectCores()
n_cores = floor(numCores/2)


y1 = as.data.frame(fread(fn_y1)) ## FID, IID, y
y2 = as.data.frame(fread(fn_y2))

## read GRMs 
grm_all = as.matrix(fread(grm_path,fill=TRUE,sep='\t')) ## takes a while


## read ind lists
grm_id = as.data.frame(fread(grm_id_path))
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
y1_std = (y_all[,1]-mean(y_all[,1]))/sd(y_all[,1])
y2_std = (y_all[,2]-mean(y_all[,2]))/sd(y_all[,2])
A = matrixA(npara,grm_all,N,n_cores) 
B = matrixB(grm_all,y1_std,y2_std,npara,N,n_cores)
est_all = solve(A,B)

setwd(dout)
write.table(A,paste0(pheno_y1,'_',pheno_y2,'_A.txt'),quote=F,row.names=F,col.names=F,sep='\t')
write.table(B,paste0(pheno_y1,'_',pheno_y2,'_B.txt'),quote=F,row.names=F,col.names=F,sep='\t')
write.table(est_all,paste0(pheno_y1,'_',pheno_y2,'_est_all.txt'),quote=F,row.names=F,col.names=F,sep='\t')
 
