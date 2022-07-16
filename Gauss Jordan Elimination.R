
## GAUSS-JORDAN ELIMINATION
## Roll: bs2016, Debabrata Sarkar


#M = function(){
#  mat[i,] <- (mat[i,]/mat[i,j]);
#}

#S = function(l){
  
#  mat[l,] <- mat[l,] - (mat[l,]*mat[l,j]);
#}

rref = function(mat){

m = nrow(mat);
n = ncol(mat);
rhead = rownames(mat);
#chead = colnames(mat);
A = mat

i=1;
j=1;
  for(t in 1:m){ # we run the process for every row
  
    # the first element is the first pivot
    # we know the pivot cannot be 0, so we adjust in the following way
    if( mat[i,j] == 0){ 
      indicator = 0;
      
      for(k in i:m){
        
        if( mat[k,i] != 0){ # next non-zero pivot
          
          # if such a pivot is present in the k-th row for the 1st time
          # swap rows k and i
          temp = mat[k,];
          mat[k,] <- mat[i,];
          mat[i,] <- temp;
          
          # also swap the row heads
          tmp = rhead[k];
          rhead[k] = rhead[i];
          rhead[i] = tmp;
          
          rownames(mat) <- rhead;
          
          break;
        }
        
        if(k == m){ # we have reached the last row
          indicator = 1; # all non zero pivotal rows are at the top now
        }
      }
      # at this stage, no more non zero elements are left under the current element
      # we move the pivot to the right
      #if(j <= n){
        if(indicator == 1){
          
          j <- j+1;
        }
      if(j > n){
        break;
      }
      
      if( mat[i,j] != 0){ # check if current pivot is nonzero
        
        # carry out the M step
        # M step: divide the pivotal row by the pivot (undefined if pivot=0)
        mat[i,] <- (mat[i,]/mat[i,j]);
      }
      
      for(l in 1:m){
        if(l != i){
          
          # carry out the S step
          # S step: sweep the pivotal column (all elements below pivot = 0)
          mat[l,] <- mat[l,] - (mat[l,]*mat[l,j]);
        }
        #else{
        #  next;
        #}
      }
    }
    
    i <- i+1;
    j <- j+1;
    
  }
  print("The row reduced echelon form (RREF) of the matrix: ");
  print(mat);
  
  # Finding the basis of row space
  # The basis vectors of the row space are the non-zero rows in the RREF of the matrix
  
  for(i in 1:m){
    #for(j in 1:n){
      #if( sum(mat[i,]^2) != 0){
      
    if( sum(mat[i,]*mat[i,]) != 0){ # row is non-zero, if norm of row vector is non zero 
        print("The basis vectors for the row space:");
        print(mat[i,]);
      }
    }
  
  # Finding the basis of column space
  # The basis vectors are those columns of the original matrix corresponding to
  # the columns of the RREF matrix containing the leading 1's of the non-zero rows
  
  #for(j in 1:n){
    for(i in 1:m){
      for(j in 1:n){ # run through each element of a row to find the first 1
        
      if( mat[i,j] == 1){
        print("The basis vectors for the column space:");
        print(A[,j]);
      break;
      }
    }
  }

  B = matrix(0, nrow = n, ncol = n);
  p = 1;
  
  for(i in 1 : m){   
    p = n;
    
    for(j in 1 : n){
      if(mat[i, j] == 1){ 
        
        p = j; 
        break;
        }
      
    }
    
    for(j in (p+1): n){
      
      if(j > n){
        break;
      }
      
      B[p, j] = mat[i, j];
      
      if(mat[i, j] != 0)
        B[j, j] = -1;
    }
  }
  
  
  print("Nullspace basis : ");
  
  indicator1 = 0;
  
  for(i in 1 : n){
    if(sum(B[,i]*B[,i]) > 0){
      
      print(B[, i]);
      
      indicator1 = 1;
      }
  }
  if(indicator1 == 0)
    print( rep(0, n));
}
