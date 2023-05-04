# get and assign are 2 important functions in R
ls_lst <- ls() # obtain the objects in the workspace
char_2_num_df <- function() {
  for(dataset in ls_lst){
    df_name <- dataset
    print(dataset)
    if(is.data.frame(get(dataset))){ # here isdataframe checks if the object is a dataframe
      dataset <- get(dataset) # get will mold the character string to object
      dataset <- as.data.frame(apply(dataset,2,as.numeric)) 
      assign(x = df_name, value = dataset,envir = globalenv())# assign is needed to assign the result to a object
    }else{
      print("not a dataframe")}
  } 
  
}




