# script to move reclustered images to separate folders
import pandas as pd
import os,sys
import re, shutil

csv_path="/home/anubratadas/Documents/GBM_BBB_LAB_GLASGOW/multiplex_IHC/reclustered images/sample_reclustered.csv"

core_folder= "/home/anubratadas/Documents/GBM_BBB_LAB_GLASGOW/multiplex_IHC/core_images/"
margin_folder= "/home/anubratadas/Documents/GBM_BBB_LAB_GLASGOW/multiplex_IHC/margin_images/" 
output_folder= "/home/anubratadas/Documents/GBM_BBB_LAB_GLASGOW/multiplex_IHC/reclustered images/"

# read the csv file
def read_csv(csv_path:str):
    
    input_csv = pd.read_csv(csv_path)
    return input_csv[['sample','newcluster']]
# loop
def extract_info():
    input_df = read_csv(csv_path)
    sample_cluster_list=[]
    for index,row in input_df.iterrows():
        sample = row['sample']
        cluster  = row['newcluster']
        sample_cluster = sample,cluster
        sample_cluster_list.append(sample_cluster)
    return sample_cluster_list    

def reframe_path():
    sample_cluster_list = extract_info()
    for sc in sample_cluster_list:
        sample = sc[0]
        cluster = sc[1]
        sample_file = re.search(r"^[uU].+_\d{5}_[0-9]A?",sample) # start with U then underscore then 5 numbers then underscore
        sample_file = sample_file.group()
        sample_file = sample_file + "_set1_512px"
        image_file = sample.split("_")[-1]+".jpg"
        # print(image_file)
        if cluster == "tumor":
            # print(f"{sample} is a tumor sample")  
            image_path = core_folder
            tf_path = os.path.join(image_path,sample_file,"images")   
            yield tf_path,cluster,image_file        
        elif cluster == "margin":
            # print(f"{sample} is a margin sample")  
            image_path =  margin_folder
            tf_path = os.path.join(image_path,sample_file,"images")
            yield tf_path,cluster,image_file
        elif cluster == "common":
            # print(f"{sample} is a common sample")    
            tf_path = sample
            yield tf_path,cluster,image_file
        else:
            print(f"not recognized cluster")    
       

def transfer_file():
   output = reframe_path()
   datafile = open(os.path.join(output_folder,"common_images.txt"),"a+")
   for tf_path in output:
       image_path = tf_path[0]
       cluster = tf_path[1]
       image_file = tf_path[2]
       if cluster == "common":
           datafile.write(image_path+"\n")
       elif cluster == "tumor":              
           shutil.copy(os.path.join(image_path,image_file),os.path.join(output_folder,"tumor"))
       elif cluster == "margin":
           shutil.copy(os.path.join(image_path,image_file),os.path.join(output_folder,"margin"))      
       else: 
           print("unable to copy")    
   datafile.close()
        
transfer_file()
