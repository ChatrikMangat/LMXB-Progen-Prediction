import numpy as np
import pickle
import time
import os.path
import sys

# Progenitor Grid
m_arrl = np.arange(4.0,7.1,0.1)
m_arrs = np.arange(0.95,4.05,0.05)
p_arrl = np.arange(1.65,4.05,0.05)
p_arrs = np.arange(-0.60,1.66,0.02)

# Black hole masses to segregate databases 
bh_masses = [7.0,10.0]

# Grid separation w.r.t density since simulations are stored in four folders
confs = {'lmlp':[m_arrl,p_arrl],'smlp':[m_arrs,p_arrl],'lmsp':[m_arrl,p_arrs],'smsp':[m_arrs,p_arrs]}

# m1 -> Donor Mass {Msol}
# m2 -> Accretor Mass {Msol}
# mt -> log10(MT Rate) {Msol/yr}
# p -> Orbital Period {days}
# teff -> log10(Donor Effective Temp) {K}
labels = ['m1','m2','mt','p','teff']

errors = []


def search(array, element):
"""! Binary search algorithm used for Donor Mass.

@param array  The array to be searched.
@param element The element to be searhed for.
@return Index of element in array.
 
"""
    mid = 0
    start = 0
    end = len(array)
    step = 0

    while (start <= end):
        step = step+1
        mid = (start + end) // 2
        
        if mid == len(array):
            return len(array)-1

        if element == array[mid]:
            return mid

        if element > array[mid]:
            end = mid - 1
        else:
            start = mid + 1
    return start


def get_query(qpath):
"""! Function to retrieve query information from .txt file and return a dictionary for the query.

@param qpath Path to query_name.txt file.
@return Query in dictionary form or -1 if the input format is incorrect.

"""    
    try:
        q = open(qpath,"r")
        bhns = (q.readline()).split(None, 1)[0]
        q.close()

        qtemp = np.loadtxt(qpath,delimiter=",",skiprows=1)
    except Exception as e:
        errors.append("Error in query input, please recheck format\n"+str(e))
        return -1

    #print(qtemp.shape)
    query = {'bhns': int(bhns)}

    for i in range(qtemp.shape[0]):
        query[labels[i]] = qtemp[i]
    
    # Various input checks
    if (query['m1'][0] < 0.0 or query['m1'][0] > query['m1'][1]):
        errors.append("Wrong Donor Mass range")
        return -1
    if (query['m2'][0] < 0.0 or query['m2'][0] > query['m2'][1]):
        errors.append("Wrong Accretor Mass range")
        return -1
    if (query['mt'][0] < -100.0 or query['mt'][0] > query['mt'][1]):
        errors.append("Wrong MT Rate range")
        return -1
    if (query['p'][0] < 0.0 or query['p'][0] > query['p'][1]):
        errors.append("Wrong Orbital Period range")
        return -1
    if query['teff'][0] > query['teff'][1]:
        errors.append("Wrong Donor Teff range")
        return -1

    return query


def match_props(data,query):
"""! Function to check whether a simulated system (data) matches all query properties. 

@param data Array of properties taken from simulation data.
@param query Dictionary form of query.
@return 1 if simulated system matches query and 0 if it does not.

"""
    try:
        if 'm1' in query:
            if (data[0] < query['m1'][0] or data[0] > query['m1'][1]):
                #print("No m1")
                return 0
        if 'm2' in query:
            if (data[1] < query['m2'][0] or data[1] > query['m2'][1]):
                #print("No m2")
                return 0
        if 'mt' in query:
            if (data[2] < query['mt'][0] or data[2] > query['mt'][1]):
                #print("No mt")
                return 0
        if 'p' in query:
            if (data[3] < query['p'][0] or data[3] > query['p'][1]):
                #print("No p")
                return 0
        if 'teff' in query:
            if (data[4] < query['teff'][0] or data[4] > query['teff'][1]):
                #print("No teff")
                return 0
    except Exception as e:
        errors.append("Error in matching properties")
        return 0

    return 1


def find_mt_start(mt_arr):
"""! Function to find index where system starts MT.

@param mt_arr Array containing MT values during simulation.
@return Index of start of MT.

"""
    idx_start = 0

    if len(mt_arr) <= 5:
        return idx_start

    for i in range(len(mt_arr)-3):
        if (mt_arr[i] >= -15 and mt_arr[i+1] >= -15 and mt_arr[i+2] >= -15 and mt_arr[i+3] >= -15):
            idx_start = i
            break
    
    return idx_start


def get_progens(query):
"""! Function to look through data files to find progenitors for the query.

@param query Dictionary form of query.
@return Array containing progenitor properties.

"""    
    progens = []
    
    data_paths = []
    # Select relevant databases to search through
    if query['bhns'] == 0:
        data_paths.append('/home/chatriks/search_test/ns_data/')
    else:
        for x in bh_masses:
            if(x > query['m2'][1]):
                continue
            if(x < query['m2'][0]-3.5):
                continue
            data_paths.append('/home/chatriks/search_test/runs'+str(int(x))+'_data/')
    print(data_paths)

    # Loop through databases
    for dpath in data_paths:
        # Loop through four folders
        for c in confs.keys():
            # Loops for initial mass and period values 
            for i in confs[c][0]:
                for j in confs[c][1]:
                    
                    # Ignore simulations where initial donor mass is less than the lower limit of donor mass in query
                    if i >= query['m1'][0]:
                        try:
                            fpath = dpath+c+'/m_'+f'{i:4.2f}'+'_p_'+f'{j:4.2f}'+'.data'
                            if os.path.exists(fpath):
                                infile = open(fpath,'rb')
                                vals = pickle.load(infile)
                                infile.close()
                                
                                #print(fpath)
                                #print(vals.shape)
                                
                                # Ignore simulations where final donor mass is greater than upper limit of donor mass in query. Similarly, ignore simulations where initial accretor mass is greater than upper limit of accretor mass in query
                                if (vals[0][-1] <= query['m1'][1] and vals[1][0] <= query['m2'][1]):
                                    # Flags to indicate if system is a progenitor
                                    flag = 0
                                    pflag = 0
                                    # Indices in array to define window of entries where donor mass satisfies query donor mass limits
                                    start_m = 0
                                    end_m = 0
                                    
                                    # If upper limit of donor mass in query is greater than initial donor mass, start window at top of file
                                    if vals[0][0] <= query['m1'][1]:
                                        start_m = 0
                                    else:
                                    # Binary search for upper donor mass limit
                                        start_m = search(vals[0],query['m1'][0])-1
                                    
                                    # If lower limit of donor mass in query is less than final donor mass end window at bottom of file
                                    if vals[0][-1] >= query['m1'][0]:
                                        end_m = vals.shape[1]-1
                                    else:
                                    # Binary search for lower donor mass limit
                                        end_m = search(vals[0],query['m1'][0])
                                    
                                    # If donor mass window is found, go through simulation entries one by one to match all system properties to query. Flag all simulations with a matching system
                                    if (start_m <= end_m):
                                        #print("Mass matched ",start_m,end_m,end_m-start_m)
                                        idx_low = -1
                                        idx_high = -1
                                        obs_time = 0.0
                                        mt_start = 0

                                        for k in range(start_m,end_m+1):
                                            if (k < 0 or k >= vals.shape[1]):
                                                break
                                            flag = match_props([vals[0][k],vals[1][k],vals[2][k],vals[3][k],vals[4][k]],query)
                                            
                                            # If a match is found, find the time spent as observed system, each traversal across the window is added to total observed time incrementally
                                            if (flag == 1 and idx_low == -1):
                                                idx_low = k

                                            if ((flag == 0 and idx_low != -1 and idx_high == -1) or (flag == 1 and k == end_m and idx_high == -1)):
                                                pflag = 1
                                                idx_high = k
                                                obs_time = obs_time + vals[5][idx_high]-vals[5][idx_low]
                                                idx_low = -1
                                                idx_high = -1

                                        # If simulated system spends time as the observed system, find start of MT information and add to list of progenitors
                                        if pflag == 1:
                                            mt_start = find_mt_start(vals[2])
                                            tot_time = vals[5][-1]-vals[5][0]
        
                                            progens.append([i,j,vals[1][0],obs_time,tot_time,vals[0][mt_start],np.log10(vals[3][mt_start])])
                                            print("Found Progenitor -> "+fpath)
                                
                            else:
                                print("No path found: "+fpath)

                        except Exception as e:
                            errors.append("Error occurred in "+fpath+"\n"+str(e)+"\n")
    return progens


def output_progens(rpath,progens):
"""! Function to store progenitor information for the given query in a .txt file.

@param rpath Path to progenitor .txt file.
@param Array containing progenitor properties.

"""
    try:
        outfile = open(rpath,'w')

        outfile.write("Number of Progenitors Found: "+str(len(progens))+"\n")
        outfile.write("Progenitors:\n") 
        for i in progens:
            s = f'{i[0]:.2f}'+' '+f'{i[1]:.2f}'+' '+f'{i[2]:.2f}'+' '+f'{i[3]:.4f}'+' '+f'{i[4]:.4f}'+' '+f'{i[5]:.4f}'+' '+f'{i[6]:.4f}'+'\n'
            outfile.write(s)

        outfile.write("\nErrors:\n")
        if (len(errinstall doxygen ubuntuors) > 0):
            for i in errors:
                outfile.write(i+'\n')
        else:
            outfile.write("None")

        outfile.close()
        print("Progenitor properties stored in: "+rpath)
    
    except Exception as e:
        errors.append("Error occurred in writing output:\n"+str(e))
        error_exit(rpath)


def error_exit(rpath):
"""! Function to output fatal errors

@param Path to progenitor .txt file.

"""
    
    outfile = open(rpath,'w')
    outfile.write("Errors:\n")
    for i in errors:
        outfile.write(i+'\n')

    print("Error in query input: Check progens file for details")
    exit()


# Driver Code
print("Start Search")
start_time = time.time()

# Query .txt filename given as command line argument, change qpath according to where file is stored
qname = sys.argv[1]
qpath = ""+qname

# Progenitor list stored in corresponding .txt file, change rpath according to preferences
rpath = "progens_"+qname

query = get_query(qpath)
if (query == -1):
    error_exit(rpath)

print("Query input taken:")
print(query)

progens = get_progens(query)
output_progens(rpath,progens)

print("Search Complete, time taken: "+str(time.time() - start_time)+" seconds")
