##################################################################################
import csv
import os
import pandas as pd
import numpy as np
import math as mt
import matplotlib.pyplot as plt

##################################################################################
#Input parameters

#interval of mean
interval = 50
#starting iteration of mean
start = 300
#criteria of convergence
CoIterations = 50
# tolerance percentage
tol = 0.5
# tolerance degree
deg = 5

##################################################################################

#Read .dat files of 2D and 6D pressure files to form final cv
#path = os.getcwd()

# read 2D .dat to a list of lists
dat_2DContent = [i.strip().split() for i in open('./postProcessing/TwoDPlane/0/surfaceFieldValue.dat').readlines()]

# write it as a new CSV file
with open("./surfaceFieldValue_2D.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(dat_2DContent)

# read 6D .dat to a list of lists
dat_6DContent = [i.strip().split() for i in open('./postProcessing/SixDPlane/0/surfaceFieldValue.dat').readlines()]

# write it as a new CSV file
with open("./surfaceFieldValue_6D.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(dat_6DContent)


# reading the csv file
colnames = ['Iterations', '2D pressure']
#read 2D .csv
df = pd.read_csv("surfaceFieldValue_2D.csv", names=colnames, comment='#')
#read 6D .csv
colnames = ['Iterations', '6D pressure']
df1 = pd.read_csv("surfaceFieldValue_6D.csv", names=colnames, comment='#')

#print(df.head())
#print(df1.head())

df["6D pressure"] = df1["6D pressure"]

#print(df.head())

##################################################################################

#extract the flow rate from the text file

with open('./massFlowRate_Outlet.log') as openfile:
    for line in openfile:
        s = line.split('=')
        for part in s:
            #print(part)
            if "    sum(Outlet) of phi" in part:
                flow_rate = float(s[1])                

##################################################################################

# Cv calculations

# volumetric flow rate
#flow_rate = 0.0037241842
# pressure difference between the 2D and 6D plane 
df['deltaP']= (df["2D pressure"] - df["6D pressure"]).astype(float)
# flow coefficent
df["cv"] = (flow_rate * 15850.32)/np.sqrt(df['deltaP']*998.98*0.000145038)

#meanNumber = df.loc[1:10, "Iterations"]
#print("the mean number is ", meanNumber)
# print(df.head())
#print(df.head())
#print(df.tail())
# writing the Cv.csv file
df.to_csv("Cv.csv", encoding='utf-8', index=False)

##################################################################################

# mean Cv
cv = pd.read_csv("Cv.csv")

rows, column = cv.shape
#print(rows)

#cv.set_index("Iterations", inplace=True)
#print(cv.head())

#meanNumber = df.loc[1:10, "cv"]
#print("the mean number is ", meanNumber)

# no. of mean divisions
divisions = (rows - start)/interval
#print(divisions)

meanList = np.arange(start, rows+interval, interval)

#meanList = [[300,349], [350,399]]
n = meanList.size
if(meanList[n-1]>rows):
    meanList[n-1]=rows

#print("\nThe mean list values are :", meanList)
#print("\nThe size of list size is :", n)

meanCv = []
for i in range(0, n-1):
    meanCv.append((df.loc[meanList[i]:(meanList[i+1] - 1), "cv"].mean()))
    #print(meanCv[i])

print("\nThe mean Cv values are :",meanCv)
m = len(meanCv)
#print("\nThe size of mean Cv values is :",m)

#meanCv = pd.DataFrame(meanCv, columns = ["Mean Cv"])
#print(meanCv.head())

##################################################################################

#percentage difference with respect to Cv
diffCv = []
for i in range(0, n-1):
    for j in range(meanList[i], meanList[i+1]):
        diffCv.append((((df["cv"][j]) - meanCv[i])/meanCv[i])*100)

meandiffCv = []
for i in range(0, n-1):
    for j in range(meanList[i], meanList[i+1]):
        meandiffCv.append((((df["cv"][j]) - meanCv[i])))

#print(diffCv)
#print(meandiffCv)

diffCv = pd.DataFrame(diffCv, columns = ["Percentage mean difference of Cv"])
diffCv.to_csv("diffCv.csv", encoding='utf-8', index=False)

meandiffCv = pd.DataFrame(meandiffCv, columns = ["Mean difference of Cv"])
meandiffCv.to_csv("meandiffCv.csv", encoding='utf-8', index=False)
#print(diffCv.tail())
#print(meandiffCv.tail())


##################################################################################

k = 0
l = 0
length = diffCv.size
#print("\nThe size of diffCv size is :",length)
for i in range(0,length):
    s = diffCv["Percentage mean difference of Cv"][i]
    if (mt.fabs(s) < tol):
        k = k+1 
        if(k > CoIterations):
            slope = (df["cv"][start + (i-1)] - df["cv"][(start + (i-1))- CoIterations])/CoIterations
            #slope = (df["cv"][start + (i-1)] - df["cv"][(start + (i-1))- 1])
            #diff_Slope = (diffCv["Percentage mean difference of Cv"][start + (i-1)] - diffCv["Percentage mean difference of Cv"][(start + (i-1))- CoIterations])/CoIterations
            #diff_Slope = (diffCv["Percentage mean difference of Cv"][start + (i-1)] - diffCv["Percentage mean difference of Cv"][(start + (i-1))- 1])
            #diff_degree_slope = mt.fabs(mt.atan(diff_Slope)*57.2958) 
            degree_slope = mt.fabs(mt.atan(slope)*57.2958)
            #degree_slope = mt.fabs(mt.atan(diff_Slope)*57.2958)
            if(degree_slope < deg):
                l= l+1
                if(l > CoIterations):
                    #diff_Slope = (diffCv["Percentage mean difference of Cv"][start + (i-1)] - diffCv["Percentage mean difference of Cv"][(start + (i-1))- CoIterations])/CoIterations
                    #diff_degree_slope = mt.fabs(mt.atan(diff_Slope)*57.2958) 

                    print("\nThe solution is converged at iterations\t", df["Iterations"][start + (i-1)])
                    print("with the  degree of slope of cv change in last\t"+str(CoIterations)+"\titerations is\t", str(degree_slope))
                    #print("with the  degree of slope of diffcv change in last\t"+str(CoIterations)+"\titerations is\t", str(diff_degree_slope))
                    print("The final value of converged Cv is\t", df["cv"][start + (i)])
                    print("The final Cv value repoted at " + str(rows)+ "\tis\t" + str(df["cv"][rows - 1]))
                    break
            else:
                l=0
    else:
        k=0
##################################################################################

# plotting of Cv curve and percentage difference of Cv

#Cv plot
df.plot(kind = 'line', x = 'Iterations', y= 'cv',color = 'blue')
plt.ylabel('Flow coefficient (Cv)')
#plt.show()
plt.savefig('Cv.png')

#Percentage mean difference of Cv plot
diffCv.plot(kind = 'line', y= 'Percentage mean difference of Cv',color = 'red')
plt.xlabel('Iterations')
plt.ylabel('Percentage mean difference of Cv')
#plt.show()
plt.savefig('diffCv.png')

##################################################################################
