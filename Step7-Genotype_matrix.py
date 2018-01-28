"""

Generate genotype matrix

"""
import sys
import glob

GLB=glob.glob('*.csv')

def getVect(csv_file):
   d={}
   with open(csv_file, 'r') as inp:
     for line in inp:
       A=line.strip().split(',')
       if len(A) == 2:
          d[A[0]]=A[1]
   return d   


megaS=set() 
tableD={}
ACC=[]
for csv_file  in GLB:
   acc=csv_file[:-4] 
   D=getVect(csv_file)
   megaS.update(D.keys())
   tableD[acc]=D
   ACC.append(acc)
   
#----------------------
if __name__ == "__main__":
  # Write both the fasta formatted and the tab separated file
  with open('Table.fas','w') as fas, open('Table.tsv','w') as out:
    out.write('\t'.join(['marker'] + ACC  ) + '\n')  # headerline  
    for u in megaS:
      pipe=[u]
      for acc in ACC:
         try:
          nt=tableD[acc][u]
         except:
          nt='N'
         pipe.append(nt)
      NN=''.join(pipe)
      nncount=NN.count('N')
      if nncount > 117: # filter out the markers, where 60 percent data is missing and 
                        # markers with minor allele frequency less than 5 percent  
         continue 

      out.write('\t'.join(pipe ) + '\n')
      fas.write('>'+ pipe[0] + '\n' )
      fas.write(''.join(pipe[1:]) + '\n' )
