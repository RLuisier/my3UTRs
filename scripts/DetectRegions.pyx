#If you modify the script, you need to recompile it as it follows:
# 1/Prepare the setup.py file as follows: 
#
#     from distutils.core import setup
#     from distutils.extension import Extension
#     from Cython.Distutils import build_ext
#     setup(
#      cmdclass = {'build_ext':build_ext},
#      ext_modules = [Extension("DetectRegions", ["DetectRegions.pyx"])]
#     )
# 2/In the same directory: python setup.py build_ext --inplace
#
import numpy as np
def DetectRegions(float lim, int Lim, int winS, toRead, int comment):
    print 'Start routine'
    #Open file
    data = np.genfromtxt(fname=toRead,dtype='int')
    #Declare C local variables
    cdef int noLines  = len(data)-1
    cdef int start    = data[0]#Read first line of the file
    cdef int Cond     = 0 #Check whether previous sliding window fullfilled the condition.
    cdef int deltaPos = 1
    cdef int tempEnd  = data[0]
    cdef int idx      = -1 #Index of the start of the sliding window
    cdef int ix       
    #Create lits that are going to be returned
    Start=[] #Pas certaine de pouvoir definir une list python en variable c
    End=[]
    while (idx < noLines): #Go as long as idx is smaller than number of lines of the file
        if idx%1000000==0:
            print(str(idx/1000000)+'M...');
        idx       = idx+1;
        ix        = idx;#Index in the sliding window (always higher or equal to idx)
        npos      = 1; #Number of positions covered by more than 10 reads in the sliding window.
        start     = data[idx];
        end       = data[idx]; #current location of the pointer
        deltaPos  = 1; #Size of the window covered
        #Start small block:
        if comment==1:
            print '(re)start block'
        while (deltaPos < winS)&(ix<noLines): #This loops over a npos number of position until deltapos<winS
            ix        = ix+1;     #Increment index of the sliding window
            npos      = npos + 1; #Increment the number of posititions covered
            end       = data[ix];
            deltaPos  = end - start;
        else:    #Goes out of the while loop as soon as deltapos is the size of the sliding window
            if  (deltaPos>=winS) & (npos >= Lim):
                if Cond==0:#Then I have to start a new block
                    Start.append(start); #Need to store the start as start new block
                    Cond=1;
                    if (float(npos)/deltaPos >= lim):
                        tempEnd    = end;#Store the end in a temporary variable
                    else:
                        tempEnd    = start+100;#Store the end in a temporary variable
                else:
                    if (float(npos)/deltaPos >= lim):
                        if (tempEnd < start): #Then I can close previous block and start new one
                            End.append(tempEnd);#Close previous block
                            Start.append(start);#Open new block
                            tempEnd=end;#Store current end
                        else: #Then I can continue previous one
                            tempEnd=end;#Store current end
                    else:
                        if (tempEnd < start): #Then I can close previous block and start new one
                            End.append(tempEnd);#Close previous block
                            Start.append(start);#Open new block
                            tempEnd=start+100;#Store current end
                        else: #Then I can continue previous one
                            tempEnd=start+100;#Store current end  
            elif (npos < Lim) & (Cond==1):#Then I need to close the previous block
                End.append(tempEnd);
                Cond=0;
    else: #I am at the end of the file               
        if (Cond==1):
            print 'enter final loop by adding last end from current window.'
            End.append(tempEnd);
        else:
            print 'end final loop by doing nothing (previous block close and current window not fullfilling conditions.)'
    out = np.array([Start, End]).transpose()
    return out 

