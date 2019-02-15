import re
import sys
import math
import numpy as np

verbose = False

def eprint(*args, **kwargs):
    if verbose:
        print(*args, file=sys.stderr, **kwargs)

def read_gap_costs(filename):
    gap_costs = {}
    with open(filename) as fin:
        for line in fin:
            fields = line.rstrip().split()
            gap_costs[fields[0]] = [ int(x) for x in fields[1:] ] 
    if ( gap_costs['tablesize'][0] != len(gap_costs['qGap']) or
         gap_costs['tablesize'][0] != len(gap_costs['tGap']) or
         gap_costs['tablesize'][0] != len(gap_costs['bothGap']) ):
             print("Wrong format for the gap cost file %s" %(filename))
             exit(1)
    return gap_costs

def slope(y2,y1,x2,x1):
    # Calculate slope of line from x1/y1 to x2/y2 
    return (y2-y1)/(x2-x1)
    
class gapCalc(object):
    """
        an axtChain gap scoring schement (see src/lib/gapCalc.c in
        kentUtils)
    """
    def __init__(self, filename):
        gap_costs = read_gap_costs(filename)
        self.__tablesize = gap_costs['tablesize'][0]
        self.__smallSize  = gap_costs['smallSize'][0]
        self.__position = gap_costs['position']
        self.__qGap = gap_costs['qGap']
        self.__tGap = gap_costs['tGap']
        self.__bGap = gap_costs['bothGap']
        self.__startlong  = self.__position.index(self.__smallSize)
        
        self.__init_small()
        
        self.__qLastPos = self.__position[-1]
        self.__tLastPos = self.__position[-1]
        self.__bLastPos = self.__position[-1]
        self.__qLastPosVal = self.__qGap[-1];
        self.__tLastPosVal = self.__tGap[-1];
        self.__bLastPosVal = self.__bGap[-1];

        self.__qLastSlope = slope(self.__qLastPosVal, 
                                      self.__qGap[-2],
                                      self.__qLastPos,
                                      self.__position[-2])
        self.__tLastSlope = slope(self.__tLastPosVal, 
                                      self.__tGap[-2],
                                      self.__tLastPos,
                                      self.__position[-2])    
        self.__bLastSlope = slope(self.__bLastPosVal, 
                                      self.__bGap[-2],
                                      self.__bLastPos,
                                      self.__position[-2])                                                                   
        
    def interpolate(self, x, v, start):
        s = self.__position
        for i in range(start,self.__tablesize):
            ss = s[i]
            if x == ss:
                return v[i]
            elif x < ss:
                ds = ss - s[i-1]
                dv = v[i] - v[i-1]
                return v[i-1] + dv * (x - s[i-1])/ds
        # If get to here extrapolate from last two values 
        ds = s[-1] - s[-2];
        dv = v[-1] - v[-2];
        return v[-2] + dv * (x - s[-2]) / ds;
    def __init_small(self):
        self.__qsmall = [0]*self.__smallSize
        self.__tsmall = [0]*self.__smallSize
        self.__bsmall = [0]*self.__smallSize
        for i in range(1,self.__smallSize):
            self.__qsmall[i] = self.interpolate(i, self.__qGap, 0)
            self.__tsmall[i] = self.interpolate(i, self.__tGap, 0)
            self.__bsmall[i] = self.interpolate(i, self.__bGap, 0)
    
    def gapCalcCost(self,dq, dt):
        # Figure out gap costs. 
        dt = 0 if dt<0 else dt
        dq = 0 if dq<0 else dq
        if dt == 0:
            if dq < self.__smallSize:
                return self.__qsmall[dq]
            elif dq > self.__qLastPos:
                return self.__qLastPosVal + self.__qLastSlope * (dq - self.__qLastPos)
            else:
                return self.interpolate(dq, self.__qGap, self.__startlong )
        elif dq == 0:
            if dt < self.__smallSize:
                return self.__tsmall[dt]
            elif dt > self.__tLastPos:
                return self.__tLastPosVal + self.__tLastSlope * (dt - self.__tLastPos)
            else:
                return self.interpolate(dt, self.__tGap, self.__startlong )
        else:
            both = int(dq + dt)
            if both < self.__smallSize:
                return self.__bsmall[both]
            elif dt > self.__tLastPos:
                return self.__bLastPosVal + self.__bLastSlope * (both - self.__bLastPos)
            else:
                return self.interpolate(both, self.__bGap, self.__startlong )

    

