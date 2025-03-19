from tests.utils.repertoireData import *
from tests.utils.cdr3SetGenerator import *
import numpy as np
import inspect
from tests.utils.repertoireDataFrameTransformer import repertoireDataFrameTransformer
from functools import wraps
from inspect import signature



class tcrGenerationPipeline:
    def __init__(self, numOfRepertoires):
        self.repertoires = np.array([ repertoireData() for _ in range(numOfRepertoires) ])

    def _validateGeneratorRepertoireRange(self,n,idx=None, rangeStart=None, rangeEnd=None ):
        if idx is not None and (rangeStart is not None or rangeEnd is not None):
            raise ValueError("Index cannot be used alongside range")
        if idx is None and ((rangeStart is not None) ^ (rangeEnd is not None)):
            raise ValueError("If using range, both range thesholds need to be provvided")
        # To do: add checks if range thresholds are in repertoires range limits

    def _validateGeneratorParameters(self, n, frequencies, tcrA=None, tcrB=None, repertoireRange=None):
        if repertoireRange is None :
            raise ValueError("repertoireRange cannot be None")

        #Does substituting both make any sense?
        for tcr, tcrName, tcrInverseName in zip([tcrA,tcrB],["tcrA","tcrB"],["tcrB","tcrA"]):
            if (tcr is not None and 
                    (not isinstance(tcr,np.ndarray) and
                     not inspect.ismethod(tcr))
                ) :
                raise TypeError(f"Alternate {tcrName} sequence has to be a list or numpy array")

            if (tcr is not None and
                not inspect.ismethod(tcr) and
                    (len(np.array(tcr).shape) != 1 or
                     np.array(tcr).shape[0] != n)
                ):
                raise ValueError(f"lenghts of generated {tcrInverseName} and alternative {tcrName} sequence must be identical")

        if (not inspect.ismethod(frequencies) and
            not isinstance(frequencies, list) and
            not isinstance(frequencies, np.ndarray)):
            raise TypeError("""Frequencies must be one of:
                                * method
                                * array of size n
                                * list of size of repertoire range of numpy arrays of size n, 
                                * list of size of repertoire range of methods
                            """)
        if (isinstance(frequencies,np.ndarray) and
                (len(frequencies.shape) != 1 or
                 frequencies.shape[0] != n)
            ):
            raise ValueError("If frequencies are in an array they have to be one dimensional and have lenght equal n")

        if (isinstance(frequencies,list) and
            len(frequencies) != len(repertoireRange)):
            raise ValueError("If frequencies are in a list, the list has to be an equal size to the repertoire range")
            

        if (isinstance(frequencies,list) and 
            not all([ isinstance(arr,type(frequencies[0])) for arr in frequencies])):
            raise ValueError("If frequencies are a list, then all elements have to be of the same type")


        if (isinstance(frequencies,list) and
            isinstance(frequencies[0],np.ndarray) and
            not all([ (len(arr.shape) == 1 and arr.shape[0] == n) for arr in frequencies])):
            raise ValueError("If frequencies are a list of arrays then all arrays must be one dimensional and have zie equal n")


    def toDataFrame(self):
        transform = repertoireDataFrameTransformer().transform
        return [ transform(repertoire) for repertoire in self.repertoires]

    def _assingWorkingRange(self, idx = None, rangeStart=None, rangeEnd=None):
        if idx is not None:
            return self.repertoires[idx:idx+1]
        elif rangeStart is not None and rangeEnd is not None:
            return self.repertoires[rangeStart:rangeEnd+1]
        else:
            return self.repertoires

    def _appendToTcr(self, repertoire, defaultGenerator, n, repIdx, repertoireWorkingRange, receptorInput, cloneStatus, repertoireEntryName):
        repertoireEntry = getattr(repertoire,repertoireEntryName)
        # Clone first entry if cloning is enable for given tcr
        if cloneStatus == True and repIdx != 0 :
            setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,getattr(repertoireWorkingRange[0],repertoireEntryName)[-n:])))
            return

        # Checking for alternate version of cdr3 sequences, they can be passed either as a set of sequences (np.array) or function used the generate a set 

        # if an alternate version is a method use it to generate cdr3
        if inspect.ismethod(receptorInput) :
            setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,receptorInput(n))))
            return

        # if an alternate version is a list or array then just concatenate it
        if receptorInput is not None:
            setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,receptorInput)))
            return

        # if no alternate was provide it generate with default function
        if isinstance(defaultGenerator,np.ndarray):
            setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,defaultGenerator)))
            return

        setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,defaultGenerator(n))))

    def _makeInclusiveClone(self, n, repertoire, cloneTcrInclusive):
            if cloneTcrInclusive == True:
                repertoire.tcrB[-n:] = repertoire.tcrA[-n:].copy()
            elif cloneTcrInclusive == False:
                repertoire.tcrA[-n:] = repertoire.tcrB[-n:].copy()

    def _generateRightFrequencies(self, n, repIdx, frequencies):
            if inspect.ismethod(frequencies) :
                return frequencies(n)
            elif isinstance(frequencies,np.ndarray) :
                return frequencies
            elif isinstance(frequencies, list) and isinstance(frequencies[0],np.ndarray):
                return frequencies[repIdx]
            elif isinstance(frequencies, list) and inspect.ismethod(frequencies[0]):
                return frequencies[repIdx](n)
            return np.array([])


    def generateRepertoire(
            self,
            n,
            frequencies,
            defaultGenerator,
            idx=None, 
            rangeStart=None, 
            rangeEnd=None,
            tcrA=None,
            tcrB=None,
            cloneTcrA=False,
            cloneTcrB=False,
            cloneTcrInclusive=None,
            split=False):

        self._validateGeneratorRepertoireRange(n, idx, rangeStart, rangeEnd)
        
        # Assign correct working range
        repertoireWorkingRange = self._assingWorkingRange(idx,rangeStart,rangeEnd)

        self._validateGeneratorParameters(n,frequencies,tcrA,tcrB, repertoireWorkingRange)

        defaultDataAlpha = np.array([])
        defaultDataBeta = np.array([])
        if split:
            defaultDataAlpha = defaultGenerator(n*len(repertoireWorkingRange))
            defaultDataAlpha = [ defaultDataAlpha[i*n:i*n+n] for i in range(len(repertoireWorkingRange))]
            defaultDataBeta = defaultGenerator(n*len(repertoireWorkingRange))
            defaultDataBeta = [ defaultDataBeta[i*n:i*n+n] for i in range(len(repertoireWorkingRange))]


        # Iterate through repertoires in range
        for repIdx, repertoire in enumerate(repertoireWorkingRange):
            # consider the entry alternative tcr entry, cloning status and entry in repertoire for given receptor type
            for receptorInput,cloneStatus, repertoireEntryName,defaultData  in zip([tcrA,tcrB], [cloneTcrA,cloneTcrB], ["tcrA", "tcrB"],[defaultDataAlpha,defaultDataBeta]):
                self._appendToTcr(repertoire=repertoire,
                                  defaultGenerator=defaultGenerator if not split else defaultData[repIdx],
                                  n=n,
                                  repIdx=repIdx,
                                  repertoireWorkingRange=repertoireWorkingRange,
                                  receptorInput=receptorInput,
                                  cloneStatus=cloneStatus,
                                  repertoireEntryName=repertoireEntryName)
                
            
            # handling incluse tcr cloning
            self._makeInclusiveClone(n, repertoire,cloneTcrInclusive) 
            
            # handling frequency generation
            addedFrequencies = self._generateRightFrequencies(n,repIdx,frequencies)
            repertoire.frequency = np.concatenate((repertoire.frequency, addedFrequencies))

        return self
    
    
    def generateRandomRepertoire(
            self,
            n,
            frequencies,
            idx=None, 
            rangeStart=None, 
            rangeEnd=None,
            tcrA=None,
            tcrB=None,
            cloneTcrA=False,
            cloneTcrB=False,
            cloneTcrInclusive=None):
        return self.generateRepertoire(
                n=n,
                frequencies=frequencies,
                defaultGenerator=cdr3SetGenerator().generateRandomSubset,
                idx=idx, 
                rangeStart=rangeStart, 
                rangeEnd=rangeEnd,
                tcrA=tcrA,
                tcrB=tcrB,
                cloneTcrA=cloneTcrA,
                cloneTcrB=cloneTcrB,
                cloneTcrInclusive=cloneTcrInclusive)

    def generateRelatedRepertoire(
            self,
            n,
            frequencies,
            distance,
            distanceFun,
            idx=None, 
            rangeStart=None, 
            rangeEnd=None,
            tcrA=None,
            tcrB=None,
            cloneTcrA=False,
            cloneTcrB=False,
            cloneTcrInclusive=None,
            split=False):
        def wrappedRelatedSubsetGenerator(n: int) -> np.ndarray:
            return cdr3SetGenerator().generateRelatedSubset(n=n,distance=distance,distanceFun=distanceFun)

        return self.generateRepertoire(
                n=n,
                frequencies=frequencies,
                defaultGenerator=wrappedRelatedSubsetGenerator,
                idx=idx, 
                rangeStart=rangeStart, 
                rangeEnd=rangeEnd,
                tcrA=tcrA,
                tcrB=tcrB,
                cloneTcrA=cloneTcrA,
                cloneTcrB=cloneTcrB,
                cloneTcrInclusive=cloneTcrInclusive,
                split=split)


