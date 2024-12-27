from tests.utils.repertoireData import *
from tests.utils.cdr3SetGenerator import *
import numpy as np
import inspect
from tests.utils.repertoireDataFrameTransformer import repertoireDataFrameTransformer


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

        self._validateGeneratorRepertoireRange(n, idx, rangeStart, rangeEnd)
        
        # Assign correct working range
        if idx is not None:
            repertoireWorkingRange = self.repertoires[idx:idx+1]
        elif rangeStart is not None and rangeEnd is not None:
            repertoireWorkingRange = self.repertoires[rangeStart,rangeEnd+1]
        else:
            repertoireWorkingRange = self.repertoires

        self._validateGeneratorParameters(n,frequencies,tcrA,tcrB, repertoireWorkingRange)

        # Iterate through repertoires in range
        for repIdx, repertoire in enumerate(repertoireWorkingRange):
            # consider the entry alternative tcr entry, cloning status and entry in repertoire for given receptor type
            for receptorInput,cloneStatus, repertoireEntryName  in zip([tcrA,tcrB], [cloneTcrA,cloneTcrB], ["tcrA", "tcrB"]):
                repertoireEntry = getattr(repertoire,repertoireEntryName)
                # Clone first entry if cloning is enable for given tcr
                if cloneStatus == True and repIdx != 0 :
                    setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,getattr(repertoireWorkingRange[0],repertoireEntryName)[-n])))
                    continue

                # Checking for alternate version of cdr3 sequences, they can be passed either as a set of sequences (np.array) or function used the generate a set 

                # if an alternate version is a method use it to generate cdr3
                if inspect.ismethod(receptorInput) :
                    setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,receptorInput(n))))
                    continue

                # if an alternate version is a list or array then just concatenate it
                if receptorInput is not None:
                    setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,receptorInput)))
                    continue

                # if no alternate was provide it generate with default function
                setattr(repertoire, repertoireEntryName, np.concatenate((repertoireEntry,cdr3SetGenerator().generateRandomSubset(n))))

            
            if cloneTcrA and cloneTcrB and repIdx != 0 :
                continue
            
            # handling incluse tcr cloning
            if cloneTcrInclusive == True:
                repertoire.tcrB[-n:] = repertoire.tcrA[-n:].copy()
            elif cloneTcrInclusive == False:
                repertoire.tcrA[-n:] = repertoire.tcrB[-n:].copy()
            
            # handling frequency generation
            addedFrequencies = np.array([])
            if inspect.ismethod(frequencies) :
                addedFrequencies = frequencies(n)
            elif isinstance(frequencies,np.ndarray) :
                addedFrequencies = frequencies
            elif isinstance(frequencies, list) and isinstance(frequencies[0],np.ndarray):
                addedFrequencies = frequencies[repIdx]
            elif isinstance(frequencies, list) and inspect.ismethod(frequencies[0]):
                addedFrequencies = frequencies[repIdx](n)
            repertoire.frequency = np.concatenate((repertoire.frequency, addedFrequencies))



        return self

