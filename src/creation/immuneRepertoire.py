class immuneRepertoire:
    def __init__(self, clones, sampleIDs):
        self.clones = clones
        self.sampleIDs = sampleIDs
        if not "sampleID" in clones.columns:
            raise ValueError('Clonotypes must have sampleID column')
    
