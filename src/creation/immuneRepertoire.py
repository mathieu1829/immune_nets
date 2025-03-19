class immuneRepertoire:
    def __init__(self, clones):
        self.clones = clones
        if not "sampleID" in clones.columns:
            raise ValueError('Clonotypes must have sampleID column')
    
