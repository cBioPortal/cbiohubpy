class GenomicVariant:
    def __init__(self, chrom, start, end, ref, alt):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.ref = ref
        self.alt = alt


class ProteinVariant:
    def __init__(self, gene, protein_change):
        self.gene = gene
        self.protein_change = protein_change if protein_change.startswith("p.") else f"p.{protein_change}"