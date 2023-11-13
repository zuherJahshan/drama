from typing import List, Dict
import cupy as cp
import os
import math
from reads import Reads


def one_hot_encode(idx, num_classes):
    one_hot = cp.zeros(num_classes, dtype=cp.int32)
    one_hot[idx] = 1
    return one_hot


def reverseComplement(sequence):    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))



class Genome:
    def __init__(self, filepath):
        if os.path.exists(filepath):
            self.filepath = filepath
        else:
            raise FileNotFoundError(f"File {filepath} not found")
        
        self.name = self.filepath.split("/")[-1].split(".")[0]
        self._build_genome()


    def kmers(self, kmer_size):
        # read file to a buffer,
        for contig in self.genome:
            for seq in [contig, reverseComplement(contig)]:
                for idx in range(len(seq) - kmer_size + 1):
                    yield seq[idx:idx+kmer_size]
    
    ########################
    #### Private Methods ###
    ########################
    def _build_genome(self):
        self.genome = []
        with open(self.filepath, "r") as f:
            contig = ""
            for idx, line in enumerate(f):
                if line[0] == ">" and idx != 0:
                    self.genome.append(contig)
                    contig = ""
                elif line[0] != ">":
                    contig += line.strip()
            self.genome.append(contig)



class DramaSimulator:
    def __init__(
        self,
        kmer_size = 32,
        et = 0,
        faulty_cells = 0,
    ):
        self.columns:                    int              = 2 ** 13
        self.banks:                      int              = 2 ** 3
        self.kmer_size:                  int              = kmer_size
        self.rows:                       int              = int(2 ** 14 / (kmer_size * 2 * 2))
        self.current_position:           Dict              = {
            "row": 0,
            "bank": 0,
            "column": 0
        }

        self.et = et

        self.kmers = {}
        """
        self.kmers = {
            kmer: {
                "probabilities": [],
                "genomes": [], // onehot encoding of their index
            }
        }
        """

        self.genomes_to_indices: Dict[str, int] = {}
        self.indices_to_genomes: Dict[int, str] = {}

        # faulty cells
        self.faulty_cells = {}
        number_of_faulty_cells = faulty_cells
        for bank in range(self.banks):
            bad_columns = cp.random.randint(0, self.columns, number_of_faulty_cells)
            bad_rows = cp.random.randint(0, self.rows, number_of_faulty_cells)
            fault_constants = cp.maximum((cp.random.randn(number_of_faulty_cells) / 10) + 0.5, 0)
            self.faulty_cells.update({bank: {
                "bad_columns": bad_columns,
                "bad_rows": bad_rows,
                "fault_constants": fault_constants
            }})

        self.reads:                     Reads            = Reads("./data/", 150)


    def insert_ref(self, genomes: List[Genome]):
        for genome_idx, genome in enumerate(genomes):
            self.genomes_to_indices[genome.name] = genome_idx
            self.indices_to_genomes[genome_idx] = genome.name
        
        for genome in genomes:
            self._insert_genome(genome)

        self.build_kmers()


    def search_reads(
        self,
        name,
        genome_filepath,
        sequencer
    ):
        results = {
            "tp": 0,
            "fp": 0,
            "fn": 0,
        }
        reads = self.reads.getReads(sequencer, genome_filepath)[:1000]
        for read in reads:
            read_idx_results = self._search_read(read, name)
            results["tp"] += read_idx_results["tp"]
            results["fp"] += read_idx_results["fp"]
            results["fn"] += read_idx_results["fn"]
        return results


    def get_cells_heatmap(self, bank, number_of_rows_to_include = None, number_of_columns_to_include = None):
        if number_of_rows_to_include is None:
            number_of_rows_to_include = self.rows
        if number_of_columns_to_include is None:
            number_of_columns_to_include = self.columns
        heatmap = cp.zeros((self.rows, self.columns))
        for row in range(number_of_rows_to_include):
            print(f"=> Started row {row}")
            for column in range(number_of_columns_to_include):
                position = {
                    "row": row,
                    "bank": bank,
                    "column": column
                }
                heatmap[row, column] = self._calculate_faulty_cell_prob(position)
            print(f"===> Finished row {row}")
        return heatmap



    ########################
    #### Private Methods ###
    ########################

    def _search_read(self, read: str, name: str):
        num_kmers = len(read) - self.kmer_size + 1
        genome_idx_results = cp.zeros(len(self.genomes_to_indices))
        for idx in range(num_kmers):
            kmer = read[idx:idx+self.kmer_size]
            genome_idx_results = genome_idx_results + self._search_kmer(kmer)
                
        # If all is zeros return None
        # print(self.indices_to_genomes)
        tp, fp, fn = 0, 0, 0
        if cp.sum(genome_idx_results) == 0:
            fn += 1
        else:
            for idx, cond in enumerate(genome_idx_results == cp.max(genome_idx_results)):
                if cond:
                    if self.indices_to_genomes[idx] == name:
                        tp += 1
                    else:
                        fp += 1
        return {
            "tp": tp,
            "fp": fp,
            "fn": fn
        }



    def _insert_genome(self, genome):
        consider_lower = False
        if math.ceil(self.kmer_size / (self.et + 1)) > math.floor(self.kmer_size / (self.et + 1)):
            consider_lower = True
        pattern_size = math.ceil(self.kmer_size / (self.et + 1))
        for pattern in genome.kmers(pattern_size):
            self._insert_kmer(pattern, genome)
            if consider_lower:
                self._insert_kmer(pattern[:-1], genome)
            self._increment_position()

    
    def _insert_kmer(self, kmer, genome):        
        if kmer in self.kmers:
            self.kmers[kmer]["probabilities"].append(self._calculate_faulty_cell_prob())
            self.kmers[kmer]["genomes"].append(
                one_hot_encode(
                    self.genomes_to_indices[genome.name],
                    len(self.genomes_to_indices)
                )
            )
        else:
            self.kmers[kmer] = {
                "probabilities": [self._calculate_faulty_cell_prob()],
                "genomes": [
                    one_hot_encode(
                        self.genomes_to_indices[genome.name],
                        len(self.genomes_to_indices)
                    )
                ]
            }


    def _increment_position(self, position = None):
        if position == None:
            pos = self.current_position
        else:
            pos = position
        pos["column"] += 1
        pos["column"] %= self.columns
        pos["bank"] += 1 if self.current_position["column"] == 0 else 0
        pos["bank"] %= self.banks
        pos["row"] += 1 if self.current_position["column"] == 0 and self.current_position["bank"] == 0 else 0
        pos["row"] %= self.rows
        return pos


    def _calculate_faulty_cell_prob(self, position = None):
        if position is None:
            pos = self.current_position
        else:
            pos = position
        denomenator = 1 +\
              (pos["row"] - self.faulty_cells[pos["bank"]]["bad_rows"]) ** 2 +\
                (pos["column"] - self.faulty_cells[pos["bank"]]["bad_columns"]) ** 2
        fault_prob_vec = self.faulty_cells[pos["bank"]]["fault_constants"] / denomenator
        fault_prob = cp.sum(fault_prob_vec)
        return cp.minimum(fault_prob, 1)
    

    def build_kmers(self):
        for kmer in self.kmers:
            self.kmers[kmer]["probabilities"] = cp.array(self.kmers[kmer]["probabilities"])
            self.kmers[kmer]["genomes"] = cp.array(self.kmers[kmer]["genomes"])


    def _search_kmer(self, kmer):
        # Devide the kmer into pigeonholes (partial kmers)
        pigeonholes = self.et + 1
        partial_kmer_lengths = [int(math.floor(self.kmer_size / (self.et + 1)))] * pigeonholes
        for i in range(self.kmer_size % pigeonholes):
            partial_kmer_lengths[i] += 1
        
        # Calculate match results for the pigeonholes
        kmer_idx_results = 0
        partial_kmer_idx = 0
        for partial_kmer_length in partial_kmer_lengths:
            partial_kmer = kmer[partial_kmer_idx: partial_kmer_idx + partial_kmer_length]
            if partial_kmer in self.kmers:
                appearing = (cp.random.rand(len(self.kmers[partial_kmer]["probabilities"])) >= self.kmers[partial_kmer]["probabilities"])[:, cp.newaxis]
                kmer_idx_results = kmer_idx_results + (cp.sum(appearing * self.kmers[partial_kmer]["genomes"], axis=0) > 0)
        kmer_idx_results = kmer_idx_results > 0
        return kmer_idx_results

            



