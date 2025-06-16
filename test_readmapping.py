import unittest
import os
from tempfile import NamedTemporaryFile
import readmapping

class TestReadMapping(unittest.TestCase):

    def setUp(self):
        # Create a temporal file
        self.test_fasta_data = """>seq1
        ATGCATGC
        >seq2
        CGTACGTA
        """
        self.test_file = NamedTemporaryFile(delete=False, mode='w', newline='', suffix= '.fa')
        self.test_file.write(self.test_fasta_data)
        self.test_file.close()

    def test_check_dna(self):
        valid_dna = "ATGCATGC"
        invalid_dna = "ATGCATXZ"
        # Check the DNA sequence
        try:
            readmapping.check_dna(valid_dna)
        except ValueError:
            self.fail("check_dna raised ValueError unexpectedly!")
            
        with self.assertRaises(ValueError):
            readmapping.check_dna(invalid_dna)

    def test_read_fasta_valid(self):
        # Check the read fasta 
        result = readmapping.read_fasta(self.test_file.name)
        self.assertEqual(len(result), 2)  # We have two sequences
        self.assertIn('seq1', result)
        self.assertEqual(result['seq1'], "ATGCATGC")

    def test_get_kmers(self):
        seq = "ATGCATGC"
        expected_kmers = {
            "ATG": [0, 4],
            "TGC": [1, 5],
            "GCA": [2],
            "CAT": [3]}
        kmers = readmapping.get_kmers(seq, 3)
        self.assertEqual(kmers, expected_kmers)
        # Test raising the error
        with self.assertRaises(ValueError):
            readmapping.get_kmers(seq, 15)
            
    def test_kmer_aligment(self):
        seq = 'ATA'
        # We test
        # - Sequences that align twice to a same reference (ref1)
        # - Sequences that do not align to a reference (ref2)
        # - Sequences that align at the very start (ref1) or end of the ref  (ref3)
        ref = {"ref1": "ATGCATAC",
               "ref2": "CGCGCGCG",
               "ref3": "CGCGCATA"}
        result = readmapping.kmers_alignment(seq, ref, 2)
        expected_result = {"ref1": [0, 4],
               "ref2": [],
               "ref3": [5]}
        self.assertEqual(result, expected_result)

    def test_hamming_distance(self):
        # Testing d = 0
        seq1 = "ATGC"
        seq2 = "ATGC"
        self.assertEqual(readmapping.hamming_distance(seq1, seq2), 0)
        # Testing d = 1
        seq1 = "ATGC"
        seq2 = "ATCC"
        self.assertEqual(readmapping.hamming_distance(seq1, seq2), 1)
        # Test raising the ValueError
        seq1 = "ATGC"
        seq2 = "ATCCA"
        with self.assertRaises(ValueError):
            readmapping.hamming_distance(seq1, seq2)
        

    def test_best_mapping_location(self):
        seq = 'ATA'
        # Test only 1 full alignment (seq2)
        refs_1 = {'ref1': 'ATGCATGC',
                'ref2': 'CGATACGC'}
        expected_1 = {'ref2' : [(2,0)]}
        res_1 = readmapping.best_mapping_location(seq, refs_1, 2)
        self.assertEqual(res_1, expected_1)
        # Test two best locations in different refs
        refs_2_dif = {'ref1': 'ATACATGC',
                'ref2': 'CGATACGC'}
        expected_2_dif = {'ref1' : [(0,0)], 'ref2': [(2,0)]}
        res_2_dif = readmapping.best_mapping_location(seq, refs_2_dif, 2)
        self.assertEqual(res_2_dif, expected_2_dif)
        # Test two best locations in the same ref
        refs_2_same = {'ref1': 'ATACATAC',
                'ref2': 'CGATACGC'}
        expected_2_same = {'ref1' : [(0,0), (4,0)], 'ref2': [(2,0)]}
        res_2_same = readmapping.best_mapping_location(seq, refs_2_same, 2)
        self.assertEqual(res_2_same, expected_2_same)
        # Test best location with errors
        refs_errors = {'ref1': 'ATTCCTCC',
                'ref2': 'CGTTACGC'}
        expected_errors = {'ref1': [(0, 1)], 'ref2' : [(2, 1)]}
        res_errors = readmapping.best_mapping_location(seq, refs_errors, 2)
        self.assertEqual(res_errors, expected_errors)
        # Test no coincide
        refs_no = {'ref' : 'CCCCCCCC'}
        expected_no = {}
        res_no = readmapping.best_mapping_location(seq, refs_no, 2)
        self.assertEqual(res_no, expected_no)

    def test_get_data_for_output(self):
        # We are testing
        # -Sequences that align twice to the same reference (read1)
        # -Sequences that have different kmer aligments but only one is the best one (read2)
        # -Sequences that do not align to any reference (read3)
        # -Refences covered by the same read more than once (ref1)
        # -References with no coverage (ref2)
        # -References covered by multiple reads (ref3)
        
        reads = {'read1': 'ATG',
                 'read2': 'ATA',
                 'read3': 'GGG'}
        refs = {'ref1': 'ATCCATCC',
                'ref2' : 'CCCCCCC',
                'ref3' : 'TTTATATT'}
        file1, file2, plot = readmapping.get_data_for_output(reads, refs, 2)
        ex_file1 = {'read1': {'ref1': [(0,1), (4,1)], 'ref3':[(3,1),(5,1)]},
                    'read2': {'ref3': [(3,0)]},
                    'read3': {}}
        ex_file2 = {'ref1': [('read1', 3, 8), ('read1', 3, 8)],
                    'ref2': [],
                    'ref3': [('read1', 3, 8),('read1', 3, 8) ,('read2', 3, 8)]}
        ex_plot = {'read1':4, 'read2':1, 'read3':0}
        self.assertEqual(file1, ex_file1)
        self.assertEqual(file2, ex_file2)
        self.assertEqual(plot, ex_plot)
        
    def test_write_output_file1(self):
        file1_d = {'read1': {'ref1': [(0,1), (4,1)], 'ref3':[(3,1),(5,1)]},
                    'read2': {'ref3': [(3,0)]},
                    'read3': {}}
        
        expected_output = """read1,ref1,0,1\nread1,ref1,4,1\nread1,ref3,3,1\nread1,ref3,5,1\nread2,ref3,3,0\n"""
        
        readmapping.write_output_file_1(file1_d, 'file')
            
        with open('file_1.txt') as f:
            actual_output = f.read()
                
        
        self.assertEqual(actual_output, expected_output)
        
    def test_calculate_coverage(self):
        file2_d = {'ref1': [('read1', 3, 8), ('read1', 3, 8)],
                    'ref2': [],
                    'ref3': [('read1', 3, 8),('read1', 3, 8) ,('read2', 3, 8)]}
        plot_d = {'read1':4, 'read2':1, 'read3':0}
        
        # Calculate the coverage
        # ref1 = 3 * (2/4)/8 = 0.1875
        # ref2 = 0
        # ref3 = 3 * (2/4)/8 + 3 * (1/1)/8 = 0.5625
        
        ex_output = "ref1,0.1875\nref2,0\nref3,0.5625\n"
        
        readmapping.calculate_coverage(file2_d, plot_d, 'file')
        with open('file_2.txt') as f:
            actual_output = f.read()
            
        self.assertEqual(ex_output, actual_output)


    def tearDown(self):
        # Remove the temporal files
        os.remove(self.test_file.name)

if __name__ == '__main__':
    unittest.main()


