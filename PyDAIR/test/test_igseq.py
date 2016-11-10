import unittest
from PyDAIR.utils.PyDAIRUtils import *
from PyDAIR.utils.PyDAIRArgs import *
from PyDAIR.seq.IgSeq import *


class Test_pydair_seq_IgSeq(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_forward_igseq(self):
        # forward query, forward subject
        ## query      TGTTCTCGGACTCCTATAACGGATGAATTTGACAAGTGGGGAAAGGGAAAAGCGCGGACGGGTTCGACCG
        ## V          TGTGCTCGAGAG----------------------------------------------------------
        ## cdr3       --------GACTCCTATAACGG------------------------------------------------
        ## J          ------------------TACTATGAATTTGACTACTGGGGAAAGGGCACATCGGTGACTGTTTCGTCC-
        q_name = 'M03400:8:000000000-ADYYJ:1:1101:14295:1984 1:N:0:1'
        q_seq  = 'GGTTAAAAGACCTGGAGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTATGGGCTGAACTGGGTCTGACAGGCTCCTTGAAAAGGACTGGCGTGGATTGCTTATATCTACATCTGCACATACTACTCTGAGTCATTCAAAGGCCGGTTTATCCTCTCCAGAGACAACAACATCGCACAGCTGAATCTGCAGATGCACACCCTGAAGACTGAAGTTTCTTCTGTTTATTATTGTTCTCGgactcctataacggATGAATTTGACAAGTGGGGAAAGGGAAAAGCGCGGACGGGTTCGACCGAAAAAACAAAAGACCGTTCTAAGTGGACATTGATAGAAGGAGGATCTGGGACTGGAAAAAAGGTAAAGCGCGGAGGTTTGGCAGACGACATCACGCCAAGGGAACTAACCTACACCTGGAGAAAAGACGGAGTGGATGGGAAAGACTTCATTCA'
        v_seq = 'GTGTTGATGCTCAGACTCTGACCCAGTCTGAACCAGTGGTTAAAAGACCTGGGGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTACTGGATGGTCTGGGTCAGACAGGCTCCTGGAAAAGGACTGGAGTGGATCGCTTATATCACCACCAGTAGCAGCCCATACTACTCTGAGTCAGTCAAAGGCCGGTTTATCATCTCCAGAGACAACAACAGAGCACAGCTGAATCTGCAGATTAACAGCCTGAAGACTGAAGATTCTGCTGTTTATTATTGTGCTCGAGAG'
        v_qaligned_seq = 'GGTTAAAAGACCTGGAGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTATGGGCTGAACTGGGTCTGACAGGCTCCTTGAAAAGGACTGGCGTGGATTGCTTATATCTACATCTGCA------CATACTACTCTGAGTCATTCAAAGGCCGGTTTATCCTCTCCAGAGACAACAACATCGCACAGCTGAATCTGCAGATGCACACCCTGAAGACTGAAGTTTCTTCTGTTTATTATTGTTCTCG'
        v_saligned_seq = 'GGTTAAAAGACCTGGGGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGCAGCTACTGGATGGTCTGGGTCAGACAGGCTCCTGGAAAAGGACTGGAGTGGATCGCTTATATCACCACCAGTAGCAGCCCATACTACTCTGAGTCAGTCAAAGGCCGGTTTATCATCTCCAGAGACAACAACAGAGCACAGCTGAATCTGCAGATTAACAGCCTGAAGACTGAAGATTCTGCTGTTTATTATTGTGCTCG'
        j_seq = 'TACTATGAATTTGACTACTGGGGAAAGGGCACATCGGTGACTGTTTCGTCC'
        j_qaligned_seq = 'ATGAATTTGACAAGTGGGGAAAGGGAAAAGCGCGGACGGGTTCGACC'
        j_saligned_seq = 'ATGAATTTGACTACTGGGGAAAGGGCACATCGGTGACTGTTTCGTCC'
        
        igseq_v_query = IgSeqAlignQuery( q_name, q_seq, v_qaligned_seq,  1, 255, '+')
        igseq_v_sbjct = IgSeqAlignSbjct('V1.14', v_seq, v_saligned_seq, 38, 298, '+')
        igseq_v = IgSeqAlign(igseq_v_query, igseq_v_sbjct)
        igseq_d = None
        igseq_j_query = IgSeqAlignQuery(q_name, q_seq, j_qaligned_seq, 270, 316, '+')
        igseq_j_sbjct = IgSeqAlignSbjct( 'Jm5', j_seq, j_saligned_seq,   5,  51, '+')
        igseq_j = IgSeqAlign(igseq_j_query, igseq_j_sbjct)
        
        igseq = IgSeq(igseq_v, igseq_d, igseq_j)
        igseq.seek_cdr3()
        print([igseq.variable_region.untemplate_region,
               igseq.query.seq[igseq.variable_region.untemplate_region[0]:igseq.variable_region.untemplate_region[1]]])
        ## => [[255, 269], 'gactcctataacgg']
        if igseq.query.orf is not None:
            print([igseq.variable_region.cdr3,
                   Seq(igseq.query.seq[igseq.variable_region.cdr3[0]:igseq.variable_region.cdr3[1]], generic_dna).translate()])
            ## => [[247, 286], Seq('CSRTPITDEFDKW', ExtendedIUPACProtein())]
        cdr3_data = igseq.get_cdr3_data()
        print([cdr3_data.nucl_seq, cdr3_data.prot_seq])
    
    def test_revcomp_igseq(self):
        # rev. comp. query and forward subjcet
        ## query      AACAGCCTGAAGACTGAAGAGTCTGGTGTTTATTATGGAAAGCGATCGACTACTTTTGACTAATGGGGGAAAGGAACA (rev. comp.)
        ## V          AACAGCCTGAAGACTGAAGATTCTGCTGTTTATTATTGTGCTCGA---------------------------------
        ## cdr3       ------------------------------------GGAAAGCGATCGACTACT------------------------
        ## J          ------------------------------------------TACTACGCATACTTTGACTACTGGGGGAAAGGAACA
        q_name = 'M03400:8:000000000-ADYYJ:1:1101:14012:1984 1:N:0:1'
        q_seq  = 'TGAATGAAGTCTTTCAGATCGACTCCGTCTTTTCTCCAGGTGTAGGTTAGGTCCGATGGCGTGAAGTCGGCGGCCCAACATCCGAGAGTGACCATGTTTCCAGTCCCAGATCCGCATTGTATCAATGTAAACAGAGAAGGGGCTTTGGGTGTTGCAGATGTTACTGTACCTGTTGTTCCTTTCCCCCATTAGTCAAAagtagtcgatcgctttccATAATAAACACCAGACTCTTCAGTCTTCAGGCTGTTCATCTTCAGATTCACCTGGGGTGTGTTTGGGAGGGTGGAGATGAAAAAGGGGCGTTGGAAAGAGTAAGAGTAGTAGGGACGGCGGTAGATAAAAGAAAGCCACGCCAGTCATGAACCAGGAGCCTGTATGAAGCAGTTCATCCCAAAGTCACTGAATGTGAATCCAGAACCTGAACAGGTCAGTGTGTGAGAGGGTCCAGGTCTTTTAAC'
        v_seq  = 'GTGTTGATGCTCAGACTCTGACCCAGTCTGAACCAGTGGTTAAAAGACCTGGAGAATCTCACACACTGACCTGTTCAGCCTCTGGATTCACATTCAGTGACTATGGGATGCACTGGGTCAGACAGGTTCCTGGAAAAGGACTGGAGTGGATCGCTTATATCTACAGCAGTGGCAGCAGCCAGTACTACTCTGAGTCAGTCAAAGGCCGGTTTATCATCTCCAGAGACAACAACAGAGCACAGCTGAATCTGCAGATGAACAGCCTGAAGACTGAAGATTCTGCTGTTTATTATTGTGCTCGA'
        v_qaligned_seq = 'ATAATAAACACCAGACTCTTCAGTCTTCAGGCTGTTCATCTTCAGATTCACCTGGGGTGTGTTTGGGAGGGTGGAGATGAAAAAGGGGCGTTGGAAAGAGTAAGAGTAGTA---------GGGACGGCGGTAGATAAAAGAAAGCCACGCCAGTCATGAACCAGGAGCCTGTATGAAGCAGTTCATCCCAAAGTCACTGAATGTGAATCCAGAACCTGAACAGGTCAGTGTGTGAGAGGGTCCAGGTCTTTTAAC'
        v_saligned_seq = 'ATAATAAACAGCAGAATCTTCAGTCTTCAGGCTGTTCATCTGCAGATTCAGCTGTGCTCTGTTGTTGTCTCTGGAGATGATAAACCGGCCTTTGACTGACTCAGAGTAGTACTGGCTGCTGCCACTGCTGTAGATATAAGCGATCCACTCCAGTCCTTTTCCAGGAACCTGTCTGACCCAGTGCATCCCATAGTCACTGAATGTGAATCCAGAGGCTGAACAGGTCAGTGTGTGAGATTCTCCAGGTCTTTTAAC'
        j_seq = 'TACTACGCATACTTTGACTACTGGGGGAAAGGAACAACAGTTACAGTAACATCT'
        j_qaligned_seq = 'AGATGTTACTGTACCTGTTGTTCCTTTCCCCCATTAGTCAAA'
        j_saligned_seq = 'AGATGTTACTGTAACTGTTGTTCCTTTCCCCCAGTAGTCAAA'
        
        igseq_v_query = IgSeqAlignQuery( q_name, q_seq, v_qaligned_seq, 216, 461, '+')
        igseq_v_sbjct = IgSeqAlignSbjct('V1.21', v_seq, v_saligned_seq,  39, 293, '-')
        igseq_v = IgSeqAlign(igseq_v_query, igseq_v_sbjct)
        igseq_d = None
        igseq_j_query = IgSeqAlignQuery(q_name, q_seq, j_qaligned_seq, 156, 197, '+')
        igseq_j_sbjct = IgSeqAlignSbjct( 'Jm5', j_seq, j_saligned_seq,  13,  54, '-')
        igseq_j = IgSeqAlign(igseq_j_query, igseq_j_sbjct)
        
        igseq = IgSeq(igseq_v, igseq_d, igseq_j)
        igseq.seek_cdr3()
        
        print([igseq.variable_region.untemplate_region,
               igseq.query.seq[igseq.variable_region.untemplate_region[0]:igseq.variable_region.untemplate_region[1]]])
        
        if igseq.query.orf is not None:
            print([igseq.variable_region.cdr3,
                   Seq(igseq.query.seq[igseq.variable_region.cdr3[0]:igseq.variable_region.cdr3[1]], generic_dna).translate()])
        cdr3_data = igseq.get_cdr3_data()
        
        print([cdr3_data.nucl_seq, cdr3_data.prot_seq])
        cdr3_data = igseq.get_cdr3_data(v_adj = -12, j_adj = 9)
        print([cdr3_data.nucl_seq, cdr3_data.prot_seq])
    
    

if __name__ == '__main__':
    unittest.main()


