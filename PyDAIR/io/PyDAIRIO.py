import re
import os
import sys
import math
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.utils.PyDAIRUtils import *
from PyDAIR.seq.IgSeq import *


class PyDAIRIO: 
    """PyDAIRIO class.
    
    This class implements the input/output methods to handle PYDAIR format file.
    """
    
    def __init__(self, pydair_file_path, open_mode = None):
        """PyDAIRIO class initialize method.
        
        Args:
            pydair_file_path (str): A file path to open or write.
            open_mode (str): A character (`r`, `w` or `a`) to specify open mode.
        
        Set up file path, open mode, and file format.
        After setting, the method opens the file hanld.
        """
        
        if open_mode is None:
            raise ValueError('Open mode is required. Please set open mode as \'r\' for reading mode, or \'w\' or \'a\' for writting mode.')
        self.__path   = pydair_file_path
        self.__format = 'pydair'
        self.__mode   = open_mode
        self.__fh     = open(self.__path, self.__mode)
        self.__utils  = PyDAIRUtils()
        
    
    def close(self):
        """Close file handle.
        
        Close file handle that is opened by `__init__`.
        """
        
        self.__fh.close()
    
    
    
    def write(self, igseq):
        """Write **IgSeq** objects to a file.
        
        Args:
            igseq (list): A list of IgSeq objects.
        
        Write a list of **IgSeq** object to a file.
        It is expected to use `close` method to close the file handle after writting.
        """
        
        if (not isinstance(igseq, list)) and (not isinstance(igseq, tuple)):
            igseq_list = [igseq]
        else:
            igseq_list = igseq
        
        for igseq_single in igseq_list:
            self.__fh.write(self.__convert_to_pydair(igseq_single))
    
    
    def __convert_to_pydair(self, igseq):
        
        tmpl = '''#BEGIN
QN %s
VN %s
DN %s
JN %s
OP %s
QA %s
VA %s
JA %s
UA %s
CA %s
VD %s
JD %s
VJ %s
#CDR3AA %s
#AL QSTART\tQEND\tSSTART\tSEND\tIDENTITY\tSCORE
AL QV %s\t%s\t%s\t%s\t%s\t%s
AL QD %s\t%s\t%s\t%s\t%s\t%s
AL QJ %s\t%s\t%s\t%s\t%s\t%s
AL QU %s\t%s\t%s\t%s\t%s\t%s
AL QC %s\t%s\t%s\t%s\t%s\t%s
#END
'''
        
        alignment = igseq.print_alignment()
        
        metadata_list = []
        
        # gene name section
        metadata_list.append(igseq.query.name if igseq.query is not None else None)
        metadata_list.append(igseq.v.sbjct.name if igseq.v is not None else None)
        metadata_list.append(igseq.d.sbjct.name if igseq.d is not None else None)
        metadata_list.append(igseq.j.sbjct.name if igseq.j is not None else None)
        
        # query ORF section
        metadata_list.extend([igseq.query.orf] if igseq.query is not None else [None, None])
        
        # alignment section
        metadata_list.extend(alignment)
        
        # indels section
        metadata_list.extend([igseq.indels.v_deletion, igseq.indels.j_deletion, igseq.indels.vj_insertion] if igseq.indels is not None else [None, None, None])
        
        # CDR3AA section
        metadata_list.append(str(Seq(alignment[4].replace('-', '').replace(' ', ''), generic_dna).translate()))
        
        # alignment stats section
        if igseq.v is not None:
            metadata_list.extend([igseq.v.query.start, igseq.v.query.end,
                                  igseq.v.sbjct.start, igseq.v.sbjct.end,
                                  igseq.v.metadata['identity'], igseq.v.metadata['score']])
        else:
            metadata_list.extend([None, None, None, None, None, None])
        if igseq.d is not None:
            metadata_list.extend([None, None, None, None,
                                 igseq.d.metadata['identity'], igseq.d.metadata['score']])
        else:
            metadata_list.extend([None, None, None, None, None, None])
        if igseq.j is not None:
            metadata_list.extend([igseq.j.query.start, igseq.j.query.end,
                                  igseq.j.sbjct.start, igseq.j.sbjct.end,
                                  igseq.j.metadata['identity'], igseq.j.metadata['score']])
        else:
            metadata_list.extend([None, None, None, None, None, None])
        if igseq.variable_region is not None and igseq.variable_region.untemplate_region is not None:
            metadata_list.extend(igseq.variable_region.untemplate_region)
            metadata_list.extend(['.', '.', '.', '.'])
        else:
            metadata_list.extend([None, None, None, None, None, None])
        if igseq.variable_region is not None and igseq.variable_region.cdr3 is not None:
            metadata_list.extend(igseq.variable_region.cdr3)
            metadata_list.extend(['.', '.', '.', '.'])
        else:
            metadata_list.extend([None, None, None, None, None, None])
        
        metadata_list = self.__utils.none_to_dot(metadata_list, convert_to_str = True, empty_to_dot = True)
        record = tmpl % tuple(metadata_list)
        return record
    
    
    
    
    
    def __iter__(self):
        return self
    
        
    def __next__(self):
        self.parse()
    
    
    def next(self):
        """Return the next IgSeq object from the iterator.
        
        Return the next **IgSeq** object during parsing PYDAIR format
        file with `parse` method.
        """
        
        return self.__next__()
    
        
    def parse(self):
        """Parse PYDAIR format file into an iterator returning IgSeq object.
        
        Parse PYDAIR format file into an iterator returning **IgSeq** object.
        Typical usage is to loop over the records with `for` statement.
        It is expected to use `close` method to close the file handle after parsing.
        
        >>> pydairio = PyDAIRIO('path_to_file', 'r')
        >>> for igseq in pydairio.parse():
        >>>     print(igseq)
        >>> pydairio.close()
        
        """
        
        igseq = self.__parse_pydair()
        return igseq
    
    
    def __parse_pydair(self):
        for buf in self.__fh:
            buf = buf.replace('\n', '')
            buf_record = buf[3:].split('\t')
            if buf[0:6] == '#BEGIN':
                q_name = None
                v_name = None
                d_name = None
                j_name = None
                q_seq  = None
                v_seq  = None
                d_seq  = None
                j_seq  = None
                u_seq  = None
                c_seq  = None
                v_del = None
                j_del = None
                vj_ins = None
                aligned_q_seq = None
                aligned_v_seq = None
                aligned_d_seq = None
                aligned_j_seq = None
                alignv_qstart   = None
                alignv_qend     = None
                alignv_sstart   = None
                alignv_send     = None
                alignv_identity = None
                alignv_socre    = None
                alignd_qstart   = None
                alignd_qend     = None
                alignd_sstart   = None
                alignd_send     = None
                alignd_identity = None
                alignd_socre    = None
                alignj_qstart   = None
                alignj_qend     = None
                alignj_sstart   = None
                alignj_send     = None
                alignj_identity = None
                alignj_socre    = None
                untmpl_start    = None
                untmpl_end      = None
                cdr3_start      = None
                cdr3_end        = None
                q_orf           = None
            if buf[0:2] == 'QN':
                q_name = buf_record[0] if buf_record[0] is not '.' else None
            if buf[0:2] == 'VN':
                v_name = buf_record[0] if buf_record[0] is not '.' else None
            if buf[0:2] == 'DN':
                d_name = buf_record[0] if buf_record[0] is not '.' else None
            if buf[0:2] == 'JN':
                j_name = buf_record[0] if buf_record[0] is not '.' else None
            if buf[0:2] == 'OP':
                q_orf = buf_record[0] if buf_record[0] is not '.' else None
            if buf[0:2] == 'QA':
                aligned_q_seq = buf[3:]
                q_seq = buf[3:].replace(' ', '').replace('-', '')
            if buf[0:2] == 'VA':
                aligned_v_seq = buf[3:]
                v_seq = buf[3:].replace(' ', '').replace('-', '')
            if buf[0:2] == 'JA':
                aligned_j_seq = buf[3:]
                j_seq = buf[3:].replace(' ', '').replace('-', '')
            if buf[0:2] == 'UA':
                u_seq = buf[3:].replace(' ', '').replace('-', '')
            if buf[0:2] == 'CA':
                c_seq = buf[3:].replace(' ', '').replace('-', '')
            if buf[0:2] == 'VD':
                v_del = buf[3:].replace('.', '')
            if buf[0:2] == 'JD':
                j_del = buf[3:].replace('.', '')
            if buf[0:2] == 'VJ':
                vj_ins = buf[3:].replace('.', '')
            if buf[0:2] == 'AL':
                buf_record = buf[6:].split('\t')
                if buf[3:5] == 'QV':
                    alignv_qstart   = int(buf_record[0]) if buf_record[0] is not '.' else None
                    alignv_qend     = int(buf_record[1]) if buf_record[1] is not '.' else None
                    alignv_sstart   = int(buf_record[2]) if buf_record[2] is not '.' else None
                    alignv_send     = int(buf_record[3]) if buf_record[3] is not '.' else None
                    alignv_identity = float(buf_record[4]) if buf_record[4] is not '.' else None
                    alignv_score    = float(buf_record[5]) if buf_record[5] is not '.' else None
                if buf[3:5] == 'QD':
                    alignd_qstart   = int(buf_record[0]) if buf_record[0] is not '.' else None
                    alignd_qend     = int(buf_record[1]) if buf_record[1] is not '.' else None
                    alignd_sstart   = int(buf_record[2]) if buf_record[2] is not '.' else None
                    alignd_send     = int(buf_record[3]) if buf_record[3] is not '.' else None
                    alignd_identity = float(buf_record[4]) if buf_record[4] is not '.' else None
                    alignd_score    = float(buf_record[5]) if buf_record[5] is not '.' else None
                if buf[3:5] == 'QJ':
                    alignj_qstart   = int(buf_record[0]) if buf_record[0] is not '.' else None
                    alignj_qend     = int(buf_record[1]) if buf_record[1] is not '.' else None
                    alignj_sstart   = int(buf_record[2]) if buf_record[2] is not '.' else None
                    alignj_send     = int(buf_record[3]) if buf_record[3] is not '.' else None
                    alignj_identity = float(buf_record[4]) if buf_record[4] is not '.' else None
                    alignj_score    = float(buf_record[5]) if buf_record[5] is not '.' else None
                if buf[3:5] == 'QU':
                    untmpl_start    = int(buf_record[0]) if buf_record[0] is not '.' else None
                    untmpl_end      = int(buf_record[1]) if buf_record[1] is not '.' else None
                if buf[3:5] == 'QC':
                    cdr3_start      = int(buf_record[0]) if buf_record[0] is not '.' else None
                    cdr3_end        = int(buf_record[1]) if buf_record[1] is not '.' else None
                    
            if buf[0:4] == '#END':
                qv_clipped, vv_clipped = self.__clip_alignment(aligned_q_seq, alignv_qstart, alignv_qend,
                                                               aligned_v_seq, alignv_sstart, alignv_send)
                qj_clipped, jj_clipped = self.__clip_alignment(aligned_q_seq, alignj_qstart, alignj_qend,
                                                               aligned_j_seq, alignj_sstart, alignj_send)
                
                igseqquery = IgSeqQuery(q_name, q_seq, '+', q_orf)
                igseqalign_v = IgSeqAlign(IgSeqAlignQuery(q_name, q_seq, qv_clipped, alignv_qstart, alignv_qend, '+'),
                                          IgSeqAlignSbjct(v_name, v_seq, vv_clipped, alignv_sstart, alignv_send, '+'),
                                          alignv_identity, alignv_socre)
                igseqalign_d = IgSeqAlign(IgSeqAlignQuery(q_name, q_seq, qv_clipped, alignd_qstart, alignd_qend, '+'),
                                          IgSeqAlignSbjct(d_name, None,  None,       alignd_sstart, alignd_send, '+'),
                                          alignd_identity, alignd_socre)
                igseqalign_j = IgSeqAlign(IgSeqAlignQuery(q_name, q_seq, qj_clipped, alignj_qstart, alignj_qend, '+'),
                                          IgSeqAlignSbjct(j_name, j_seq, jj_clipped, alignj_sstart, alignj_send, '+'),
                                          alignj_identity, alignj_socre)
                igseqvariable_region = IgSeqVariableRegion(q_name, q_seq, untmpl_start, untmpl_end, cdr3_start, cdr3_end)
                igseqindels = IgSeqIndels(v_del, j_del, vj_ins)
                
                igseq = IgSeq(igseqalign_v, igseqalign_d, igseqalign_j, igseqvariable_region, igseqindels)
                yield igseq
        
        
        
    def __clip_alignment(self, qseq, qstart, qend, sseq, sstart, send):
        qseg = None
        sseg = None
        base_count_q = 0
        base_count_s = 0
        if qstart != None and sstart != None:
            qseg = ''
            sseg = ''
            for base_index in range(max(len(qseq), len(sseq))):
                qbase = qseq[base_index:(base_index + 1)]
                sbase = sseq[base_index:(base_index + 1)]
                if qbase != ' ' and qbase != '-':
                    base_count_q += 1
                if sbase != ' ' and sbase != '-':
                    base_count_s += 1
                if qstart <= base_count_q:
                    qseg += qbase
                if sstart <= base_count_s:
                    sseg += sbase
                if base_count_q == qend and base_count_s == send:
                    break
        return [qseg, sseg]
    
    


