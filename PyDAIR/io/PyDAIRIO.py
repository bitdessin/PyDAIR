import re
import os
import sys
import math
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from PyDAIR.utils.PyDAIRUtils import *
from PyDAIR.seq.IgSeq import *


class PyDAIRIO: 
    """PyDAIR Input/Output interface class.
    
    This class implements the input/output methods to handle PYDAIR format file.
    """
    
    def __init__(self, pydair_file_path, open_mode = None, pydair_file_format = 'pydair'):
        """PyDAIRIO class initialize method.
        
        Args:
            pydair_file_path (str): A file path to open or write.
            open_mode (str): A character (``r``, ``w`` or ``a``) to specify open mode.
            pydair_file_format (str): A file format to open or write.
                                      One of ``pydair`` or ``simple`` should be specified.
        
        Set up file path, open mode, and file format.
        After setting, the method opens the file hanld.
        """
        
        if open_mode is None:
            raise ValueError('Open mode is required. Please set open mode as \'r\' for reading mode, or \'w\' or \'a\' for writting mode.')
        self.__path   = pydair_file_path
        self.__format = pydair_file_format
        self.__mode   = open_mode
        self.__fh     = open(self.__path, self.__mode)
        self.__utils  = PyDAIRUtils()
        
    
    def close(self):
        """Close file handle.
        
        Close file handle that is opened by `__init__`.
        """
        
        self.__fh.close()
    
    
    
    def write(self, igseq):
        """Write ``IgSeq`` objects to a file.
        
        Args:
            igseq (list): A list of IgSeq objects.
        
        Write a list of IgSeq object to a file.
        It is expected to use ``close`` method to close the file handle after writting.
        """
        
        if (not isinstance(igseq, list)) and (not isinstance(igseq, tuple)):
            igseq_list = [igseq]
        else:
            igseq_list = igseq
        
        if self.__format == 'pydair':
            for igseq_single in igseq_list:
                record = self.__convert_to_pydair(igseq_single)
                self.__fh.write(record)
        elif self.__format == 'simple':
            for igseq_single in igseq_list:
                record = self.__convert_to_simplePyDAIR(igseq_single)
                self.__fh.write(record)
        elif self.__format == 'tsv':
            for igseq_single in igseq:
                record = self.__convert_to_tsv(igseq_single)
                self.__fh.write(record)
        else:
            raise ValueError('PyDAIR only accepts TSV format or PyDAIR flat file format. Choose one of \'PyDAIR\', \'tsv\', and \'simple\'')
    
    
    def __convert_to_tsv(self, igseq):
        r = igseq.get_record()
        r = self.__utils.none_to_dot(r, True)
        return '\t'.join(r) + '\n'
    
    
    def __convert_to_pydair(self, igseq):
        alignment = igseq.print_alignment()
        r = igseq.get_record()
        r = self.__utils.none_to_dot(r, True)
        ValueError()
        for a in range(len(r)):
            r[a] = str(r[a])
        
        tmpl = '''#BEGIN
QN %s
VN %s
DN %s
JN %s
OP %s
OC %s
QA %s
VA %s
JA %s
UA %s
CA %s
#CDR3AA %s
#AL QSTART\tQEND\tSSTART\tSEND\tIDENTITY\tSCORE
AL QV %s\t%s\t%s\t%s\t%s\t%s
AL QD %s\t%s\t%s\t%s\t%s\t%s
AL QJ %s\t%s\t%s\t%s\t%s\t%s
AL QU %s\t%s\t%s\t%s\t%s\t%s
AL QC %s\t%s\t%s\t%s\t%s\t%s
#END
'''
        metadata_list = [r[0], r[5], r[8], r[11], r[3], r[4]]
        metadata_list.extend(alignment)
        metadata_list.append(str(Seq(alignment[4].replace('-', '').replace(' ', ''), generic_dna).translate()))
        metadata_list.extend([r[15], r[16], r[18], r[19], r[20], r[21]])  # V
        metadata_list.extend([  '.',   '.',   '.',   '.', r[28], r[29]])  # D
        metadata_list.extend([r[31], r[32], r[34], r[35], r[36], r[37]])  # J
        metadata_list.extend([r[38], r[39],   '.',   '.', '.',   '.'])  # UNALIGNED
        metadata_list.extend([r[40], r[41],   '.',   '.', '.',   '.'])  # CDR3
        record = tmpl % tuple(metadata_list)
        return record
    
    
    def __convert_to_simplePyDAIR(self, igseq):
        alignment = igseq.print_alignment()
        r = igseq.get_record()
        cdr3_data = igseq.get_cdr3_data()
        
        #         qname orf   orfcode vname dname jname cdr3_nucl          CDR3_prot
        record = [r[0], r[3], r[4], r[5], r[8], r[11], cdr3_data.nucl_seq, cdr3_data.prot_seq]
        record = self.__utils.none_to_dot(record, convert_to_str = True)
        return '\t'.join(record) + '\n'
    
    
    
    def __iter__(self):
        return self
    
        
    def __next__(self):
        self.parse()
    
    
    def next(self):
        """Return the next IgSeq object from the iterator.
        
        Return the next IgSeq object during parsing PYDAIR format file with ``parse`` method.
        """
        
        return self.__next__()
    
        
    def parse(self):
        """Parse PYDAIR format file into an iterator returning IgSeq object.
        
        >>> pydairio = PyDAIRIO('path_to_file', 'r', 'pydair')
        >>> for igseq in pydairio.parse():
        >>> print(igseq)
        >>> pydairio.close()
        
        Parse PYDAIR format file into an iterator returning IgSeq object.
        Typical usage is to loop over the records with ``for`` statement.
        It is expected to use ``close`` method to close the file handle after parsing.
        """
        
        igseq = None
        if self.__format == 'pydair':
            igseq = self.__parse_pydair()
        elif self.__format == 'tsv':
            igseq = self.__parse_tsv()
        elif self.__format == 'simple':
            raise ValueError('PyDAIR dose not support parse Simple-PyDAIR(\'simple\' ) format.')
        else:
            raise ValueError('PyDAIR only accepts TSV format or PyDAIR flat file format. Choose one from \'PyDAIR\' and \'tsv\'.')
        return igseq
    
    
    def __parse_tsv(self):
        for buf in self.__fh:
            buf = buf.replace('\n', '')
            if buf != '':
                buf_record = PyDAIRUtils.dot_to_none(buf.split('\t'))
                igseq = IgSeq()
                igseq.set_record(buf_record)
                yield igseq
    
    
    def __parse_pydair(self):
        for buf in self.__fh:
            buf = buf.replace('\n', '')
            buf_record = buf[3:].split('\t')
            if buf[0:6] == '#BEGIN':
                q_name = '.'
                v_name = '.'
                d_name = '.'
                j_name = '.'
                q_seq  = '.'
                v_seq  = '.'
                d_seq  = '.'
                j_seq  = '.'
                u_seq  = '.'
                c_seq  = '.'
                aligned_q_seq = '.'
                aligned_v_seq = '.'
                aligned_d_seq = '.'
                aligned_j_seq = '.'
                alignv_qstart   = '.'
                alignv_qend     = '.'
                alignv_sstart   = '.'
                alignv_send     = '.'
                alignv_identity = '.'
                alignv_socre    = '.'
                alignd_qstart   = '.'
                alignd_qend     = '.'
                alignd_sstart   = '.'
                alignd_send     = '.'
                alignd_identity = '.'
                alignd_socre    = '.'
                alignj_qstart   = '.'
                alignj_qend     = '.'
                alignj_sstart   = '.'
                alignj_send     = '.'
                alignj_identity = '.'
                alignj_socre    = '.'
                untmpl_start    = '.'
                untmpl_end      = '.'
                cdr3_start      = '.'
                cdr3_end        = '.'
                q_orf           = '.'
                q_orfcode       = '.'
            if buf[0:2] == 'QN':
                q_name = buf_record[0]
            if buf[0:2] == 'VN':
                v_name = buf_record[0]
            if buf[0:2] == 'DN':
                d_name = buf_record[0]
            if buf[0:2] == 'JN':
                j_name = buf_record[0]
            if buf[0:2] == 'OP':
                q_orf = buf_record[0]
            if buf[0:2] == 'OC':
                q_orfcode = buf_record[0]
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
            if buf[0:2] == 'AL':
                buf_record = buf[6:].split('\t')
                if buf[3:5] == 'QV':
                    alignv_qstart   = buf_record[0]
                    alignv_qend     = buf_record[1]
                    alignv_sstart   = buf_record[2]
                    alignv_send     = buf_record[3]
                    alignv_identity = buf_record[4]
                    alignv_score    = buf_record[5]
                if buf[3:5] == 'QD':
                    alignd_qstart   = buf_record[0]
                    alignd_qend     = buf_record[1]
                    alignd_sstart   = buf_record[2]
                    alignd_send     = buf_record[3]
                    alignd_identity = buf_record[4]
                    alignd_score    = buf_record[5]
                if buf[3:5] == 'QJ':
                    alignj_qstart   = buf_record[0]
                    alignj_qend     = buf_record[1]
                    alignj_sstart   = buf_record[2]
                    alignj_send     = buf_record[3]
                    alignj_identity = buf_record[4]
                    alignj_score    = buf_record[5]
                if buf[3:5] == 'QU':
                    untmpl_start    = buf_record[0]
                    untmpl_end      = buf_record[1]
                if buf[3:5] == 'QC':
                    cdr3_start      = buf_record[0]
                    cdr3_end        = buf_record[1]
            if buf[0:4] == '#END':
                r = [None] * 42
                # query, v, d, j
                r[0] = q_name if q_name != '.' else None
                r[1] = q_seq  if q_seq  != '.' else None
                r[2] = '+'
                r[5] = v_name if v_name != '.' else None
                r[6] = v_seq  if v_seq  != '.' else None
                r[7] = '+'
                r[8] = d_name if d_name != '.' else None
                r[9] = None
                r[10] = '+'
                r[11] = j_name if j_name != '.' else None
                r[12] = j_seq if j_seq  != '.' else None
                r[13] = '+'
                r[3] = q_orf if q_orf != '.' else None
                r[4] = q_orfcode if q_orfcode != '.' else None
                
                # alignment
                #r[13] = aligned_q_seq      if aligned_q_seq != '.' else None
                r[15] = int(alignv_qstart) if alignv_qstart != '.' else None
                r[16] = int(alignv_qend)   if alignv_qend   != '.' else None
                #r[16] = aligned_v_seq      if aligned_v_seq != '.' else None
                r[18] = int(alignv_sstart) if alignv_sstart != '.' else None
                r[19] = int(alignv_send)   if alignv_send   != '.' else None
                r[20] = float(alignv_identity)  if alignv_identity  != '.' else None
                r[21] = float(alignv_socre)     if alignv_socre     != '.' else None
                
                #r[21] = aligned_q_seq      if aligned_q_seq != '.' else None
                r[23] = int(alignd_qstart) if alignd_qstart != '.' else None
                r[24] = int(alignd_qend)   if alignd_qend   != '.' else None
                #r[24] = aligned_d_seq      if aligned_d_seq != '.' else None
                r[26] = int(alignd_sstart) if alignd_sstart != '.' else None
                r[27] = int(alignd_send)   if alignd_send   != '.' else None
                r[28] = float(alignd_identity)  if alignd_identity  != '.' else None
                r[29] = float(alignd_socre)     if alignd_socre     != '.' else None
    
                #r[29] = aligned_q_seq      if aligned_q_seq != '.' else None
                r[31] = int(alignj_qstart) if alignj_qstart != '.' else None
                r[32] = int(alignj_qend)   if alignj_qend    != '.' else None
                #r[32] = aligned_j_seq      if aligned_j_seq != '.' else None
                r[34] = int(alignj_sstart) if alignj_sstart != '.' else None
                r[35] = int(alignj_send)   if alignj_send   != '.' else None
                r[36] = float(alignj_identity)  if alignj_identity  != '.' else None
                r[37] = float(alignj_socre)     if alignj_socre  != '.' else None
                
                r[38] = int(untmpl_start) if untmpl_start != '.' else None
                r[39] = int(untmpl_end)   if untmpl_end   != '.' else None
                r[40] = int(cdr3_start)   if cdr3_start   != '.' else None
                r[41] = int(cdr3_end)     if cdr3_end     != '.' else None
    
                # clip alignment into segments
                r[14], r[17] = self.__clip_alignment(
                             aligned_q_seq, r[15], r[16],
                             aligned_v_seq, r[18], r[19])
                r[30], r[33] = self.__clip_alignment(
                             aligned_q_seq, r[31], r[32],
                             aligned_j_seq, r[34], r[35])
                if r[14] is not None:
                    r[22] = r[14]
                elif r[30] is not None:
                    r[22] = r[30]
                igseq = IgSeq()
                igseq.set_record(r)
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
    
    


