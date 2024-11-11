#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import SeqIO
import numpy as np

def parse_fasta_and_substitution_matrix(fasta_path, matrix_path):
    sequences = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]
    sub_matrix = {}
    with open(matrix_path, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    headers = lines[0].split()
    for line in lines[1:]:
        values = line.split()
        row_label = values[0]
        for col_label, score in zip(headers, values[1:]):
            sub_matrix[(row_label, col_label)] = int(score)
    return sequences, sub_matrix

def initialize_matrices(seq1, seq2, mode, gap=0):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    if mode == "global":
        for i in range(1, m + 1):
            score_matrix[i][0] = gap * i
        for j in range(1, n + 1):
            score_matrix[0][j] = gap * j
    return score_matrix

def calculate_global_alignment(seq1, seq2, sub_matrix, gap):
    score_matrix = initialize_matrices(seq1, seq2, "global", gap)
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i-1][j-1] + sub_matrix[(seq1[i-1], seq2[j-1])]
            delete = score_matrix[i-1][j] + gap
            insert = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match, delete, insert)
    return score_matrix

def calculate_local_alignment(seq1, seq2, sub_matrix):
    m, n = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)
    max_score = 0
    max_positions = []

    # Calculate the score matrix and find the max score
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i-1][j-1] + sub_matrix[(seq1[i-1], seq2[j-1])]
            score_matrix[i][j] = max(0, match)  # Ensure no negative scores
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_positions = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_positions.append((i, j))

    return score_matrix, max_positions

def traceback(seq1, seq2, score_matrix, sub_matrix, gap, max_positions=None):
    aligned_seq1, aligned_seq2 = "", ""
    
    if max_positions is None:  # Global alignment
        i, j = len(seq1), len(seq2)
    else:  # Local alignment (start from max score position)
        i, j = max_positions[0]

    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        
        # Local alignment: stop when we reach a score of 0 only if all neighbors are 0
        if max_positions is not None and current_score == 0:
            if score_matrix[i-1][j] == 0 and score_matrix[i][j-1] == 0 and score_matrix[i-1][j-1] == 0:
                break
        
        if current_score == score_matrix[i-1][j-1] + sub_matrix[(seq1[i-1], seq2[j-1])]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + gap:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif current_score == score_matrix[i][j-1] + gap:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            # For local alignment, allow continuation if score is 0
            if max_positions is not None:
                i -= 1
                j -= 1
            else:
                break  # For global alignment, stop on mismatch

    return aligned_seq1, aligned_seq2

def alignment(input_fasta, substitution_matrix, output_fasta, mode, gap):
    seqs, sub_matrix = parse_fasta_and_substitution_matrix(input_fasta, substitution_matrix)
    seq1, seq2 = seqs[0], seqs[1]

    # 使用 fasta 檔案中的 ID 作為輸出標識符
    records = list(SeqIO.parse(input_fasta, "fasta"))
    seq1_id, seq2_id = records[0].id, records[1].id

    if mode == "global":
        score_matrix = calculate_global_alignment(seq1, seq2, sub_matrix, gap)
        aligned_seq1, aligned_seq2 = traceback(seq1, seq2, score_matrix, sub_matrix, gap)
    elif mode == "local":
        score_matrix, max_positions = calculate_local_alignment(seq1, seq2, sub_matrix)
        aligned_seq1, aligned_seq2 = traceback(seq1, seq2, score_matrix, sub_matrix, gap, max_positions)

    # 將對齊結果寫入 output_fasta 文件
    with open(output_fasta, "w") as out:
        out.write(f">{seq1_id}\n{aligned_seq1}\n")
        out.write(f">{seq2_id}\n{aligned_seq2}\n")

