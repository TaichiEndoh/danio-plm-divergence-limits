import gzip
import math
from pathlib import Path
import pandas as pd
import numpy as np
import torch
import esm

ROOT = Path('/Users/d/AI/zebrafish-aescallii-go-plm')
Z_FASTA = ROOT / 'zebrafish_aescallii_GO_pipeline/data_raw/zebrafish_np_proteins.fasta.gz'
A_FASTA = ROOT / 'zebrafish_aescallii_GO_pipeline/data_raw/Danio_aesculapii_protein.faa.gz'
ORTHO = ROOT / 'zebrafish_aescallii_GO_pipeline/ortholog_table_relio_aesculapii.csv'
DIAMOND = ROOT / 'diamond_results/rerio_vs_aes_raw.tsv'
OUT = ROOT / 'reports/kcnj13_local_plm_summary.tsv'
FIG = ROOT / 'reports/fig_kcnj13_local_plm_bar.png'

TARGETS = {
    'Kcnj13': 'NP_001039014',
    'Cx39.4': 'NP_001038288',
    'Cx41.8': 'NP_001030160',
    'Igsf11': 'NP_001265751',
}


def read_fasta_gz(path):
    seqs = {}
    with gzip.open(path, 'rt') as f:
        key = None
        buf = []
        for line in f:
            if line.startswith('>'):
                if key:
                    seqs[key] = ''.join(buf)
                key = line[1:].strip().split()[0]
                key = key.split('.')[0]
                buf = []
            else:
                buf.append(line.strip())
        if key:
            seqs[key] = ''.join(buf)
    return seqs


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1):
    n, m = len(a), len(b)
    score = np.zeros((n+1, m+1), dtype=int)
    trace = np.zeros((n+1, m+1), dtype=np.int8)  # 0 diag, 1 up, 2 left
    for i in range(1, n+1):
        score[i,0] = score[i-1,0] + gap
        trace[i,0] = 1
    for j in range(1, m+1):
        score[0,j] = score[0,j-1] + gap
        trace[0,j] = 2
    for i in range(1, n+1):
        ai = a[i-1]
        for j in range(1, m+1):
            bj = b[j-1]
            diag = score[i-1,j-1] + (match if ai == bj else mismatch)
            up = score[i-1,j] + gap
            left = score[i,j-1] + gap
            best = max(diag, up, left)
            score[i,j] = best
            if best == diag:
                trace[i,j] = 0
            elif best == up:
                trace[i,j] = 1
            else:
                trace[i,j] = 2
    # backtrack
    i, j = n, m
    a_aln = []
    b_aln = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i,j] == 0:
            a_aln.append(a[i-1])
            b_aln.append(b[j-1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace[i,j] == 1):
            a_aln.append(a[i-1])
            b_aln.append('-')
            i -= 1
        else:
            a_aln.append('-')
            b_aln.append(b[j-1])
            j -= 1
    return ''.join(reversed(a_aln)), ''.join(reversed(b_aln))


def mismatch_positions(a_aln, b_aln):
    apos = -1
    bpos = -1
    diffs = []
    for aa, bb in zip(a_aln, b_aln):
        if aa != '-':
            apos += 1
        if bb != '-':
            bpos += 1
        if aa != '-' and bb != '-' and aa != bb:
            diffs.append((apos, bpos))
    return diffs


def windows_from_positions(seq, positions, w=10):
    windows = []
    seen = set()
    for pos in positions:
        if pos in seen:
            continue
        seen.add(pos)
        start = max(0, pos - w)
        end = min(len(seq), pos + w + 1)
        windows.append(seq[start:end])
    return windows


def embed_seqs(seqs, model, alphabet):
    batch_converter = alphabet.get_batch_converter()
    data = [(str(i), s) for i, s in enumerate(seqs)]
    _, _, tokens = batch_converter(data)
    with torch.no_grad():
        out = model(tokens, repr_layers=[model.num_layers], return_contacts=False)
    reps = out["representations"][model.num_layers]
    embs = []
    for i, s in enumerate(seqs):
        emb = reps[i, 1:len(s)+1].mean(0)
        embs.append(emb)
    return embs


def cosine_distance(a, b):
    return 1 - torch.nn.functional.cosine_similarity(a, b, dim=0).item()


def main():
    zseqs = read_fasta_gz(Z_FASTA)
    aseqs = read_fasta_gz(A_FASTA)

    ortho = pd.read_csv(ORTHO)
    ortho['Danio_rerio_base'] = ortho['Danio_rerio'].str.replace(r'\\.\\d+$','', regex=True)
    ortho['Danio_aesculapii_base'] = ortho['Danio aesculapii'].str.replace(r'\\.\\d+$','', regex=True)

    diamond = pd.read_csv(DIAMOND, sep='\t', names=[
        'qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','qcovhsp','scovhsp'
    ])
    diamond['q_base'] = diamond['qseqid'].str.replace(r'\\.\\d+$','', regex=True)
    diamond['s_base'] = diamond['sseqid'].str.replace(r'\\.\\d+$','', regex=True)

    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()

    rows = []
    for gene, z_id in TARGETS.items():
        z_base = z_id
        a_base = None
        src = None
        match = ortho[ortho['Danio_rerio_base'] == z_base]
        if not match.empty:
            a_base = match.iloc[0]['Danio_aesculapii_base']
            src = 'ortholog_table'
        else:
            dmatch = diamond[diamond['q_base'] == z_base]
            if not dmatch.empty:
                top = dmatch.sort_values('bitscore', ascending=False).iloc[0]
                a_base = top['s_base']
                src = 'diamond_top1'
        z_seq = zseqs.get(z_base)
        a_seq = aseqs.get(a_base) if a_base else None
        if not z_seq or not a_seq:
            rows.append({
                'gene': gene,
                'z_id': z_base,
                'a_id': a_base,
                'source': src,
                'status': 'missing_sequence'
            })
            continue

        # full distance
        emb_full = embed_seqs([z_seq, a_seq], model, alphabet)
        full_dist = cosine_distance(emb_full[0], emb_full[1])

        # local mismatch windows
        a_aln, b_aln = needleman_wunsch(z_seq, a_seq)
        diffs = mismatch_positions(a_aln, b_aln)
        if not diffs:
            rows.append({
                'gene': gene,
                'z_id': z_base,
                'a_id': a_base,
                'source': src,
                'status': 'no_mismatch',
                'n_mismatch': 0,
                'full_plm': full_dist,
                'local_mean_plm': math.nan,
                'local_median_plm': math.nan,
                'windows_used': 0,
            })
            continue
        # cap windows to avoid too many
        max_windows = 50
        z_windows = windows_from_positions(z_seq, [p[0] for p in diffs])[:max_windows]
        a_windows = windows_from_positions(a_seq, [p[1] for p in diffs])[:max_windows]
        # align window counts
        n = min(len(z_windows), len(a_windows))
        z_windows = z_windows[:n]
        a_windows = a_windows[:n]
        if n == 0:
            rows.append({
                'gene': gene,
                'z_id': z_base,
                'a_id': a_base,
                'source': src,
                'status': 'no_windows',
                'n_mismatch': len(diffs),
                'full_plm': full_dist,
                'local_mean_plm': math.nan,
                'local_median_plm': math.nan,
                'windows_used': 0,
            })
            continue
        embs = embed_seqs(z_windows + a_windows, model, alphabet)
        local_dists = []
        for i in range(n):
            local_dists.append(cosine_distance(embs[i], embs[i+n]))
        rows.append({
            'gene': gene,
            'z_id': z_base,
            'a_id': a_base,
            'source': src,
            'status': 'ok',
            'n_mismatch': len(diffs),
            'full_plm': full_dist,
            'local_mean_plm': float(np.mean(local_dists)),
            'local_median_plm': float(np.median(local_dists)),
            'windows_used': n,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep='\t', index=False)

    # figure
    try:
        import matplotlib.pyplot as plt
        plot = out[out['status']=='ok'].copy()
        if not plot.empty:
            x = np.arange(len(plot))
            plt.figure(figsize=(7,4))
            plt.bar(x-0.15, plot['full_plm'], width=0.3, label='Full PLM')
            plt.bar(x+0.15, plot['local_mean_plm'], width=0.3, label='Local mismatch windows')
            plt.xticks(x, plot['gene'], rotation=30, ha='right')
            plt.ylabel('PLM distance')
            plt.title('Kcnj13 vs conserved genes: full vs local PLM')
            plt.legend(frameon=False)
            plt.tight_layout()
            plt.savefig(FIG, dpi=180)
            plt.close()
    except Exception as e:
        print('Plot error', e)


if __name__ == '__main__':
    main()
