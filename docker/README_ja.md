# VSCode + Docker で GPU 解析（RTX A5000 / Windows 11）

コンテナで環境を固定して、A5000をそのまま使う手順。**`fair-esm` と `esm`(ESM-3) の衝突を、
イメージを2つに分けることで解決**しています（`zae-esm2` と `zae-esm3`）。

## 0. 前提（Windows 側で一度だけ）
1. **NVIDIA の Windows ドライバ**を最新に（これだけでWSL2/Docker へGPUが渡ります。CUDA Toolkitの別途導入は不要）。
2. **Docker Desktop** をインストールし、Settings で **WSL2 backend** を有効化。
3. **VSCode** ＋ 拡張 **「Dev Containers」**（`ms-vscode-remote.remote-containers`）を入れる。
4. リポジトリを取得：
   ```bash
   git clone https://github.com/TaichiEndoh/zebrafish-aescallii-esm3-alphafold.git
   cd zebrafish-aescallii-esm3-alphafold
   git checkout claude/esm3-alphafold-mapping-s3qn70
   ```

## A. VSCode Dev Containers で開く（日常作業＝ESM-2/ESMFold）
1. VSCodeでこのフォルダを開く → 右下の通知、または `F1` →
   **「Dev Containers: Reopen in Container」** を選ぶ。
2. `.devcontainer/devcontainer.json` に従って `zae-esm2` が自動ビルドされ、
   コンテナ内でリポジトリが `/workspace` にマウントされる。
3. 起動後、ターミナルに `CUDA available: True | NVIDIA RTX A5000` と出れば成功
   （`postCreateCommand` が自動チェック）。
4. あとは通常どおり実行（GPUを使う）：
   ```bash
   python analysis/variant_rediscovery_benchmark.py \
     --fasta data/reference/sequences/positive_controls.fasta \
     --controls "HRAS_HUMAN:G12V,G12D;TP53_HUMAN:R175H" \
     --models esm2_t6_8M_UR50D,esm2_t12_35M_UR50D,esm2_t30_150M_UR50D,esm2_t33_650M_UR50D \
     --scoring wt-marginal --device cuda --out reports/rediscovery.csv
   ```

## B. コマンドライン（compose）で使う
```bash
docker compose -f docker/docker-compose.yml build            # 両イメージをビルド
docker compose -f docker/docker-compose.yml run --rm esm2    # ESM-2 環境のシェル
docker compose -f docker/docker-compose.yml run --rm esm3    # ESM-3 環境のシェル
```
コンテナ内でGPU確認：
```bash
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# -> True NVIDIA RTX A5000
```

## ESM-3 (1.4B) を回すとき（`zae-esm3`）
```bash
docker compose -f docker/docker-compose.yml run --rm esm3
# コンテナ内で：
huggingface-cli login          # HFトークンを入力（EvolutionaryScale/esm3-sm-open-v1 のライセンス同意も必要）
python analysis/esm3_embed.py --pairs reports/esm3_subset_pairs.csv \
   --fasta-a data/sequences/rerio.fasta --fasta-b data/sequences/aesculapii.fasta \
   --out reports/esm3_pair_distances.csv        # ※ embed_esm3() の実装確定後
```
HFキャッシュは名前付きボリューム `hf-cache` に残るので、毎回のダウンロードは不要。

## ESMFold（全長AHR2）について
`esmfold_predict.py --backend api` は短い標的（≤400aa）ならコンテナ内でGPU不要で動く。
**全長AHR2(1027aa)のlocal予測**は openfold 依存が重いので、必要になったら
`Dockerfile.esm2` に openfold を追加するか、transformers 版 ESMFold を使う（別途対応可）。

## データの置き場所
- リポジトリが `/workspace` にマウントされるので、ホスト側の `data/` `reports/` がそのまま見える。
- 全長プロテオームFASTA等の大きな入力は `data/sequences/` に置く（git管理外）。
- 生成した図・表は `reports/` に置いてホスト側でコミット＆push。

## うまくいかないとき
- `torch.cuda.is_available()` が False：NVIDIA Windowsドライバを更新／Docker DesktopのWSL2 backend確認／
  `docker run --rm --gpus all nvidia/cuda:12.1.1-base-ubuntu22.04 nvidia-smi` でGPUパススルー単体確認。
- ビルドが遅い：初回のみイメージDL＋pip。2回目以降はキャッシュされる。
